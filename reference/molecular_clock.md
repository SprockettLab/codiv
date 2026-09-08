# Molecular-clock corroboration of co-diversification

For the most strongly co-diversifying clades, regresses symbiont
divergence against host clade age. A positive relationship corroborates
that symbiont clades diversified contemporaneously with their hosts, and
the slope is an estimate of the symbiont molecular rate (substitutions
per unit host time). This follows the approach of Sanders et al. (2023)
and Sprockett et al. (2025).

## Usage

``` r
molecular_clock(
  codiv_results,
  host_time_tree,
  statistic = NULL,
  stat_threshold = 0.95,
  pvalue = NULL,
  p_threshold = 0.01,
  symbiont_divergence = NULL
)
```

## Arguments

- codiv_results:

  A `codiv` result (or data frame) from
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)

- host_time_tree:

  A time-calibrated host tree (phylo, ultrametric) whose branch lengths
  are in time units (e.g. millions of years)

- statistic:

  Name of the test-statistic column. NULL (default) uses `Hommola_r`,
  the collapsed Hommola correlation that
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  always writes when the Hommola method runs

- stat_threshold:

  Only clades with a statistic above this value are used; default 0.95
  (strong co-diversification)

- pvalue:

  Name of the p-value column. NULL (default) resolves alongside
  `statistic`

- p_threshold:

  Only clades with a p-value below this value are used; default 0.01

- symbiont_divergence:

  Optional override for the symbiont-divergence summary. Either a named
  numeric vector keyed by `Node_ID`, or a function taking a symbiont
  subtree (phylo) and returning a single number. If NULL (default), the
  clade's crown depth is used (mean root-to-tip distance in the symbiont
  tree's branch-length units)

## Value

A list with:

- `model`: the fitted `lm` of symbiont divergence on host age

- `slope`, `intercept`, `r_squared`, `p_value`: regression summaries

- `slope_ci`: 95% confidence interval for the slope

- `per_clade`: data frame with one row per clade used (`Node_ID`,
  `n_hosts`, `host_age`, `symbiont_divergence`)

## Details

Host clade age is the crown age of the clade's hosts (the age of their
most recent common ancestor in `host_time_tree`). The default symbiont
divergence is a simplification: it is the clade crown depth in
substitution units, not a dedicated absolute date. For rigorous absolute
dating supply your own estimates via `symbiont_divergence`. Nested
clades are non-independent; the regression treats each passing node as a
point, so consider restricting to non-nested clades for formal
inference.

## Examples

``` r
# \donttest{
sim <- simulate_codiv_data(n_hosts = 14, clades = list(
  list(codiversifying = TRUE, congruence = 1, n_per_host = 3)), seed = 1)
res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links, span_fraction = 0.9,
             methods = "hommola", permutations = 199, verbose = FALSE)
#> Warning: One or both trees are not fully bifurcating; resolve polytomies with ape::multi2di() for correct tree-shape metrics
# sim$host_tree is ultrametric, so it stands in for a time-calibrated tree
mc <- molecular_clock(res, sim$host_tree, stat_threshold = 0.75)
mc$slope
#> [1] 1.010627
mc$per_clade
#>   Node_ID n_hosts  host_age symbiont_divergence
#> 1 Node_29       5 0.1622165           0.1642994
#> 2  Node_4       6 0.4981807           0.5164868
#> 3  Node_3       9 0.7521906           0.7612640
#> 4  Node_2      14 1.9897942           2.0152759
# }
```
