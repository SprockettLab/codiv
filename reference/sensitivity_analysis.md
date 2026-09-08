# Parameter Sensitivity Analysis

Reruns
[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
across a grid of `span_fraction` and `min_symbiont_tips` values and
summarizes how the number of scanned and significant nodes changes, so
you can see whether findings are stable or threshold-dependent.

## Usage

``` r
sensitivity_analysis(
  Host_tree,
  Symbiont_tree,
  Host_to_Symbiont_df,
  span_fraction_range = c(0.05, 0.1, 0.2, 0.3),
  min_symbiont_tips_range = c(5, 7, 10, 15),
  permutations = 99,
  min_hosts = 3,
  methods = c("hommola", "paco", "parafit"),
  alpha = 0.05,
  cores = 1,
  seed = 8675309,
  verbose = TRUE
)
```

## Arguments

- Host_tree:

  A phylo object; binary phylogenetic tree of hosts

- Symbiont_tree:

  A phylo object; binary phylogenetic tree of symbionts

- Host_to_Symbiont_df:

  A data frame with columns "Host" and "Symbiont"

- span_fraction_range:

  Numeric vector of `span_fraction` values to test; default c(0.05, 0.1,
  0.2, 0.3)

- min_symbiont_tips_range:

  Integer vector of `min_symbiont_tips` values; default c(5, 7, 10, 15)

- permutations:

  Integer; number of permutations, applied to every run

- min_hosts:

  Integer; minimum hosts, applied to every run (\>= 3)

- methods:

  Character vector of methods to use

- alpha:

  Significance threshold for counting significant nodes; default 0.05

- cores:

  Number of cores passed to
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)

- seed:

  Random seed for reproducibility

- verbose:

  If TRUE (default), report progress across the grid

## Value

A list with elements:

- `sensitivity_data`: one row per parameter combination, with the number
  of nodes scanned and, per method, the number of significant nodes and
  the mean of the primary statistic

- `results_by_params`: named list of the full
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  result for each parameter combination

- `heatmap_data`: long-format data (one row per combination and method)
  with `n_sig`, `n_nodes`, and `prop_sig`, ready for plotting

## Details

Runtime is roughly
`baseline_time x length(span_fraction_range) x length(min_symbiont_tips_range)`.
Tree-distance metrics are turned off for speed, since sensitivity is
judged from node counts and significance.

## Examples

``` r
# \donttest{
sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
sens <- sensitivity_analysis(sim$host_tree, sim$symbiont_tree, sim$links,
                             span_fraction_range = c(0.3, 0.6),
                             min_symbiont_tips_range = c(7, 10),
                             permutations = 99, methods = "hommola",
                             verbose = FALSE)
sens$sensitivity_data
#>   span_fraction min_symbiont_tips n_nodes hommola_n_sig hommola_mean_stat
#> 1           0.3                 7      10             0       -0.05226943
#> 2           0.6                 7      12             0       -0.03900951
#> 3           0.3                10       8             0        0.01857405
#> 4           0.6                10      10             0        0.02031726
# }
```
