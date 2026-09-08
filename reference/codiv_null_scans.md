# Second-order permutation test for scan-wide co-diversification

Assesses whether a co-diversification scan detects more significant
clades than expected by chance, accounting for the phylogenetic
non-independence and pseudoreplication of nodes that per-node FDR
corrections (e.g. Benjamini- Hochberg) ignore. The host tree's tip
labels are permuted and the entire
[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md) scan
is rerun; repeating this builds a null distribution for the *number* of
significant clades. Because each permutation reruns the whole pipeline,
the null captures exactly the dependence structure that a per-node
correction does not.

## Usage

``` r
codiv_null_scans(
  Host_tree,
  Symbiont_tree,
  Host_to_Symbiont_df,
  n_permutations = 100,
  statistic = NULL,
  stat_threshold = 0.75,
  pvalue = NULL,
  p_threshold = 0.01,
  observed = NULL,
  seed = 8675309,
  verbose = TRUE,
  ...
)
```

## Arguments

- Host_tree:

  A phylo object; binary host tree

- Symbiont_tree:

  A phylo object; binary symbiont tree

- Host_to_Symbiont_df:

  A data frame with "Host" and "Symbiont" columns

- n_permutations:

  Number of host-label permutations (outer null); the analyses in the
  source papers used 100

- statistic:

  Name of the test-statistic column used to call a node significant.
  NULL (default) uses `Hommola_r`, the collapsed Hommola correlation
  that
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  always writes when the Hommola method runs

- stat_threshold:

  A node counts as significant when its statistic exceeds this value;
  default 0.75

- pvalue:

  Name of the p-value column used to call a node significant; NULL
  (default) resolves alongside `statistic`

- p_threshold:

  A node counts as significant when its p-value is below this value;
  default 0.01

- observed:

  Optional precomputed
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  result for the real (unpermuted) data; if NULL it is computed

- seed:

  Random seed for reproducibility

- verbose:

  If TRUE (default), report progress across permutations

- ...:

  Additional arguments passed to
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  (e.g. min_hosts, min_symbiont_tips, span_fraction, permutations,
  methods, cores)

## Value

A list with:

- `observed`: number of significant clades in the real data

- `null_counts`: number of significant clades in each permuted scan

- `global_pvalue`: (permutations with count \>= observed + 1) / (n + 1)

- `empirical_fdr`: expected false-discovery proportion at the chosen
  thresholds, `mean(null_counts) / observed` (capped at 1)

- `thresholds`: the statistic/p-value columns and cutoffs used

## Details

This is the recommended way to gauge false discoveries for a whole scan
(Sanders et al. 2023; Sprockett et al. 2025).
[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
deliberately does not apply a per-node multiple-testing correction,
which would be invalid for non-independent nested clades.

Runtime is roughly `(n_permutations + 1) x` a single scan, so use
`cores` and consider fewer inner `permutations` for a first pass. Each
rerun uses its own checkpoint file and does not resume.

## Examples

``` r
# \donttest{
sim <- simulate_codiv_data(n_hosts = 8, n_clades = 2, seed = 1)
null <- codiv_null_scans(sim$host_tree, sim$symbiont_tree, sim$links,
                         n_permutations = 10, min_hosts = 3,
                         min_symbiont_tips = 7, span_fraction = 0.6,
                         permutations = 49, methods = "hommola",
                         verbose = FALSE)
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
#> Warning: Recommend setting `permutations` to 99 or more for statistical validity
null$observed
#> [1] 0
null$empirical_fdr
#> [1] NA
# }
```
