# Summarize a codiv result

Per-method summary of a codiversification scan: the range of the primary
statistic and the number of nodes with a per-node permutation p-value
below `alpha`. This p-value count is descriptive and not corrected
across nodes; use
[`codiv_null_scans()`](https://www.sprockettlab.com/codiv/reference/codiv_null_scans.md)
for a scan-level false-discovery assessment.

## Usage

``` r
# S3 method for class 'codiv'
summary(object, alpha = 0.05, ...)
```

## Arguments

- object:

  A `codiv` object from
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md).

- alpha:

  Significance threshold; default 0.05.

- ...:

  Ignored.

## Value

A data frame with one row per method (invisibly returned by print),
containing the statistic range and the per-node significant count.

## Examples

``` r
# \donttest{
sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
             methods = c("hommola", "paco"), permutations = 99)
#> Creating unique node labels
#> Of the 59 internal nodes in the symbiont tree,
#>    49 (83%) have span > 0 and <= 10% of max (0 dropped for zero span)
#>    5 (8%) have 7-500 symbiont tips (0 dropped as too large)
#>    5 (8%) have >= 3 hosts
#> Scanning 5 nodes for codiversification.
#> Scanning 5 nodes across 3 cores ...
summary(res)
#> codiv summary: 5 nodes, per-node p < 0.05 
#> 
#>   method statistic   stat_min stat_median  stat_max n_sig_p n_nodes
#>  hommola Hommola_r -0.4922685  -0.2492913 0.5514694       0       5
#>     paco   PACo_ss  0.8180839   6.2805921 9.0588433       0       5
summary(res, alpha = 0.01)   # stricter significance threshold
#> codiv summary: 5 nodes, per-node p < 0.01 
#> 
#>   method statistic   stat_min stat_median  stat_max n_sig_p n_nodes
#>  hommola Hommola_r -0.4922685  -0.2492913 0.5514694       0       5
#>     paco   PACo_ss  0.8180839   6.2805921 9.0588433       0       5
# }
```
