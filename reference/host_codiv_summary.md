# Summarize which hosts have co-diversifying symbionts

Reports, for each host, how many co-diversifying clades it participates
in. A node counts as co-diversifying when its statistic exceeds
`stat_threshold` and its p-value is below `p_threshold` (the defaults,
`Hommola_r > 0.75` and `p < 0.01`, follow Sanders et al. 2023 and
Sprockett et al. 2025). The hosts in each node are read from the
per-node host tree stored by
[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md).

## Usage

``` r
host_codiv_summary(
  codiv_results,
  statistic = NULL,
  stat_threshold = 0.75,
  pvalue = NULL,
  p_threshold = 0.01
)
```

## Arguments

- codiv_results:

  A `codiv` result (or data frame) from
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)

- statistic:

  Name of the test-statistic column. NULL (default) uses `Hommola_r`,
  the collapsed Hommola correlation that
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  always writes when the Hommola method runs

- stat_threshold:

  A node is co-diversifying when its statistic exceeds this value;
  default 0.75

- pvalue:

  Name of the p-value column. NULL (default) resolves alongside
  `statistic`

- p_threshold:

  A node is co-diversifying when its p-value is below this value;
  default 0.01

## Value

A data frame with one row per host (ordered by `n_codiversifying_nodes`,
descending):

- `Host`: host tip label

- `n_codiversifying_nodes`: co-diversifying nodes that include this host

- `n_scanned_nodes`: all scanned nodes that include this host

- `prop_codiversifying`: `n_codiversifying_nodes / n_scanned_nodes`

- `max_statistic`: the largest statistic value among this host's
  co-diversifying nodes (NA if none)

## Details

Scanned nodes are nested, so a host appears in several overlapping
nodes; the counts are a descriptive ranking of which hosts carry the
most co-diversifying symbiont lineages, not independent tallies.

## Examples

``` r
# \donttest{
sim <- simulate_codiv_data(n_hosts = 10, clades = list(
  list(codiversifying = TRUE,  congruence = 1, n_per_host = 2),
  list(codiversifying = FALSE, n_per_host = 2)), seed = 3)
results <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
                 methods = "hommola", permutations = 99, verbose = FALSE)
host_codiv_summary(results)
#>   Host n_codiversifying_nodes n_scanned_nodes prop_codiversifying max_statistic
#> 1   H3                      0               1                   0            NA
#> 2   H4                      0               1                   0            NA
#> 3   H5                      0               2                   0            NA
#> 4   H6                      0               3                   0            NA
#> 5   H7                      0               3                   0            NA
#> 6   H8                      0               3                   0            NA
#> 7   H9                      0               2                   0            NA
# }
```
