# Leave-One-Host-Out Analysis

Reruns
[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md) with
each host removed in turn and compares each run to the full analysis,
showing how much each host contributes to the codiversification signal.

## Usage

``` r
loo_host_analysis(
  Host_tree,
  Symbiont_tree,
  Host_to_Symbiont_df,
  hosts_to_test = NULL,
  alpha = 0.05,
  verbose = TRUE,
  ...
)
```

## Arguments

- Host_tree:

  A phylo object; binary phylogenetic tree of hosts

- Symbiont_tree:

  A phylo object; binary phylogenetic tree of symbionts

- Host_to_Symbiont_df:

  A data frame with columns "Host" and "Symbiont" linking symbionts to
  hosts

- hosts_to_test:

  Character vector of host names to test. If NULL (default), tests all
  hosts. Use a subset to reduce computation time.

- alpha:

  Significance threshold for counting significant nodes; default 0.05

- verbose:

  If TRUE (default), report progress across hosts

- ...:

  Additional arguments passed to
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  (e.g. min_hosts, min_symbiont_tips, span_fraction, permutations,
  methods, cores, seed)

## Value

A list with elements:

- `full_results`: the
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  result on the complete data

- `loo_results`: named list of the
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  result with each host removed

- `host_importance`: one row per tested host (ordered by
  `impact_score`), with the number of nodes scanned, per-method
  significant-node counts and their change from the full analysis
  (`<method>_delta_n_sig`), and an `impact_score` (mean absolute change
  in significant-node count across methods when the host is removed)

## Details

Each run uses its own isolated checkpoint file so removing a host
genuinely re-computes the affected nodes. Tree-distance metrics are
disabled for speed. Runtime is roughly
`baseline_time x length(hosts_to_test)`, so consider `hosts_to_test`,
fewer permutations, or a single method for a first pass, and `cores` to
parallelize each run.

## Examples

``` r
# \donttest{
sim <- simulate_codiv_data(n_hosts = 8, n_clades = 2, seed = 1)
loo <- loo_host_analysis(sim$host_tree, sim$symbiont_tree, sim$links,
                         permutations = 99, methods = "hommola",
                         verbose = FALSE)
loo$host_importance
#>   host n_nodes hommola_n_sig hommola_delta_n_sig hommola_mean_stat impact_score
#> 1   H8       1             0                   0        -0.5240938            0
#> 2   H1       1             0                   0        -0.4636440            0
#> 3   H6       2             0                   0         0.5264145            0
#> 4   H4       1             0                   0         0.9987138            0
#> 5   H3       1             0                   0        -0.4120341            0
#> 6   H7       2             0                   0         0.5264145            0
#> 7   H5       2             0                   0         0.5264145            0
#> 8   H2       2             0                   0         0.5264145            0
# }
```
