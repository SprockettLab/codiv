# Simulate host and symbiont trees with known co-diversification

Builds a host tree and a symbiont tree assembled from several clades,
each with known ground truth, all mapping to the same set of hosts. Use
it for a fixed sanity-check or to generate many random datasets to study
how tree features affect the scan.

## Usage

``` r
simulate_codiv_data(
  n_hosts = 12,
  clades = NULL,
  n_clades = NULL,
  prop_codiversifying = 0.5,
  congruence_range = c(0.5, 1),
  n_per_host_range = c(1, 4),
  bl_noise = 0.1,
  host_tree = NULL,
  seed = NULL
)
```

## Arguments

- n_hosts:

  Number of host tips (ignored if `host_tree` is supplied)

- clades:

  Optional list of clade specifications, each a list with
  `codiversifying` (logical) and optionally `congruence` and
  `n_per_host`. If NULL, `n_clades` random clades are generated.

- n_clades:

  Number of random clades to generate when `clades` is NULL

- prop_codiversifying:

  Probability a random clade is co-diversifying

- congruence_range:

  Range to draw random co-diversifying `congruence` from

- n_per_host_range:

  Integer range to draw random `n_per_host` from

- bl_noise:

  Multiplicative branch-length jitter (0 = none)

- host_tree:

  Optional host tree (phylo); if NULL a random coalescent tree of
  `n_hosts` tips is used

- seed:

  Optional random seed

## Value

A list with:

- `host_tree`: the host tree (phylo)

- `symbiont_tree`: the combined symbiont tree (phylo)

- `links`: data frame of Host-Symbiont associations

- `truth`: per-symbiont ground truth (clade, codiversifying, congruence,
  n_per_host)

## Details

Each clade is one monophyletic block of the symbiont tree, controlled
by:

- `codiversifying`: whether its topology tracks the host tree

- `congruence`: for co-diversifying clades, how closely it tracks the
  host (1 = identical branching order; lower = more scrambled)

- `n_per_host`: symbionts sampled per host (1 = no within-host
  pseudoreplication; \> 1 introduces it)

## Examples

``` r
# fixed four-clade sanity check
sim <- simulate_codiv_data(
  n_hosts = 10,
  clades = list(
    list(codiversifying = TRUE,  congruence = 1.0, n_per_host = 1),
    list(codiversifying = TRUE,  congruence = 0.6, n_per_host = 4),
    list(codiversifying = FALSE, n_per_host = 1),
    list(codiversifying = FALSE, n_per_host = 4)
  ),
  seed = 1
)
sim$host_tree
#> 
#> Phylogenetic tree with 10 tips and 9 internal nodes.
#> 
#> Tip labels:
#>   H1, H2, H3, H4, H5, H6, ...
#> 
#> Rooted; includes branch length(s).
head(sim$links)
#>   Host Symbiont
#> 1   H1    c1_h1
#> 2   H2    c1_h2
#> 3   H3    c1_h3
#> 4   H4    c1_h4
#> 5   H5    c1_h5
#> 6   H6    c1_h6

# many random datasets (vary seed) to study tree-feature effects
rand <- simulate_codiv_data(n_hosts = 15, n_clades = 6, seed = 42)
```
