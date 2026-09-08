# Phylogenetic Association with Co-diversification (PACo) analysis

Performs PACo (Phylogenetic Association with Co-diversification)
analysis to test host-symbiont co-evolution. Uses principal coordinates
to reduce phylogenetic distance matrices and computes a goodness-of-fit
statistic.

## Usage

``` r
paco_wf(
  i_host_subtree,
  i_symbiont_subtree,
  i_host_to_symbiont_mat,
  permutations,
  seed
)
```

## Arguments

- i_host_subtree:

  A phylo object for the host subtree

- i_symbiont_subtree:

  A phylo object for the symbiont subtree

- i_host_to_symbiont_mat:

  A binary matrix with hosts as rows and symbionts as columns, where 1
  indicates an association

- permutations:

  Number of randomizations for significance testing; at least 99
  recommended

- seed:

  Random seed for reproducibility; if NA, no seed is set

## Value

A PACo object containing:

- gof: Goodness-of-fit statistics including ss (squared sum) and p-value
  from permutation test

- Additional fields from paco::PACo()

## References

Balbuena JA, Míguez-Argüello R, Blasco-Costa I, Llopis-Beltrán A (2013).
PACo: A novel procedure to estimate the packedness of a given
association matrix. Journal of Biogeography 40: 948-960.

## Examples

``` r
# \donttest{
h_tree <- ape::rtree(5)
s_tree <- ape::rtree(15)
links <- data.frame(Host = rep(h_tree$tip.label, 3),
                    Symbiont = s_tree$tip.label)
hs_mat <- host_symbiont_links(links)   # binary host-by-symbiont matrix
result <- paco_wf(h_tree, s_tree, hs_mat, permutations = 99, seed = 123)
# }
```
