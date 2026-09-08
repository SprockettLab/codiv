# Calculate Hommola correlation coefficient

Computes Pearson correlation between pairwise phylogenetic distances in
host and symbiont trees. Tests the null hypothesis of no association
between host and symbiont topology (Hommola et al. 2011).

## Usage

``` r
hommola(i_host_subtree_dist, i_symbiont_subtree_dist, i_host_to_symbiont_df)
```

## Arguments

- i_host_subtree_dist:

  A distance matrix (e.g., from
  [`dist()`](https://rdrr.io/r/stats/dist.html) or
  [`adephylo::distTips()`](https://rdrr.io/pkg/adephylo/man/distTips.html))
  for the host subtree tips

- i_symbiont_subtree_dist:

  A distance matrix for the symbiont subtree tips

- i_host_to_symbiont_df:

  A data frame with "Host" and "Symbiont" columns linking each symbiont
  to its host

## Value

A numeric correlation coefficient (Pearson's r) ranging from -1 to 1.
Higher positive values indicate symbiont topology more closely matches
host topology. Returns 0 if all host distances are identical (degenerate
case).

## Examples

``` r
# \donttest{
library(adephylo)
#> Loading required package: ade4
h_tree <- ape::rtree(5)
s_tree <- ape::rtree(15)
hs_df <- data.frame(
  Host = rep(h_tree$tip.label, 3),
  Symbiont = s_tree$tip.label
)
h_dist <- distTips(h_tree, method = "patristic")
s_dist <- distTips(s_tree, method = "patristic")
r <- hommola(h_dist, s_dist, hs_df)
# }
```
