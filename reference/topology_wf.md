# Topology-congruence test with permutation significance

Tests whether a symbiont tree's branching order matches a host tree's,
ignoring branch lengths. Each host is reduced to one representative tip
(its largest single-host clade, tie-broken by smallest total clade
branch length); hosts whose symbionts are entirely scattered have no
defined position and are excluded. Both trees are then pruned to the
shared hosts and compared with normalized mutual clustering information
(1 = identical branching order). The null distribution is built by
shuffling the host labels.

## Usage

``` r
topology_wf(
  host_tree,
  symbiont_subtree,
  host_to_symbiont_df,
  permutations,
  seed
)
```

## Arguments

- host_tree:

  A phylo object; the host (sub)tree

- symbiont_subtree:

  A phylo object; the symbiont (sub)tree

- host_to_symbiont_df:

  A data frame with "Host" and "Symbiont" columns

- permutations:

  Number of permutations for the null; at least 99 recommended

- seed:

  Random seed for reproducibility; if NA, no seed is set

## Value

A list with:

- `stat`: normalized mutual clustering information (0-1; 1 = identical
  branching order), or NA when fewer than 3 hosts have a defined
  position

- `pvalue`: one-tailed permutation p-value

- `hosts_used`: number of hosts with a defined position

- `fraction_kept`: fraction of symbiont tips represented after reducing
  each host to its dominant clade (low values flag heavy pruning)

## Details

This complements the distance-based methods (Hommola, PACo, ParaFit): it
is insensitive to branch-length differences, so it can detect congruent
topology even when evolutionary rates differ between the trees.

## Examples

``` r
# \donttest{
h <- ape::rcoal(8)
s <- h; s$tip.label <- paste0("s", seq_along(s$tip.label))
hs <- data.frame(Host = h$tip.label, Symbiont = s$tip.label)
topology_wf(h, s, hs, permutations = 99, seed = 1)
#> $stat
#> [1] 1
#> 
#> $pvalue
#> [1] 0.01
#> 
#> $hosts_used
#> [1] 8
#> 
#> $fraction_kept
#> [1] 1
#> 
# }
```
