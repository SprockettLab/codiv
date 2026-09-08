# Collapse monophyletic groups in a phylogenetic tree

Identifies groups of symbiont tips that are monophyletic (all descended
from the same host) and reduces each such group to a single
representative tip. Used during codiversification analysis to remove
within-host variation and focus on between-host branching patterns.

## Usage

``` r
collapse_monophyletic(symbiont_subtree, Host_to_Symbiont_df)
```

## Arguments

- symbiont_subtree:

  A phylo object representing a subtree of the symbiont phylogeny

- Host_to_Symbiont_df:

  A data frame with "Host" and "Symbiont" columns linking symbionts to
  their source hosts

## Value

A phylo object in which every monophyletic group of tips from a single
host has been reduced to one representative tip. Tips from a host that
do not form a single monophyletic group (polyphyly) are kept, one per
monophyletic group. When a group carries branch lengths the tip furthest
from the group's root is kept; otherwise the first tip is kept.

## Details

Collapsing does not depend on branch lengths. A monophyletic group of
tips from one host is reduced to a single tip whether its members sit on
long branches, on branches of equal length, on zero-length branches
(identical sequences), or on a tree with no branch lengths at all.

## Examples

``` r
host_tree <- ape::rtree(5)
symbiont_tree <- ape::rtree(30)
hs_df <- data.frame(
  Host = rep(host_tree$tip.label, 6),
  Symbiont = symbiont_tree$tip.label
)
collapsed <- collapse_monophyletic(symbiont_tree, hs_df)
```
