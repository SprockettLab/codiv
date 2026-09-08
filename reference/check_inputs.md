# Validate inputs to codiv()

Checks that phylogenetic trees are valid, parameters are in acceptable
ranges, and the host-to-symbiont linking data frame has the correct
structure.

## Usage

``` r
check_inputs(
  Host_tree,
  Symbiont_tree,
  Host_to_Symbiont_df,
  min_hosts,
  min_symbiont_tips,
  span_fraction,
  permutations
)
```

## Arguments

- Host_tree:

  A phylo object; binary phylogenetic tree of hosts

- Symbiont_tree:

  A phylo object; rooted, binary phylogenetic tree of symbionts. Rooting
  is required: see Details.

- Host_to_Symbiont_df:

  A data frame with columns "Host" and "Symbiont" that links each
  symbiont tip to the host it was isolated from

- min_hosts:

  Minimum number of hosts required for a node to be included in the
  scan; must be \>= 3

- min_symbiont_tips:

  Minimum number of symbiont tips in a subtree; default 7 is recommended

- span_fraction:

  Fraction of total tree span; filters nodes with span \<= this fraction
  of maximum span (0-1)

- permutations:

  Number of permutations for significance testing; default 99
  recommended, higher values more accurate

## Value

Invisibly returns TRUE if all checks pass; stops with error message if
validation fails

## Details

`Symbiont_tree` must be rooted.
[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md) tests
one clade per internal node, so the rooting decides which clades exist
to be tested. An unrooted tree read by
[`ape::read.tree()`](https://rdrr.io/pkg/ape/man/read.tree.html) still
has a node sitting in the root position, picked by the order of the
Newick string, so a scan would run and report an arbitrary set of clades
rather than failing. Root the tree before resolving polytomies:
[`ape::multi2di()`](https://rdrr.io/pkg/ape/man/multi2di.html) applied
to an unrooted tree splits the basal polytomy with a zero-length branch,
after which [`ape::is.rooted()`](https://rdrr.io/pkg/ape/man/root.html)
returns TRUE for a tree nobody rooted.

## Examples

``` r
Host_tree <- ape::rtree(10)
Symbiont_tree <- ape::rtree(100)
Host_to_Symbiont_df <- data.frame(
  "Host" = rep(Host_tree$tip.label, 10),
  "Symbiont" = Symbiont_tree$tip.label
)
check_inputs(Host_tree, Symbiont_tree, Host_to_Symbiont_df, 4, 7, 0.25, 99)
```
