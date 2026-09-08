# Identify nodes to scan and apply filtering criteria

Builds list of symbiont tree subtrees and filters by span, tips, and
hosts.

## Usage

``` r
.identify_nodes_to_scan(
  Symbiont_tree,
  Host_to_Symbiont_df,
  min_hosts,
  min_symbiont_tips,
  span_fraction,
  max_symbiont_tips = NULL,
  verbose = TRUE
)
```

## Arguments

- Symbiont_tree:

  A phylo object with unique node labels

- Host_to_Symbiont_df:

  Host-symbiont linkage data frame

- min_hosts:

  Minimum hosts threshold

- min_symbiont_tips:

  Minimum symbiont tips threshold

- span_fraction:

  Span cutoff as fraction of max span

- verbose:

  If TRUE, print filtering statistics

## Value

A list with elements:

- nodes_to_scan: Character vector of node labels to analyze

- subtree_list: Named list of subtrees

- Symbiont_df: Data frame with filtering metadata
