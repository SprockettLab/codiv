# Prepare symbiont tree with unique node labels

Ensures symbiont tree has unique, non-null node labels for tracking
nodes through the codiversification scan.

## Usage

``` r
.prepare_symbiont_tree(tree, verbose = TRUE)
```

## Arguments

- tree:

  A phylo object to prepare

- verbose:

  If TRUE, print status messages

## Value

A phylo object with guaranteed unique node labels
