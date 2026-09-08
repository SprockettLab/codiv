# Print a codiv result

Compact overview of a codiversification scan: the run settings, how many
nodes are significant per method, and a short preview of the strongest
nodes (Newick and tree-distance columns are hidden).

## Usage

``` r
# S3 method for class 'codiv'
print(x, alpha = 0.05, n = 6, ...)
```

## Arguments

- x:

  A `codiv` object from
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md).

- alpha:

  Significance threshold for the summary counts; default 0.05.

- n:

  Number of preview rows to show; default 6.

- ...:

  Ignored.

## Value

`x`, invisibly.

## Examples

``` r
# \donttest{
sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
             methods = "hommola", permutations = 99)
#> Creating unique node labels
#> Of the 59 internal nodes in the symbiont tree,
#>    49 (83%) have span > 0 and <= 10% of max (0 dropped for zero span)
#>    5 (8%) have 7-500 symbiont tips (0 dropped as too large)
#>    5 (8%) have >= 3 hosts
#> Scanning 5 nodes for codiversification.
#> Scanning 5 nodes across 3 cores ...
res                 # calls print.codiv
#> <codiv> codiversification scan
#>   Nodes scanned: 5 
#>   Methods: hommola 
#>   Permutations: 99 
#>   Nodes with p < 0.05: hommola 0/5
#> 
#> Top nodes:
#>  Node_ID N_Symbionts N_Hosts  Hommola_r Hommola_pvalue
#>  Node_49          12       4  0.5514694           0.18
#>  Node_36           9       3 -0.4922685           0.67
#>  Node_15          12       6 -0.2673787           0.72
#>  Node_17          10       5 -0.2492913           0.80
#>  Node_18           8       4 -0.1790181           0.73
print(res, n = 10)  # show more preview rows
#> <codiv> codiversification scan
#>   Nodes scanned: 5 
#>   Methods: hommola 
#>   Permutations: 99 
#>   Nodes with p < 0.05: hommola 0/5
#> 
#> Top nodes:
#>  Node_ID N_Symbionts N_Hosts  Hommola_r Hommola_pvalue
#>  Node_49          12       4  0.5514694           0.18
#>  Node_36           9       3 -0.4922685           0.67
#>  Node_15          12       6 -0.2673787           0.72
#>  Node_17          10       5 -0.2492913           0.80
#>  Node_18           8       4 -0.1790181           0.73
# }
```
