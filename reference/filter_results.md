# Filter a codiversification scan, optionally collapsing nested calls

Subsets the output of
[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md) on
any combination of its columns, on which hosts a clade contains, and on
nesting. Scanned nodes are nested by construction, so a single
codiversifying lineage is usually reported many times over: once for the
whole clade and once for most of its subclades. Counting those rows
directly overstates how many distinct lineages were found. The `nested`
argument reduces each nested series to one call.

## Usage

``` r
filter_results(
  codiv_results,
  ...,
  hosts = NULL,
  hosts_match = c("all", "any"),
  nested = c("keep_all", "outermost", "innermost", "best"),
  statistic = "Hommola_r",
  sort_by_filters = FALSE
)
```

## Arguments

- codiv_results:

  A `codiv` result (or data frame) from
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)

- ...:

  Filtering expressions evaluated inside `codiv_results`, in the style
  of [`subset()`](https://rdrr.io/r/base/subset.html), for example
  `Hommola_pvalue < 0.05` and `Hommola_r > 0.7`. Multiple expressions
  are combined with AND. Rows where an expression is NA are dropped,
  matching [`subset()`](https://rdrr.io/r/base/subset.html).

- hosts:

  Optional character vector of host names to require.
  [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  writes a `<host>_PRESENT` column for every name passed to its
  `focus_hosts` argument, and those columns are what this reads, so a
  host can only be filtered on if it was tracked during the scan.

- hosts_match:

  Whether a clade must contain all of `hosts` ("all", the default) or at
  least one of them ("any")

- nested:

  How to treat nested calls, applied after every other filter:

  - `"keep_all"` (default) leaves them alone, so the result still counts
    one lineage many times

  - `"outermost"` keeps the largest call in each nested series

  - `"innermost"` keeps the calls that contain no other surviving call.
    A call can hold two disjoint subclades, so this can return more rows
    than `"outermost"`, not fewer

  - `"best"` keeps whichever call in each series has the highest
    `statistic`. The strongest signal often sits inside a series rather
    than at its outermost call, so this frequently differs from
    `"outermost"`

- statistic:

  Column used to choose a winner when `nested = "best"`; default
  "Hommola_r"

- sort_by_filters:

  If TRUE, order the surviving rows by the columns used in `...`, in the
  order those conditions were written, taking the direction from each
  comparison: `Hommola_pvalue < 0.01` sorts that column ascending,
  `Hommola_r > 0.5` sorts it descending. Only simple
  `column op constant` comparisons give a direction, so compound
  conditions and `==` are ignored for sorting. Default FALSE, which
  leaves the scan's own order alone.

## Value

A `codiv` object holding the surviving rows, carrying the same
`codiv_params` attribute as the input plus a `codiv_filter` attribute
recording what was applied.

## Details

Filtering always runs before de-nesting. The other order would pick a
representative for a series and then possibly discard it, silently
losing the whole series even when other members passed the filter.

De-nesting reads the per-node symbiont trees stored by
[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md), so
`Symbiont_Tree` must still be present.

## See also

[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md),
[`host_codiv_summary()`](https://www.sprockettlab.com/codiv/reference/host_codiv_summary.md)

## Examples

``` r
# \donttest{
sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)
res <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
             methods = "hommola", permutations = 99, verbose = FALSE)

# every significant call, including nested duplicates
filter_results(res, Hommola_pvalue < 0.05)
#> <codiv> codiversification scan
#>   Nodes scanned: 0 
#>   Methods: hommola 
#>   Permutations: 99 
#>   Nodes with p < 0.05: hommola 0/0
#> 
#> (no nodes met the filtering thresholds)

# one row per distinct lineage
filter_results(res, Hommola_pvalue < 0.05, nested = "outermost")
#> <codiv> codiversification scan
#>   Nodes scanned: 0 
#>   Methods: hommola 
#>   Permutations: 99 
#>   Nodes with p < 0.05: hommola 0/0
#> 
#> (no nodes met the filtering thresholds)
# }
```
