# Convert host-symbiont associations to binary presence-absence matrix

Transforms a long-format data frame of host-symbiont associations into a
wide-format binary matrix suitable for co-phylogenetic analysis methods
like PACo and ParaFit.

## Usage

``` r
host_symbiont_links(Host_to_Symbiont_df)
```

## Arguments

- Host_to_Symbiont_df:

  A data frame with columns "Host" and "Symbiont" where each row
  represents one association. Symbionts must be unique (one symbiont per
  host per row).

## Value

A binary matrix with hosts as rows and symbionts as columns, filled with
0s and 1s where 1 indicates the symbiont was isolated from that host.

## Examples

``` r
hs_df <- data.frame(
  Host = c("H1", "H1", "H2", "H2", "H3"),
  Symbiont = c("S1", "S2", "S3", "S4", "S5")
)
mat <- host_symbiont_links(hs_df)
```
