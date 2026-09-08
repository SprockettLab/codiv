
<!-- README.md is generated from README.Rmd. Please edit that file -->

# codiv <img src="man/figures/logo.png" align="right" height="139" alt="codiv hex logo" />

Host-Microbe Codiversification Scans in R

<!-- badges: start -->

[![R-CMD-check](https://github.com/SprockettLab/codiv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SprockettLab/codiv/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

------------------------------------------------------------------------

`codiv` scans pairs of phylogenetic trees — a host tree and a symbiont
tree — and asks, at every internal node of the symbiont tree, whether
the symbionts’ branching pattern mirrors their hosts’. Nodes where
symbiont topology tracks host topology are evidence of
**codiversification**: hosts and their microbes diversifying together.

The main function, `codiv()`, takes a host tree, a symbiont tree, and a
data frame linking each symbiont to the host it was isolated from. For
each qualifying node it runs one or more statistical methods (Hommola’s
test, PACo, ParaFit, and an opt-in topology-only test) and reports a
test statistic and a per-node permutation p-value. Helper functions
cover scan-wide false-discovery control, leave-one-host-out analysis,
parameter sensitivity sweeps, host-level summaries, molecular-clock
regression, tanglegram plotting, and ground-truth simulation.

## Installation

Install the development version from
[GitHub](https://github.com/SprockettLab/codiv):

``` r
# install.packages("pak")
pak::pak("SprockettLab/codiv")
```

## Quick example

``` r
library(codiv)

# simulate a small dataset with known co-diversification
sim <- simulate_codiv_data(n_hosts = 12, n_clades = 20, seed = 1)

# scan every qualifying node of the symbiont tree
results <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
                 methods = "hommola", permutations = 99)

results            # compact overview: settings + significant nodes per method
summary(results)   # per-method statistic ranges and significant-node counts
```

The package also ships the rodent gut-microbiome dataset from Sprockett
et al. (2025) in `inst/extdata`, so you can run a real scan immediately
— see `vignette("codiv")`.

## Citation

> Daniel D. Sprockett, Brian A. Dillard, Abigail A. Landers, Jon G.
> Sanders, Andrew H. Moeller. (2025) **Recent genetic drift in the
> co-diversified gut bacterial symbionts of laboratory mice.** *Nature
> Communications.* <https://doi.org/10.1038/s41467-025-57435-z>
