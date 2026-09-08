## Submission

This is the first submission of codiv (0.1.0), an R package for host-microbe
codiversification analysis. It implements the microbiota-wide co-diversification
scan described in Sanders et al. (2023) and Sprockett et al. (2025)
<doi:10.1038/s41467-025-57435-z>.

## Test environments

- local: macOS 26 (aarch64), R 4.5.1 -- 0 errors | 0 warnings | 1 note
- win-builder (devel and release) -- <fill in before submitting>
- macOS builder (mac-builder.r-project.org) -- <fill in before submitting>
- R-hub: linux (R-devel), windows (R-devel), macos (R-devel) -- <fill in>

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission (this is the first CRAN release of codiv).

## Notes

- All imported packages (adephylo, ape, castor, paco, parallel, pbmcapply,
  scales, TreeDist) are on CRAN.
- Examples that run a full permutation scan are wrapped in \donttest{} to keep
  each individual example under a few seconds; they are still executed by
  `R CMD check --run-donttest` and complete without error.
- The vignette runs a real (small) codiversification scan on the bundled
  example data in a few seconds.
- The scan is parallelised with pbmcapply::pbmclapply(); it honours
  `_R_CHECK_LIMIT_CORES_` and never uses more than 2 cores during checks.
- ggplot2 is used only by plot_codiv_trees(), which is guarded with
  requireNamespace().
