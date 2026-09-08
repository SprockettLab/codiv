# Changelog

## codiv 0.1.0

- Initial release.
- [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
  scans each qualifying internal node of a symbiont tree for
  codiversification with a host tree using the Hommola, PACo, and
  ParaFit methods, plus an opt-in branch-length-free topology-congruence
  test.
- [`codiv_null_scans()`](https://www.sprockettlab.com/codiv/reference/codiv_null_scans.md)
  provides a scan-level (second-order permutation) false-discovery
  assessment that respects the non-independence of nested clades.
- Secondary analyses:
  [`loo_host_analysis()`](https://www.sprockettlab.com/codiv/reference/loo_host_analysis.md),
  [`sensitivity_analysis()`](https://www.sprockettlab.com/codiv/reference/sensitivity_analysis.md),
  [`host_codiv_summary()`](https://www.sprockettlab.com/codiv/reference/host_codiv_summary.md),
  [`molecular_clock()`](https://www.sprockettlab.com/codiv/reference/molecular_clock.md),
  [`filter_results()`](https://www.sprockettlab.com/codiv/reference/filter_results.md),
  and
  [`plot_codiv_trees()`](https://www.sprockettlab.com/codiv/reference/plot_codiv_trees.md)
  tanglegrams.
- [`simulate_codiv_data()`](https://www.sprockettlab.com/codiv/reference/simulate_codiv_data.md)
  generates host/symbiont datasets with known ground truth.
- Ships the rodent gut-microbiome dataset from Sprockett et al. (2025)
  in `inst/extdata`.
