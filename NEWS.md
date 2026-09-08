# codiv 0.1.0

* Initial release.
* `codiv()` scans each qualifying internal node of a symbiont tree for
  codiversification with a host tree using the Hommola, PACo, and ParaFit
  methods, plus an opt-in branch-length-free topology-congruence test.
* `codiv_null_scans()` provides a scan-level (second-order permutation)
  false-discovery assessment that respects the non-independence of nested
  clades.
* Secondary analyses: `loo_host_analysis()`, `sensitivity_analysis()`,
  `host_codiv_summary()`, `molecular_clock()`, `filter_results()`, and
  `plot_codiv_trees()` tanglegrams.
* `simulate_codiv_data()` generates host/symbiont datasets with known
  ground truth.
* Ships the rodent gut-microbiome dataset from Sprockett et al. (2025) in
  `inst/extdata`.
