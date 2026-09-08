# Package index

## Main Functions

Core codiversification analysis

- [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md) :
  Perform codiversification scan between host and symbiont trees
- [`check_inputs()`](https://www.sprockettlab.com/codiv/reference/check_inputs.md)
  : Validate inputs to codiv()
- [`codiv_null_scans()`](https://www.sprockettlab.com/codiv/reference/codiv_null_scans.md)
  : Second-order permutation test for scan-wide co-diversification
- [`host_codiv_summary()`](https://www.sprockettlab.com/codiv/reference/host_codiv_summary.md)
  : Summarize which hosts have co-diversifying symbionts

## Working with Results

Methods for the codiv result object

- [`print(`*`<codiv>`*`)`](https://www.sprockettlab.com/codiv/reference/print.codiv.md)
  : Print a codiv result
- [`summary(`*`<codiv>`*`)`](https://www.sprockettlab.com/codiv/reference/summary.codiv.md)
  : Summarize a codiv result
- [`as.data.frame(`*`<codiv>`*`)`](https://www.sprockettlab.com/codiv/reference/as.data.frame.codiv.md)
  : Coerce a codiv result to a plain data frame
- [`filter_results()`](https://www.sprockettlab.com/codiv/reference/filter_results.md)
  : Filter a codiversification scan, optionally collapsing nested calls

## Codiversification Methods

Statistical methods used at each node

- [`hommola()`](https://www.sprockettlab.com/codiv/reference/hommola.md)
  : Calculate Hommola correlation coefficient
- [`hommola_wf()`](https://www.sprockettlab.com/codiv/reference/hommola_wf.md)
  : Hommola correlation with permutation significance testing
- [`paco_wf()`](https://www.sprockettlab.com/codiv/reference/paco_wf.md)
  : Phylogenetic Association with Co-diversification (PACo) analysis
- [`topology_wf()`](https://www.sprockettlab.com/codiv/reference/topology_wf.md)
  : Topology-congruence test with permutation significance

## Helper Functions

Tree and data manipulation utilities

- [`collapse_monophyletic()`](https://www.sprockettlab.com/codiv/reference/collapse_monophyletic.md)
  : Collapse monophyletic groups in a phylogenetic tree
- [`host_symbiont_links()`](https://www.sprockettlab.com/codiv/reference/host_symbiont_links.md)
  : Convert host-symbiont associations to binary presence-absence matrix

## Secondary Analyses

Robustness, interpretation, and visualization

- [`loo_host_analysis()`](https://www.sprockettlab.com/codiv/reference/loo_host_analysis.md)
  : Leave-One-Host-Out Analysis
- [`sensitivity_analysis()`](https://www.sprockettlab.com/codiv/reference/sensitivity_analysis.md)
  : Parameter Sensitivity Analysis
- [`molecular_clock()`](https://www.sprockettlab.com/codiv/reference/molecular_clock.md)
  : Molecular-clock corroboration of co-diversification
- [`plot_codiv_trees()`](https://www.sprockettlab.com/codiv/reference/plot_codiv_trees.md)
  : Visualize Host-Symbiont Codiversification

## Simulation

Generate host/symbiont data with known co-diversification

- [`simulate_codiv_data()`](https://www.sprockettlab.com/codiv/reference/simulate_codiv_data.md)
  : Simulate host and symbiont trees with known co-diversification
