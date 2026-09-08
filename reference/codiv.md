# Perform codiversification scan between host and symbiont trees

Scans a symbiont phylogenetic tree for evidence of codiversification
with a host tree. For each internal node of the symbiont tree meeting
size thresholds, calculates correlation coefficients and p-values using
one or more methods (Hommola, PACo, ParaFit) to test if the symbiont
topology matches the host topology at that node.

## Usage

``` r
codiv(
  Host_tree,
  Symbiont_tree,
  Host_to_Symbiont_df,
  min_hosts = 3,
  min_symbiont_tips = 7,
  max_symbiont_tips = NULL,
  span_fraction = 0.1,
  uncollapsed_hommola = FALSE,
  permutations = 99,
  seed = 8675309,
  verbose = TRUE,
  Save_fp = NULL,
  methods = c("hommola", "paco", "parafit"),
  focus_hosts = NULL,
  subtree_features = TRUE,
  cores = NULL,
  batch_size = NULL,
  n_batches = NULL,
  continue = TRUE
)
```

## Arguments

- Host_tree:

  A phylo object; binary phylogenetic tree of hosts

- Symbiont_tree:

  A phylo object; binary phylogenetic tree of symbionts

- Host_to_Symbiont_df:

  A data frame with columns "Host" and "Symbiont" linking each symbiont
  tip to its source host

- min_hosts:

  Minimum number of hosts required at a node for inclusion in the scan;
  default 3

- min_symbiont_tips:

  Minimum number of symbiont tips in a subtree for inclusion; default 7
  (higher values recommended for statistical power)

- max_symbiont_tips:

  Maximum number of symbiont tips in a subtree. If NULL (default),
  `max(500, 2% of the tree's tips)`. This is not the same control as
  `span_fraction`: that limits how deep a clade is, this limits how many
  genomes it holds, and the two are only loosely related. A densely
  sampled recent radiation is enormous but shallow, so it passes any
  reasonable span cutoff. Very large clades are also where the
  permutation test reaches significance at trivial effect sizes simply
  because n is large.

- span_fraction:

  Filters nodes to only those with phylogenetic span \<= this fraction
  of total tree span (0-1); default 0.1 focuses on distal nodes

- uncollapsed_hommola:

  If TRUE, also run the Hommola test on the uncollapsed subtree and
  report `Uncollapsed_Hommola_r`/`Uncollapsed_Hommola_pvalue`. Default
  FALSE. The uncollapsed test counts every MAG from a host that formed a
  monophyletic group as an independent observation, so a deeply sampled
  host is pseudoreplicated; PACo and ParaFit already use the collapsed
  tree. It is also by far the most expensive statistic on large clades.

- permutations:

  Number of randomizations for significance testing; default 99; higher
  values (999+) recommended for publication

- seed:

  Random seed for reproducible permutations; default 8675309

- verbose:

  If TRUE (default), prints progress and filtering statistics

- Save_fp:

  File path where results are saved after each node (enables resumable
  runs via `continue`); if NULL (default), a unique temporary file is
  used, so separate calls never share a checkpoint. Supply an explicit
  path to enable resuming a specific run

- methods:

  Character vector of methods to use: any combination of "hommola",
  "paco", "parafit" (default uses these three, all distance-based), and
  the opt-in "topology" (a branching-order congruence test that ignores
  branch lengths; each host is reduced to its largest single-host clade
  and hosts with no coherent position are excluded, reported via the
  Topology_hosts_used and Topology_fraction_kept columns)

- focus_hosts:

  Optional character vector of host names. For each name a logical
  `<host>_PRESENT` column is added to the results, TRUE at nodes whose
  subtree contains that host. Useful for sorting/filtering the scan when
  you care about co-diversifying symbionts in particular hosts. Names
  not found in `Host_to_Symbiont_df` trigger a warning and are FALSE
  everywhere.

- subtree_features:

  If TRUE (default), calculates multiple tree distance metrics
  (TreeDistance, SharedPhylogeneticInfo, etc.) for each node

- cores:

  Number of cores for parallelizing the node scan. If NULL (default),
  uses all detected cores minus one. Forking is unavailable on Windows,
  where the scan always runs on a single core.

- batch_size:

  Number of nodes to scan per checkpoint batch. Results are written to
  `Save_fp` after each batch, so smaller values checkpoint more often
  (more restart safety, slightly more disk I/O). If NULL (default), a
  size of `max(200, cores * 50)` is used. Mutually exclusive with
  `n_batches`.

- n_batches:

  Alternative to `batch_size`: split the scanned nodes into this many
  roughly equal checkpoint batches. Mutually exclusive with
  `batch_size`.

- continue:

  If TRUE (default), resumes an interrupted scan by skipping completed
  nodes; requires Save_fp file to exist

## Value

A `codiv` object: a data frame (with one row per symbiont tree node
scanned) carrying the run parameters as a `codiv_params` attribute, with
`print` and `summary` methods. It behaves like a data frame for
subsetting and column access. Columns include:

- Node_ID: Node identifier

- N_Symbionts, N_Hosts: Tip counts in each subtree

- Symbiont_Colless, Symbiont_Sackin: Tree shape metrics

- Hommola_r, Hommola_pvalue: Hommola Pearson correlation on the
  collapsed subtree and its permutation p-value (if the hommola method
  is used)

- Uncollapsed_Hommola_r, Uncollapsed_Hommola_pvalue: the same test on
  the uncollapsed subtree (only when `uncollapsed_hommola = TRUE`)

- PACo_ss, PACo_pvalue: PACo results (if method used)

- ParaFitGlobal, ParaFit_pvalue: ParaFit results (if method used)

- Topology_congruence, Topology_pvalue, Topology_hosts_used,
  Topology_fraction_kept: topology-congruence results (if method used)

- Tree distance metrics if subtree_features=TRUE

P-values are the per-node permutation p-values. They are not corrected
for multiple testing across nodes: nested clades are non-independent, so
a per-node correction (e.g. Benjamini-Hochberg) would be invalid. Use
[`codiv_null_scans()`](https://www.sprockettlab.com/codiv/reference/codiv_null_scans.md)
for a scan-level false-discovery assessment.

## Examples

``` r
# \donttest{
# simulate a small host/symbiont dataset with known co-diversification
sim <- simulate_codiv_data(n_hosts = 10, n_clades = 3, seed = 1)

results <- codiv(sim$host_tree, sim$symbiont_tree, sim$links,
                 methods = "hommola", permutations = 99, verbose = FALSE)
results            # compact overview
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
head(as.data.frame(results))
#>         Node_ID
#> Node_18 Node_18
#> Node_36 Node_36
#> Node_17 Node_17
#> Node_49 Node_49
#> Node_15 Node_15
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Symbiont_Tree
#> Node_18                                                                                                                                                                                        ((c2_h3_1:0.004661357195,c2_h3:0.004661357195)Node_24:0.01087650012,((c2_h10_1:0.001186443053,c2_h10:0.001186443053)Node_23:0.002768367124,((c2_h9_1:0.002584480569,c2_h9:0.002584480569)Node_22:0.006030454661,(c2_h2_1:0.002765316812,c2_h2:0.002765316812)Node_21:0.006452405895)Node_20:0.004002682002)Node_19:0.004072273792)Node_18;
#> Node_36                                                                                                                                                              ((c3_h4_1:0.07625152206,(c3_h4_2:0.02287545662,c3_h4:0.02287545662)Node_43:0.05337606544)Node_42:0.1779202181,((c3_h3_1:0.05618300743,(c3_h3_2:0.01685490223,c3_h3:0.01685490223)Node_41:0.0393281052)Node_40:0.131093684,(c3_h2_1:0.0565443785,(c3_h2_2:0.01696331355,c3_h2:0.01696331355)Node_39:0.03958106495)Node_38:0.1319368832)Node_37:0.05085388256)Node_36;
#> Node_17                                                                                                ((c2_h4_1:0.01835908762,c2_h4:0.01835908762)Node_25:0.04283787112,((c2_h3_1:0.004661357195,c2_h3:0.004661357195)Node_24:0.01087650012,((c2_h10_1:0.001186443053,c2_h10:0.001186443053)Node_23:0.002768367124,((c2_h9_1:0.002584480569,c2_h9:0.002584480569)Node_22:0.006030454661,(c2_h2_1:0.002765316812,c2_h2:0.002765316812)Node_21:0.006452405895)Node_20:0.004002682002)Node_19:0.004072273792)Node_18:0.04606688055)Node_17;
#> Node_49 (((c3_h10_1:0.007972939947,(c3_h10_2:0.002391881984,c3_h10:0.002391881984)Node_59:0.005581057963)Node_58:0.01860352654,(c3_h9_1:0.008837663511,(c3_h9_2:0.002651299053,c3_h9:0.002651299053)Node_57:0.006186364458)Node_56:0.02062121486)Node_55:0.08991485234,((c3_h8_1:0.004766214497,(c3_h8_2:0.001429864349,c3_h8:0.001429864349)Node_54:0.003336350148)Node_53:0.01112116716,(c3_h7_1:0.004478456821,(c3_h7_2:0.001343537046,c3_h7:0.001343537046)Node_52:0.003134919774)Node_51:0.01044973258)Node_50:0.1188468404)Node_49;
#> Node_15        (((c2_h4_1:0.01835908762,c2_h4:0.01835908762)Node_25:0.04283787112,((c2_h3_1:0.004661357195,c2_h3:0.004661357195)Node_24:0.01087650012,((c2_h10_1:0.001186443053,c2_h10:0.001186443053)Node_23:0.002768367124,((c2_h9_1:0.002584480569,c2_h9:0.002584480569)Node_22:0.006030454661,(c2_h2_1:0.002765316812,c2_h2:0.002765316812)Node_21:0.006452405895)Node_20:0.004002682002)Node_19:0.004072273792)Node_18:0.04606688055)Node_17:0.02968134993,(c2_h1_1:0.02884987917,c2_h1:0.02884987917)Node_16:0.06731638474)Node_15;
#>         N_Symbionts Symbiont_Colless Symbiont_Sackin
#> Node_18           8                6              26
#> Node_36           9                6              30
#> Node_17          10               12              38
#> Node_49          12                4              44
#> Node_15          12               20              52
#>                                                                                                                                                                      Host_Tree
#> Node_18                                                                    (((H2:0.01678181851,H3:0.01678181851):0.04468413943,H4:0.06146595794):1.659956878,H10:1.721422836);
#> Node_36                                                                                                     ((H1:0.09053719966,H4:0.09053719966):1.630885636,H10:1.721422836);
#> Node_17                                  (((H2:0.01678181851,H3:0.01678181851):0.04468413943,(H4:0.05480904071,H5:0.05480904071):0.006656917232):1.659956878,H10:1.721422836);
#> Node_49                                                                       (((H2:0.06146595794,H5:0.06146595794):0.5234951044,H6:0.5849610623):0.17989428,H9:0.7648553423);
#> Node_15 ((H1:0.09053719966,((H2:0.01678181851,H3:0.01678181851):0.04468413943,(H4:0.05480904071,H5:0.05480904071):0.006656917232):0.02907124172):1.630885636,H10:1.721422836);
#>         N_Hosts Host_Colless Host_Sackin  Hommola_r Hommola_pvalue TreeDistance
#> Node_18       4            3           9 -0.1790181           0.73    1.0000000
#> Node_36       3            1           5 -0.4922685           0.67          NaN
#> Node_17       5            3          13 -0.2492913           0.80    0.7734456
#> Node_49       4            3           9  0.5514694           0.18    0.0000000
#> Node_15       6            7          19 -0.2673787           0.72    0.8199910
#>         SharedPhylogeneticInfo DifferentPhylogeneticInfo NyeSimilarity
#> Node_18              0.0000000                  3.169925     0.3333333
#> Node_36              0.0000000                  0.000000     0.0000000
#> Node_17              0.7369656                  7.813781     1.0000000
#> Node_49              1.5849625                  0.000000     1.0000000
#> Node_15              0.9708537                 15.639388     1.2500000
#>         JaccardRobinsonFoulds MatchingSplitDistance MatchingSplitInfoDistance
#> Node_18              1.333333                     2                  3.169925
#> Node_36              0.000000                     0                  0.000000
#> Node_17              2.000000                     3                  6.117787
#> Node_49              0.000000                     0                  0.000000
#> Node_15              3.500000                     7                 11.241245
#>         MutualClusteringInfo
#> Node_18            0.0000000
#> Node_36            0.0000000
#> Node_17            0.4399462
#> Node_49            1.0000000
#> Node_15            0.5032583

# Bundled example trees (real data) load with system.file():
# host <- ape::read.tree(system.file("extdata", "Host_Tree.nwk",
#                                    package = "codiv"))
# }
```
