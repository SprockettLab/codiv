# Getting Started with codiv

## codiv: Host-Microbe Codiversification Analysis

Welcome to **codiv**! This package performs codiversification scans
between pairs of phylogenetic trees to identify which branches show
evidence of coevolution between hosts and their microbial symbionts.

------------------------------------------------------------------------

### Quick Start

#### In 30 seconds

The [`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
function takes three inputs and returns codiversification statistics:

``` r

library(codiv)

# Inputs:
# 1. Host phylogenetic tree (phylo object)
# 2. Symbiont phylogenetic tree (phylo object)
# 3. Host-to-symbiont linking table (data frame)

codiv_results <- codiv(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  permutations = 99
)

# Output: one row per scanned symbiont node, with statistics and p-values
head(codiv_results)
```

#### What you get

[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md)
returns a `codiv` object, which is a data frame with one row per scanned
symbiont node plus `print` and `summary` methods for a quick overview:

``` r

codiv_results          # compact overview: settings and significant nodes per method
summary(codiv_results) # per-method statistic ranges and significant-node counts
```

The columns include: - **Node identifiers** from the symbiont tree -
**Per-method test statistics** (Hommola’s correlation, PACo
goodness-of-fit, ParaFit global statistic) measuring how well symbiont
topology matches host topology - **Per-node permutation p-values** -
**Tree shape metrics** (Colless, Sackin indices)

Nodes with **high correlation + low p-value** are evidence of
codiversification at that branch. To assess false discoveries for a
whole scan, use
[`codiv_null_scans()`](https://www.sprockettlab.com/codiv/reference/codiv_null_scans.md).

------------------------------------------------------------------------

### Preparing Your Data

Before running
[`codiv()`](https://www.sprockettlab.com/codiv/reference/codiv.md), you
need to prepare three pieces:

1.  **Host phylogenetic tree**
2.  **Symbiont phylogenetic tree**
3.  **Metadata linking hosts to symbionts**

#### 1. Phylogenetic Trees

Trees must be in R’s `phylo` format (from the `ape` package). Examples:

``` r

library(ape)

# Option A: Read from Newick/Nexus files
Host_Tree <- read.tree("host_phylogeny.nwk")
Symbiont_Tree <- read.tree("symbiont_phylogeny.treefile")

# Option B: Read from other formats
library(phangorn)
Host_Tree <- read.nexus("host_phylogeny.nex")

# Option C: Create from distance matrix
library(phangorn)
dist_matrix <- read.csv("distances.csv", row.names = 1)
Host_Tree <- upgma(as.dist(dist_matrix))
```

**Requirements:** - Both trees must be binary (fully resolved) - Tip
labels must be unique within each tree - Node labels are optional (codiv
will create them if missing)

#### 2. Host-to-Symbiont Linking Table

Create a data frame with exactly two columns: “Host” and “Symbiont”

``` r

# Example structure
Host_to_Symbiont <- data.frame(
  Host = c("Mus_musculus", "Rattus_norvegicus", "Peromyscus_leucopus",
           "Mus_musculus", "Rattus_norvegicus"),
  Symbiont = c("MAG_0001__Bacteroides", "MAG_0002__Bacteroides",
               "MAG_0003__Bacteroides", "MAG_0004__Lactobacillus",
               "MAG_0005__Lactobacillus")
)

# Each row = one host-symbiont association
# Hosts are tips from the host tree
# Symbionts are tips from the symbiont tree
# For MAG data: each symbiont appears in only ONE host

head(Host_to_Symbiont)
#>                  Host                Symbiont
#> 1        Mus_musculus   MAG_0001__Bacteroides
#> 2   Rattus_norvegicus   MAG_0002__Bacteroides
#> 3 Peromyscus_leucopus   MAG_0003__Bacteroides
#> 4        Mus_musculus MAG_0004__Lactobacillus
#> 5   Rattus_norvegicus MAG_0005__Lactobacillus
```

**Tips:** - For MAG/isolate data: each symbiont appears in exactly one
host (one-to-one mapping) - Column names MUST be exactly “Host” and
“Symbiont” - No missing values allowed - Save as: CSV, TSV, Excel, or R
object

#### 3. Example Data Included with codiv

You don’t need your own data to try codiv. The package ships example
host and symbiont trees and a linking table in `inst/extdata` — the
rodent gut-microbiome codiversification dataset from [Sprockett et
al. (2025)](https://www.nature.com/articles/s41467-025-57435-z), 21
rodent host species and several thousand bacterial MAGs — loaded with
[`system.file()`](https://rdrr.io/r/base/system.file.html):

``` r

library(codiv)
library(ape)

Host_Tree <- read.tree(
  system.file("extdata", "Host_Tree.nwk", package = "codiv"))
Symbiont_Tree <- read.tree(
  system.file("extdata", "Symbiont_Tree.treefile", package = "codiv"))
Host_to_Symbiont <- read.csv(
  system.file("extdata", "Host_to_Symbiont.txt", package = "codiv"),
  sep = "\t")

# The bundled symbiont tree is unrooted, as most tree inference programs
# leave it. codiv() requires a rooted symbiont tree, because the scan tests
# one clade per internal node and the rooting decides which clades exist.
# Choose the rooting deliberately for real analyses (an outgroup is best);
# midpoint rooting is a reasonable default here.
Symbiont_Tree <- castor::root_at_midpoint(Symbiont_Tree)

# Inspect the data
Host_Tree
#> 
#> Phylogenetic tree with 21 tips and 20 internal nodes.
#> 
#> Tip labels:
#>   Peromyscus_eremicus, Peromyscus_californicus, Peromyscus_polionotus, Peromyscus_maniculatus, Peromyscus_sonoriensis, Peromyscus_leucopus, ...
#> Node labels:
#>   , '24', '26', '28', '8', '13', ...
#> 
#> Rooted; includes branch length(s).
Symbiont_Tree
#> 
#> Phylogenetic tree with 5567 tips and 5566 internal nodes.
#> 
#> Tip labels:
#>   g__CAG-475___Peromyscus_californicus___4921, g__UBA3818___Peromyscus_californicus___1749, f__Ruminococcaceae___Hydrochoerus_hydrochaeris___6040, f__Ruminococcaceae___Rattus_norvegicus___6091, f__Ruminococcaceae___Microtus_ochrogaster___6055, s__Anaerotruncus_colihominis_A___Mus_musculus_domesticus___816, ...
#> Node labels:
#>   , 100, 100, 100, 100, 63, ...
#> 
#> Rooted; includes branch length(s).
head(Host_to_Symbiont)
#>                        Host
#> 1   Peromyscus_californicus
#> 2   Peromyscus_californicus
#> 3 Hydrochoerus_hydrochaeris
#> 4         Rattus_norvegicus
#> 5      Microtus_ochrogaster
#> 6   Mus_musculus_domesticus
#>                                                         Symbiont
#> 1                    g__CAG-475___Peromyscus_californicus___4921
#> 2                    g__UBA3818___Peromyscus_californicus___1749
#> 3          f__Ruminococcaceae___Hydrochoerus_hydrochaeris___6040
#> 4                  f__Ruminococcaceae___Rattus_norvegicus___6091
#> 5               f__Ruminococcaceae___Microtus_ochrogaster___6055
#> 6 s__Anaerotruncus_colihominis_A___Mus_musculus_domesticus___816
```

With 21 hosts and several thousand symbionts, this is a realistic scan.
Running just the Hommola method and skipping the tree-distance metrics
keeps it to a few seconds:

``` r

codiv_results <- codiv(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  methods = "hommola",
  permutations = 99,
  subtree_features = FALSE,
  verbose = FALSE
)

codiv_results
#> <codiv> codiversification scan
#>   Nodes scanned: 801 
#>   Methods: hommola 
#>   Permutations: 99 
#>   Nodes with p < 0.05: hommola 222/801
#> 
#> Top nodes:
#>   Node_ID N_Symbionts N_Hosts Hommola_r Hommola_pvalue
#>  2268_100           7       3 0.9999826           0.37
#>  5173_100           8       3 0.9999740           0.14
#>   1287_94          12       3 0.9998988           0.13
#>  2945_100           7       3 0.9994826           0.18
#>  4038_100           7       5 0.9982678           0.04
#>  4327_100           7       3 0.9964061           0.33
#> ... 795 more nodes
```

codiv filters the symbiont tree’s internal nodes (by span, tip count,
and number of hosts), scans the survivors, and returns one row per
scanned node. Turn on the other methods
(`methods = c("hommola", "paco", "parafit")`) and
`subtree_features = TRUE` for a full analysis; those add runtime.

------------------------------------------------------------------------

### Running codiv()

#### Basic usage

``` r

codiv_results <- codiv(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  permutations = 99,
  Save_fp = "codiv_results.tsv"
)
```

#### Important parameters

| Parameter | Default | Meaning |
|----|----|----|
| `min_hosts` | 3 | Minimum number of different hosts in a symbiont subtree to include it in the scan (must be ≥ 3) |
| `min_symbiont_tips` | 7 | Minimum number of symbiont tips in a subtree (higher = more statistical power) |
| `span_fraction` | 0.1 | Only scan nodes with span ≤ 10% of max tree span (focuses on distal nodes) |
| `permutations` | 99 | Number of permutations for p-value calculation (higher = more precise) |
| `methods` | `c("hommola", "paco", "parafit")` | Codiversification methods (distance-based defaults), plus opt-in “topology” (branching-order test) |
| `cores` | all cores − 1 | Number of cores; nodes are scanned in parallel (single core on Windows) |
| `subtree_features` | TRUE | Calculate tree distance metrics? |
| `focus_hosts` | NULL | Optional host names; adds a `<host>_PRESENT` TRUE/FALSE column per host, for sorting the results by hosts of interest |
| `batch_size` / `n_batches` | auto | Control checkpoint frequency: nodes per batch, or number of batches (results are saved after each batch) |
| `verbose` | TRUE | Print progress messages? |

Symbiont tips with no host association are dropped automatically, and
every node is scanned independently, so one problematic node cannot
abort the run.

#### Distance methods vs. the topology method

The three default methods (Hommola, PACo, ParaFit) are all
**distance-based**: they correlate pairwise phylogenetic *distances*, so
they are sensitive to branch lengths. Two trees with the same branching
order but different evolutionary rates will score less than a perfect
match.

The optional `"topology"` method asks a different question: **do the two
trees have the same branching order?**, ignoring branch lengths
entirely. It can flag congruence that the distance methods understate
when rates differ between hosts and symbionts.

``` r

codiv_results <- codiv(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  methods = c("hommola", "topology"),   # topology is opt-in
  permutations = 99
)
```

Because a tree-shape comparison needs one tip per host, the topology
method reduces each host to its **largest single-host clade**
(tie-broken by the tightest branch lengths). A host whose symbionts are
completely scattered has no coherent position and is excluded. Two
diagnostic columns make this visible:

- `Topology_hosts_used` — how many hosts contributed to the node’s test
- `Topology_fraction_kept` — the fraction of symbiont tips retained
  after reducing each host to its dominant clade

A low `Topology_fraction_kept` means the congruence rests on heavy
pruning; treat those nodes with caution, and note that nodes with fewer
than three usable hosts return `NA`.

#### Resuming interrupted runs

If your analysis stops before completing, resume it:

``` r

# Continue from where it stopped
codiv_results <- codiv(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  Save_fp = "codiv_results.tsv",  # Same filename
  continue = TRUE,                 # Pick up where we left off
  permutations = 99
)
```

------------------------------------------------------------------------

### Interpreting Results

#### Output structure

``` r

# Compact overview and per-method summary
codiv_results
summary(codiv_results)

# The full table (a data frame) for filtering / export
head(as.data.frame(codiv_results))

# Column meanings:
# Node_ID              - Symbiont tree node identifier
# N_Symbionts          - Number of tips in this subtree
# N_Hosts              - Number of different hosts in this subtree
# Symbiont_Colless     - Tree imbalance metric (higher = more imbalanced)
# Symbiont_Sackin      - Tree shape metric (related to tip distances)
# Hommola_r            - Correlation coefficient (-1 to 1)
# Hommola_pvalue       - Permutation p-value
# PACo_ss              - PACo goodness-of-fit sum of squares (if method = "paco")
# PACo_pvalue          - PACo permutation p-value (if method = "paco")
# ParaFitGlobal        - ParaFit global test statistic (if method = "parafit")
# ParaFit_pvalue       - ParaFit permutation p-value (if method = "parafit")
```

The object is still a data frame, so `$`, `[`, and `dplyr` verbs all
work; use `as.data.frame(codiv_results)` when you want a plain data
frame.

#### Collapsing within-host variation

By default, every method runs on a **collapsed** subtree. Before the
tests run, codiv collapses each monophyletic clade of symbionts from the
**same host species** down to a single representative tip. This prevents
inflated correlations that would otherwise arise when near-identical
MAGs — for example genomes assembled from different individuals of the
same host species — are counted as independent observations. If MAGs
from one host species are **not** monophyletic (they fall in separate
parts of the symbiont tree, reflecting genuinely distinct lineages), all
of them are retained.

`Hommola_r` / `Hommola_pvalue`, `PACo_ss` / `PACo_pvalue`, and
`ParaFitGlobal` / `ParaFit_pvalue` are all computed on this collapsed
subtree. Setting `uncollapsed_hommola = TRUE` additionally runs
Hommola’s test on the *uncollapsed* subtree and adds
`Uncollapsed_Hommola_r` / `Uncollapsed_Hommola_pvalue`; this is opt-in
because it treats deeply sampled hosts as pseudoreplicated and is the
most expensive statistic on large clades.

#### Identifying codiversified nodes

Following [Sanders et
al. (2023)](https://www.nature.com/articles/s41564-023-01388-w) and
[Sprockett et
al. (2025)](https://www.nature.com/articles/s41467-025-57435-z), a clade
is called co-diversifying when it clears both a statistic and a p-value
threshold, e.g. `Hommola_r > 0.75` and `Hommola_pvalue < 0.01`:

``` r

library(dplyr)

significant_nodes <- as.data.frame(codiv_results) %>%
  filter(Hommola_r > 0.75,             # strong topological association
         Hommola_pvalue < 0.01) %>%
  arrange(desc(Hommola_r))

head(significant_nodes)
```

#### Scan-wide false-discovery control

To gauge false discoveries for a whole scan, codiv implements a
**second-order permutation test**: permute the host tree’s tip labels,
rerun the entire scan, and build a null distribution for the *number* of
significant clades. This is the appropriate approach because nested
clades in a tree are not independent, so a per-node multiple-testing
correction (e.g. Benjamini-Hochberg) would be invalid here.
[`codiv_null_scans()`](https://www.sprockettlab.com/codiv/reference/codiv_null_scans.md)
does this:

``` r

null <- codiv_null_scans(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  n_permutations = 100,          # host-label permutations (outer null)
  stat_threshold = 0.75,         # r > 0.75 ...
  p_threshold = 0.01,            # ... and p < 0.01 define "co-diversifying"
  permutations = 999, methods = "hommola", cores = 8
)

null$observed        # significant clades in the real data
mean(null$null_counts)  # expected under permuted host labels
null$global_pvalue   # is the scan more co-diversifying than chance?
null$empirical_fdr   # expected false-discovery proportion at these thresholds
```

Because each permutation reruns the full scan, this captures exactly the
dependence structure among nested clades. It is the recommended way to
report significance for a scan.

#### Tree shape metrics

The Colless and Sackin indices describe tree balance: - **Colless
index:** Higher values = more imbalanced tree (bifurcations uneven) -
**Sackin index:** Higher values = longer average path to tips

**Why it matters:** Different tree shapes can affect codiversification
patterns. A perfectly balanced host tree with a highly imbalanced
symbiont subtree might show weak correlation for reasons unrelated to
coevolution.

------------------------------------------------------------------------

### Secondary Analyses

The codiv package provides functions for deeper investigation of
codiversification patterns:

#### Leave-One-Host-Out Analysis

Assess which hosts drive the codiversification signal:

``` r

loo_results <- loo_host_analysis(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  permutations = 99
)

# View importance scores
loo_results$host_importance

# Interpret: Hosts with high impact_score contribute disproportionately
# to the codiversification signal
```

**Faster option:** Test specific hosts only:

``` r

# Just test 3 hosts (much faster)
loo_results <- loo_host_analysis(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  hosts_to_test = c("host_A", "host_B", "host_C"),
  permutations = 99
)
```

#### Parameter Sensitivity Analysis

How robust are your results to parameter choices?

``` r

sensitivity <- sensitivity_analysis(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  span_fraction_range = c(0.05, 0.1, 0.2),
  min_symbiont_tips_range = c(5, 7, 10),
  permutations = 99
)

# View sensitivity across parameters
sensitivity$sensitivity_data

# Strong codiversification should be robust across reasonable parameter ranges
```

#### Visualizing Host-Symbiont Trees

Draw a tanglegram: the two trees face to face with a line for each
association. This gives an intuitive sense of how congruent the trees
are. Requires the `ggplot2` package and returns a ggplot object you can
customize or save with `ggsave()`.

``` r

# Lines colored by host (default)
plot_codiv_trees(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont
)

# Highlight symbionts in significant nodes (needs a codiv result)
plot_codiv_trees(
  Host_tree = Host_Tree,
  Symbiont_tree = Symbiont_Tree,
  Host_to_Symbiont_df = Host_to_Symbiont,
  codiv_results = codiv_results,
  color_by = "significance",
  significance_threshold = 0.05,
  title = "Codiversification Patterns"
)
```

By default a cladogram layout is used (tips aligned, nodes spread by
topology) and the symbiont tree is rotated to reduce line crossings. Set
`use_branch_lengths = TRUE` to position nodes by branch length, or
`show_symbiont_labels = TRUE` to label symbiont tips.

------------------------------------------------------------------------

### Getting Help

#### Troubleshooting

**“Tree must be binary”**

``` r

# Solution: Collapse polytomies
library(ape)
Host_Tree <- multi2di(Host_Tree)
Symbiont_Tree <- multi2di(Symbiont_Tree)
```

**“Symbiont not found in tree”**

``` r

# Check spelling and exact names
head(Symbiont_Tree$tip.label)    # Actual names in tree
head(Host_to_Symbiont$Symbiont)  # Names in linking table
# They must match exactly!
```

**Analysis running very slowly**

Runtime scales with the number of nodes scanned and permutations. The
biggest levers, roughly in order of impact:

``` r

codiv(..., subtree_features = FALSE)   # skip the tree-distance metrics (often the
                                       #   most expensive step)
codiv(..., methods = "hommola")        # run one method instead of all three
codiv(..., permutations = 9)           # fewer permutations while exploring
codiv(..., min_symbiont_tips = 10,     # scan fewer, larger nodes ...
      span_fraction = 0.05)            #   ... focused on distal clades
codiv(..., cores = 8)                  # more cores (nodes scan in parallel;
                                       #   single core on Windows)
```

For very large trees (thousands of symbiont tips) a full run of all
three methods with `subtree_features = TRUE` can take hours and several
GB of RAM; consider an HPC node, or start with the hommola-only,
metrics-off configuration above.

**Repeated `Found more than one class "phylo" in cache` messages**

This harmless message comes from R’s `methods` package when two of
codiv’s dependencies both register an S4 `phylo` class. It does not
affect results. To silence it for good, add this to your `~/.Rprofile`:

``` r

.First <- function() {
  globalCallingHandlers(message = function(m) {
    txt <- conditionMessage(m)
    if (grepl("more than one class .+ in cache", txt) ||
        startsWith(txt, "Also defined by")) invokeRestart("muffleMessage")
  })
}
```

#### Documentation & References

- **Function help:**
  [`?codiv`](https://www.sprockettlab.com/codiv/reference/codiv.md),
  [`?loo_host_analysis`](https://www.sprockettlab.com/codiv/reference/loo_host_analysis.md),
  [`?plot_codiv_trees`](https://www.sprockettlab.com/codiv/reference/plot_codiv_trees.md)
- **Package documentation:** <https://sprockettlab.github.io/codiv/>
- **Methods papers:**
  - Hommola et al. 2009 (Biology Letters)
  - Balbuena et al. 2013 (Journal of Biogeography)

#### Reporting Issues

Found a bug? Have a feature request?

Visit: <https://github.com/SprockettLab/codiv/issues>

Please include: - Minimal reproducible example - Error message (full
traceback) - Session info:
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)

------------------------------------------------------------------------

### Session Information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ape_5.8-1   codiv_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] ade4_1.7-24        tidyselect_1.2.1   dplyr_1.2.1        castor_1.8.7      
#>  [5] farver_2.1.2       S7_0.2.2           fastmap_1.2.0      promises_1.5.0    
#>  [9] PlotTools_0.4.0    shinyjs_2.1.1      XML_3.99-0.24      digest_0.6.39     
#> [13] mime_0.13          lifecycle_1.0.5    cluster_2.1.8.2    magrittr_2.0.5    
#> [17] compiler_4.6.1     rlang_1.3.0        sass_0.4.10        progress_1.2.3    
#> [21] tools_4.6.1        TreeTools_2.4.0    igraph_2.3.3       yaml_2.3.12       
#> [25] knitr_1.52         prettyunits_1.2.0  bit_4.6.0          plyr_1.8.9        
#> [29] xml2_1.6.0         RColorBrewer_1.1-3 RNeXML_2.4.11      purrr_1.2.2       
#> [33] desc_1.4.3         grid_4.6.1         xtable_1.8-8       colorspace_2.1-3  
#> [37] ggplot2_4.0.3      scales_1.4.0       MASS_7.3-65        cli_3.6.6         
#> [41] vegan_2.7-6        rmarkdown_2.32     crayon_1.5.3       ragg_1.5.2        
#> [45] generics_0.1.4     otel_0.2.0         RSpectra_0.16-2    httr_1.4.9        
#> [49] reshape2_1.4.5     naturalsort_0.1.3  cachem_1.1.0       stringr_1.6.0     
#> [53] splines_4.6.1      parallel_4.6.1     vctrs_0.7.3        Matrix_1.7-5      
#> [57] TreeDist_2.14.1    jsonlite_2.0.0     hms_1.1.4          bit64_4.8.6       
#> [61] phylobase_0.8.12   pbmcapply_1.5.1    seqinr_4.2-44      systemfonts_1.3.2 
#> [65] jquerylib_0.1.4    tidyr_1.3.2        glue_1.8.1         pkgdown_2.2.1     
#> [69] stringi_1.8.9      gtable_0.3.6       later_1.4.8        paco_0.5.0        
#> [73] tibble_3.3.1       pillar_1.11.1      htmltools_0.5.9    R6_2.6.1          
#> [77] textshaping_1.0.5  Rdpack_2.6.6       evaluate_1.0.5     shiny_1.14.0      
#> [81] lattice_0.22-9     rbibutils_2.4.1    httpuv_1.6.17      bslib_0.12.0      
#> [85] adegenet_2.1.11    Rcpp_1.1.2         rncl_0.8.10        uuid_1.2-2        
#> [89] fastmatch_1.1-8    permute_0.9-10     nlme_3.1-169       mgcv_1.9-4        
#> [93] xfun_0.60          fs_2.1.0           adephylo_1.1-17    pkgconfig_2.0.3
```
