# Input/Output Guide for codiv

This guide details the format requirements for codiv inputs and explains how to interpret the outputs.

---

## Part 1: Input Data Requirements

### Overview

codiv requires three inputs:

1. **Host phylogenetic tree** (phylo object)
2. **Symbiont phylogenetic tree** (phylo object)
3. **Host-to-Symbiont linking data** (data frame with Host and Symbiont columns)

### 1. Phylogenetic Trees

#### Format Requirements

Trees must be in R's `phylo` class (from the `ape` package). Standard formats:

```r
library(ape)

# From Newick format (.nwk, .tre)
tree <- read.tree("tree.nwk")

# From Nexus format (.nex, .nexus)
tree <- read.nexus("tree.nex")

# From other formats (phangorn, treeio)
library(phangorn)
tree <- read.phyDat("alignment.fasta", format = "fasta", type = "DNA")
tree <- nj(dist.ml(tree))  # Construct tree from alignment
```

#### Requirements

- **Binary/fully resolved:** Every internal node must have exactly 2 descendants
  - Non-binary trees will cause errors
  - Solution: `ape::multi2di(tree)` to resolve polytomies
  
- **Unique tip labels:** All tip names must be unique within each tree
  - Used to match host/symbiont records
  - Must match exactly (case-sensitive) in host-symbiont linking data
  
- **No missing branch lengths:** Recommended but not required
  - Used for patristic distance calculations
  - If missing, all edges treated as length 1
  
- **Node labels (optional):** codiv will auto-generate if missing
  - If provided, must be unique
  - Used as identifiers in output

#### Checking Your Trees

```r
# Check if binary
ape::is.binary(host_tree)  # Should be TRUE
ape::is.binary(symbiont_tree)  # Should be TRUE

# Check for unique tip labels
length(unique(host_tree$tip.label)) == length(host_tree$tip.label)
length(unique(symbiont_tree$tip.label)) == length(symbiont_tree$tip.label)

# Repair if needed
host_tree <- ape::multi2di(host_tree)
symbiont_tree <- ape::multi2di(symbiont_tree)
```

### 2. Host-to-Symbiont Linking Data

#### Format Specification

A data frame with exactly two columns: **"Host"** and **"Symbiont"**

```r
# Correct format
host_symbiont_links <- data.frame(
  Host = c("human_001", "human_002", "human_003", "human_001"),
  Symbiont = c("Bacteroides_sp1", "Bacteroides_sp2", "Faecalibacterium_sp1", "Bacteroides_sp3")
)

head(host_symbiont_links)
#        Host         Symbiont
# 1 human_001  Bacteroides_sp1
# 2 human_002  Bacteroides_sp2
# 3 human_003 Faecalibacterium_sp1
# 4 human_001  Bacteroides_sp3
```

#### Column Requirements

| Column | Requirements | Notes |
|--------|--------------|-------|
| **Host** | One per row; must match tree tip label | Each host can appear multiple times (multiple symbionts) |
| **Symbiont** | One per row; must match tree tip label | For MAG data: each symbiont appears exactly once (one-to-one to host) |

#### Data Types

- Column names: character (exact match: "Host", "Symbiont")
- Values: character (tip labels from trees)
- No missing values allowed

#### Common Formats to Convert

**From CSV file:**
```r
host_symbiont_links <- read.csv("associations.csv")
colnames(host_symbiont_links)  # Check names match "Host" and "Symbiont"
```

**From Excel:**
```r
library(readxl)
host_symbiont_links <- read_excel("associations.xlsx", sheet = 1)
```

**From tab-separated file:**
```r
host_symbiont_links <- read.delim("associations.txt")
```

#### Validation

```r
# Check structure
str(host_symbiont_links)

# Check column names
all(c("Host", "Symbiont") %in% colnames(host_symbiont_links))

# Check for missing values
!any(is.na(host_symbiont_links))

# Check all hosts are in tree
all(host_symbiont_links$Host %in% host_tree$tip.label)

# Check all symbionts are in tree
all(host_symbiont_links$Symbiont %in% symbiont_tree$tip.label)

# For MAG data, check each symbiont appears once
sum(duplicated(host_symbiont_links$Symbiont)) == 0
```

---

## Part 2: Understanding Outputs

### Output Structure

The main output is a data frame with one row per symbiont tree node and numerous columns:

```r
results <- codiv(Host_tree, Symbiont_tree, Host_to_Symbiont_df, ...)
head(results)
ncol(results)  # Many columns; see below for explanation
```

### Column Descriptions

#### Node Identification

| Column | Type | Meaning |
|--------|------|---------|
| `Node_ID` | character | Unique identifier for this node in symbiont tree |
| `N_Symbionts` | integer | Number of tip descendants in this subtree |
| `N_Hosts` | integer | Number of unique hosts represented in this subtree |

#### Correlation Coefficients (if method = "hommola")

| Column | Type | Range | Interpretation |
|--------|------|-------|-----------------|
| `Hommola_r` | numeric | -1 to 1 | Correlation between symbiont and host distances |
| `Hommola_pvalue` | numeric | 0 to 1 | Statistical significance from permutation test |
| `Collapsed_Hommola_r` | numeric | -1 to 1 | Correlation after collapsing monophyletic symbionts |
| `Collapsed_Hommola_pvalue` | numeric | 0 to 1 | P-value for collapsed analysis |

**Interpretation:**
- **r > 0.5:** Strong codiversification signal
- **0.3 < r < 0.5:** Moderate codiversification
- **0.1 < r < 0.3:** Weak codiversification
- **r < 0.1:** No clear signal

#### PACo Results (if method = "paco")

| Column | Type | Meaning |
|--------|------|---------|
| `Collapsed_PACo_ss` | numeric | Sum of squared residuals (goodness-of-fit) |
| `Collapsed_PACo_pvalue` | numeric | Significance from permutation test |

**Interpretation:** Lower values = better fit (more codiversification)

#### ParaFit Results (if method = "parafit")

| Column | Type | Meaning |
|--------|------|---------|
| `Collapsed_ParaFitGlobal` | numeric | Global ParaFit statistic |
| `Collapsed_ParaFit_pvalue` | numeric | Significance from permutation test |

#### Tree Shape Metrics

| Column | Type | Meaning |
|--------|------|---------|
| `Symbiont_Colless` | numeric | Colless index (tree imbalance) |
| `Symbiont_Sackin` | numeric | Sackin index (tip distance variance) |
| `Host_Colless` | numeric | Colless index for host subtree |
| `Host_Sackin` | numeric | Sackin index for host subtree |

**Why it matters:**
- Unbalanced trees can show apparent correlations even without coevolution
- Use these to assess whether correlations are biological or structural artifacts

#### Tree Distance Metrics (if subtree_features = TRUE)

Multiple tree comparison metrics comparing host and symbiont topologies:

- `Collapsed_TreeDistance`: Standard Robinson-Foulds distance
- `Collapsed_SharedPhylogeneticInfo`: Information shared between trees
- `Collapsed_DifferentPhylogeneticInfo`: Information unique to each tree
- `Collapsed_NyeSimilarity`: Similarity measure (0-1)
- `Collapsed_JaccardRobinsonFoulds`: Jaccard index of splits
- `Collapsed_MatchingSplitDistance`: Distance based on matching splits
- `Collapsed_MatchingSplitInfoDistance`: Information-theoretic split distance
- `Collapsed_MutualClusteringInfo`: Mutual information of clusters

#### Focus Hosts (if focus_hosts specified)

| Column | Type | Values | Meaning |
|--------|------|--------|---------|
| `[HostName]_PRESENT` | logical | TRUE/FALSE | Whether specified host is in this subtree |

Example: If `focus_hosts = c("human_001", "human_002")`, you get columns:
- `human_001_PRESENT`
- `human_002_PRESENT`

#### Tree Newick Strings (if subtree_features = TRUE)

| Column | Type | Contents |
|--------|------|----------|
| `Symbiont_Tree` | character | Newick format of symbiont subtree |
| `Host_Tree` | character | Newick format of host subtree (pruned to relevant hosts) |

### Interpreting Results

#### Identifying Codiversified Nodes

```r
# Find statistically significant codiversification
significant <- results %>%
  filter(Hommola_pvalue < 0.05,        # Significant
         Hommola_r > 0.3,              # Biologically meaningful
         N_Symbionts >= 7,             # Adequate sample size
         N_Hosts >= 3)                 # Multiple hosts

# Sort by strength of signal
significant <- significant %>%
  arrange(Hommola_pvalue, desc(Hommola_r))

head(significant)
```

#### Filtering by Node Size

```r
# Only large subtrees (more statistical power)
large_nodes <- results %>%
  filter(N_Symbionts >= 15,
         N_Hosts >= 5)

# Only small subtrees (recent divergence)
small_nodes <- results %>%
  filter(N_Symbionts >= 5,
         N_Symbionts <= 10)
```

#### Consensus Across Methods

If running multiple methods (hommola + paco + parafit), look for agreement:

```r
# Extract p-values from all methods
results_filtered <- results %>%
  filter(Hommola_pvalue < 0.05,
         Collapsed_PACo_pvalue < 0.05,
         Collapsed_ParaFit_pvalue < 0.05)  # All methods agree

# These are your most robust codiversification signals
```

#### Comparing Correlation Strength

```r
# Visualize distribution
hist(results$Hommola_r, breaks = 50, main = "Distribution of r values")

# Summary statistics
summary(results$Hommola_r[results$Hommola_pvalue < 0.05])
```

---

## Part 3: Common Issues & Solutions

### Input Issues

#### Error: "Symbiont_tree contains tips that are not in Host_to_Symbiont_df"

**Problem:** Tree tip labels don't match linking data

**Solution:**
```r
# Check what's in each
setdiff(symbiont_tree$tip.label, host_symbiont_links$Symbiont)  # In tree but not data
setdiff(host_symbiont_links$Symbiont, symbiont_tree$tip.label)  # In data but not tree

# Common causes:
# - Typos in names
# - Different naming conventions (underscores vs dots)
# - Extra characters (whitespace, quotes)

# Fix: Clean up names
host_symbiont_links$Symbiont <- trimws(host_symbiont_links$Symbiont)  # Remove whitespace
```

#### Error: "Tree must be binary"

**Problem:** Tree has polytomies (nodes with >2 descendants)

**Solution:**
```r
# Check polytomies
any(ape::branching.times(symbiont_tree) == 0)  # TRUE if polytomies present

# Resolve polytomies (randomly)
symbiont_tree <- ape::multi2di(symbiont_tree)
```

#### Error: "min_hosts must be >= 3" or similar parameter error

**Problem:** Parameter values outside valid range

**Solution:**
```r
# Valid ranges:
# min_hosts: >= 3 (default 3)
# min_symbiont_tips: >= 1 (default 7, higher recommended)
# span_fraction: > 0 and <= 1 (default 0.1)
# permutations: >= 1 (default 99)
```

### Output Issues

#### Too Few Nodes Scanned

**Problem:** Only a handful of nodes meet filtering criteria

**Causes:**
- `min_hosts` is too high (many nodes don't have enough host diversity)
- `min_symbiont_tips` is too high (few large clades)
- `span_fraction` is too small (only scanning very distal nodes)

**Solution:**
```r
# Relax criteria
results <- codiv(...,
  min_hosts = 2,              # Lower threshold
  min_symbiont_tips = 5,      # Lower threshold
  span_fraction = 0.2)        # Scan more nodes
```

#### All P-values = 1

**Problem:** No significant codiversification detected (could be real or parameter issue)

**Causes:**
- Symbionts truly don't codiversify with hosts
- `permutations` is too low (use 99 minimum)
- Sample size too small

**Solution:**
```r
# Try more permutations for better precision
results <- codiv(..., permutations = 999)
```

#### Negative Correlations

**Problem:** Some nodes show r < 0

**Meaning:** Symbiont topology is *opposite* to host topology (polyphyletic groups in one tree scattered in the other)

**Interpretation:** This is biologically valid; indicates:
- No coevolution at that node
- Or possible host-switching events (for negative r)

---

## Part 4: Data Preparation Examples

### From Alignment to Analysis

```r
library(ape)
library(phangorn)

# 1. Read alignment
alignment <- read.phyDat("host_sequences.fasta", format = "fasta", type = "DNA")

# 2. Construct tree
dist <- dist.ml(alignment)
tree <- nj(dist)

# 3. Verify it's binary
tree <- ape::multi2di(tree)

# 4. Root tree (if needed)
tree <- ape::root(tree, outgroup = "outgroup_species")

# 5. Use in codiv analysis
results <- codiv(host_tree, symbiont_tree, host_symbiont_links, ...)
```

### From QIIME2 Output

```r
library(qiime2R)

# Load QIIME2 artifacts
physeq <- qza_to_phyloseq("table.qza", "rooted_tree.qza", "taxonomy.qza")
host_symbiont_links <- read.csv("sample_metadata.csv")

# Extract trees
host_tree <- phy_tree(physeq)
symbiont_tree <- refseq(physeq)  # Or separate ASV tree

# Prepare linking data
colnames(host_symbiont_links) <- c("Host", "Symbiont")

# Analyze
results <- codiv(host_tree, symbiont_tree, host_symbiont_links, ...)
```

### From RDP or SILVA Classification

```r
# After creating OTU table and taxonomic assignment
library(phyloseq)

# Create phyloseq object
physeq <- phyloseq(otu_table, tax_table, sample_data, phy_tree)

# Extract components
host_tree <- sample_data(physeq)  # Or metadata tree
symbiont_tree <- phy_tree(physeq)

# Format associations
host_symbiont_links <- data.frame(
  Host = colnames(otu_table),
  Symbiont = rownames(otu_table)
)

results <- codiv(...)
```

---

## Part 5: Saving and Exporting Results

### Save Results as CSV

```r
# Save full results
write.csv(results, "codiv_results.csv", row.names = FALSE)

# Save significant results only
significant <- results %>%
  filter(Hommola_pvalue < 0.05, Hommola_r > 0.3)
write.csv(significant, "codiv_significant_nodes.csv", row.names = FALSE)
```

### Extract for Visualization

```r
# Prepare for plotting
plot_data <- results %>%
  select(Node_ID, N_Symbionts, N_Hosts, 
         Hommola_r, Hommola_pvalue,
         Symbiont_Colless, Host_Colless) %>%
  mutate(Significant = Hommola_pvalue < 0.05)

# Correlation plot
plot(results$Symbiont_Colless, results$Hommola_r,
     col = ifelse(results$Hommola_pvalue < 0.05, "red", "gray"),
     xlab = "Symbiont Tree Imbalance (Colless)",
     ylab = "Correlation Coefficient")
```

### Export Newick Trees

```r
# Extract and save individual trees
for (i in seq_len(nrow(results))) {
  tree_str <- results$Symbiont_Tree[i]
  node_id <- results$Node_ID[i]
  cat(tree_str, "\n", file = paste0("trees/", node_id, ".nwk"), append = TRUE)
}
```

---

## Additional Resources

- **R phylogenetics:** ape, phangorn, phyloseq packages
- **File formats:** http://ape-package.ird.fr/
- **Tree visualization:** ggtree, phytools
- **Advanced:** TreeDist, castor packages used by codiv internally
