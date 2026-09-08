# Beginner's Guide to codiv

Welcome! This guide will walk you through installing and using the codiv package step-by-step, even if you're new to R packages.

---

## What is codiv?

**codiv** analyzes how well the evolutionary history of microbes (symbionts) matches the evolutionary history of their hosts. It answers the question: *"Did these microbes and hosts evolve together?"*

For example:
- Do your patient's gut bacteria have a family tree that mirrors the evolutionary relationships between patients?
- Do your plant's fungal endophytes diversify similarly to how the plants diversify?

---

## Part 1: Installation

### Step 1: Install R (if you haven't already)

If you don't have R installed, download it from: https://cran.r-project.org/

Choose your operating system (Windows, macOS, or Linux) and follow the installer.

**Verify installation:** Open R and type:
```r
R.version
```

You should see version information. If not, R isn't installed correctly.

### Step 2: Install RStudio (Recommended)

RStudio is a user-friendly interface for R. Download from: https://www.rstudio.com/products/rstudio/download/

Choose "RStudio Desktop (Free)" and select your operating system.

### Step 3: Install codiv

Open RStudio (or R) and paste this command:

```r
# Install dependencies first
install.packages(c("ape", "castor", "tidyverse", "data.table"))

# Then install codiv from GitHub
remotes::install_github("SprockettLab/codiv")
```

This downloads and installs codiv and everything it needs.

**Note:** The first time you run this, R may ask: *"Do you want to install from sources the packages which need compilation?"* 
- Answer **yes** — this is normal and necessary.

### Step 4: Verify Installation

Paste this into R to confirm codiv is installed:

```r
library(codiv)
?codiv
```

If a help page appears describing the `codiv()` function, you're ready to go!

---

## Part 2: Finding Documentation

### Online Website

The full documentation lives at: **https://sprockettlab.github.io/codiv/**

On that site you'll find:
- **Getting Started** — quick introduction
- **Reference** — detailed function descriptions
- **Articles** — full workflows and examples

### In R

You can access documentation directly in R:

```r
# View main function documentation
?codiv

# View specific function help
?hommola
?check_inputs

# List all available functions
help(package = "codiv")
```

---

## Part 3: Prepare Your Data

Before running codiv, you need three things:

### 1. Host Tree (Newick format)

A file describing evolutionary relationships between hosts. Looks like:
```
(Host_A:0.5,(Host_B:0.3,Host_C:0.3):0.2):0;
```

Common sources:
- **From sequence alignment:** Use BEAST, RAxML, or FastTree
- **From database:** Download from TreeBASE, PhyloTree.org, or SILVA
- **For testing:** Generate random trees in R

**File formats accepted:** `.nwk`, `.newick`, `.tre`, `.txt`

### 2. Symbiont Tree (Newick format)

A file describing evolutionary relationships between symbionts (microbes). Same format as host tree.

### 3. Host-to-Symbiont Data File

A table linking which hosts have which symbionts. Create a CSV or TSV file like this:

```
Host_ID    Symbiont_ID
Host_A     Microbe_1
Host_A     Microbe_3
Host_B     Microbe_2
Host_B     Microbe_3
Host_C     Microbe_1
```

**Column names matter:** Use exactly `Host_ID` and `Symbiont_ID` (you can change these in the function call, but these are defaults).

### Test Data (Easy Option)

codiv includes test data. Load it with:

```r
library(codiv)

# Load example data
data(host_tree_example)
data(symbiont_tree_example)
data(host_symbiont_links_example)

# This gives you three objects ready to use:
# - host_tree_example (a phylogenetic tree)
# - symbiont_tree_example (a phylogenetic tree)
# - host_symbiont_links_example (a data frame)
```

---

## Part 4: Run Your First Analysis

### Easiest Way (30 seconds)

```r
library(codiv)

# Load example data
data(host_tree_example)
data(symbiont_tree_example)
data(host_symbiont_links_example)

# Run analysis
results <- codiv(
  Host_tree = host_tree_example,
  Symbiont_tree = symbiont_tree_example,
  Host_to_Symbiont_df = host_symbiont_links_example
)

# Look at results
head(results)
```

### With Your Own Data

```r
library(codiv)
library(ape)

# Read your trees from files
host_tree <- read.tree("path/to/your/host_tree.nwk")
symbiont_tree <- read.tree("path/to/your/symbiont_tree.nwk")

# Read your linking data
links_df <- read.csv("path/to/your/links.csv")

# Run analysis
results <- codiv(
  Host_tree = host_tree,
  Symbiont_tree = symbiont_tree,
  Host_to_Symbiont_df = links_df
)
```

---

## Part 5: Understanding Your Results

### The Results Table

The `results` object is a data frame (spreadsheet-like table). View it with:

```r
# See first few rows
head(results)

# See full summary
summary(results)

# See all columns
names(results)
```

### Key Columns Explained

| Column | Meaning |
|--------|---------|
| `Node_ID` | Which symbiont node (internal branch point) we measured |
| `hommola_correlation` | How well this symbiont clade matches host diversification (0 to 1). **Higher = better match** |
| `hommola_p_value` | Statistical significance (0 to 1). **Lower = more significant** (p < 0.05 is "significant") |
| `N_hosts` | How many hosts have symbionts in this clade |
| `N_symbiont_tips` | How many symbiont species in this clade |

### Interpreting Correlation Coefficients

**Hommola correlation ranges from 0 to 1:**
- **0.8-1.0**: Strong codiversification (symbionts evolved tightly with hosts)
- **0.5-0.8**: Moderate codiversification
- **0.2-0.5**: Weak codiversification
- **0.0-0.2**: Little to no codiversification

**P-values:**
- **p < 0.05**: Statistically significant (likely real pattern, not random chance)
- **p ≥ 0.05**: Not statistically significant (could be random)

### Example Interpretation

```r
# Filter for significant results
significant_results <- results[results$hommola_p_value < 0.05, ]

# See strongest codiversification signals
best_matches <- results[order(-results$hommola_correlation), ][1:5, ]

# These 5 symbiont clades showed the best host-symbiont coevolution
```

---

## Part 6: Save and Export Results

### Save to a File

```r
# Save as CSV (spreadsheet format)
write.csv(results, "my_codiv_results.csv", row.names = FALSE)

# Save as TSV (tab-separated)
write.table(results, "my_codiv_results.tsv", 
            sep = "\t", row.names = FALSE)
```

### Save as R Object (for later use)

```r
# Save everything
saveRDS(results, "my_codiv_results.RDS")

# Load it later
results <- readRDS("my_codiv_results.RDS")
```

---

## Part 7: Run Tests (Verify Everything Works)

### Quick Test

Run this to make sure codiv is working correctly:

```r
# Load the package
library(codiv)

# Run quick check
data(host_tree_example)
data(symbiont_tree_example)
data(host_symbiont_links_example)

test_results <- codiv(
  Host_tree = host_tree_example,
  Symbiont_tree = symbiont_tree_example,
  Host_to_Symbiont_df = host_symbiont_links_example,
  permutations = 9  # Quick test: only 9 permutations instead of 99
)

# If you see results with no errors, codiv is working!
print(test_results)
```

### Full Test Suite

If you want to run comprehensive tests:

```r
# In R console
devtools::test()

# Or in RStudio:
# Keyboard shortcut: Ctrl+Shift+T (Windows/Linux) or Cmd+Shift+T (macOS)
```

You should see output like:
```
✓ test-check_inputs.R
✓ test-codiv.R
✓ test-hommola.R
...
```

If all tests pass (✓), the package is working correctly!

---

## Part 8: Common Questions

### "How long will analysis take?"

**Typical times** (on a laptop):
- Small trees (100 tips): 1-5 minutes
- Medium trees (500 tips): 10-60 minutes
- Large trees (1000+ tips): Hours to days

**Speed up your analysis:**
```r
# Option 1: Fewer permutations (less statistical power, but faster)
results <- codiv(..., permutations = 9)

# Option 2: Fewer methods (faster, but fewer correlation types)
results <- codiv(..., methods = "hommola")

# Option 3: Larger span_fraction (analyzes fewer nodes)
results <- codiv(..., span_fraction = 0.2)
```

### "I got an error about 'phylo object'"

This means your tree isn't in the right format.

**Fix:**
```r
# Make sure you used read.tree()
library(ape)
tree <- read.tree("your_tree.nwk")

# Check it's correct
class(tree)  # Should print "phylo"
```

### "My data file isn't being read"

**Common mistakes:**
1. Wrong file path — use full path or check working directory with `getwd()`
2. Wrong separator — CSV uses commas, TSV uses tabs
3. Column names don't match — must be exactly `Host_ID` and `Symbiont_ID`

**Fix:**
```r
# Check your working directory
getwd()

# Read with explicit settings
links <- read.csv("links.csv", stringsAsFactors = FALSE)

# Verify column names
colnames(links)  # Should show "Host_ID" and "Symbiont_ID"
```

### "Should I use hommola, paco, or parafit?"

**Hommola** (default): Standard codiversification metric, good all-purpose choice

**PACo**: Focuses on tree topology matching, good for host-specific patterns

**ParaFit**: Combines branch lengths and topology, most stringent test

**Best practice:** Try all three and compare:
```r
results <- codiv(..., methods = c("hommola", "paco", "parafit"))
```

---

## Part 9: Troubleshooting

### Installation Problems

**Problem:** `Error: package 'X' could not be found`

**Solution:** Install missing package:
```r
install.packages("X")  # Replace X with package name
library(codiv)
```

**Problem:** `remotes::install_github()` doesn't work

**Solution:** Install remotes first:
```r
install.packages("remotes")
remotes::install_github("SprockettLab/codiv")
```

### Runtime Problems

**Problem:** `Error: Host_tree is not a phylo object`

**Solution:** Your tree isn't loaded correctly:
```r
library(ape)
tree <- read.tree("your_file.nwk")
class(tree)  # Must print "phylo"
```

**Problem:** Analysis is running very slowly

**Solution:** Try these:
```r
# Reduce computational load
results <- codiv(
  ...,
  permutations = 9,      # Fewer permutations
  span_fraction = 0.15,  # Fewer nodes to scan
  methods = "hommola"    # Single method
)
```

---

## Part 10: Next Steps

### Learn More

1. **Read the full vignette:** In R, type `vignette("codiv")`
2. **Visit the website:** https://sprockettlab.github.io/codiv/
3. **Check input requirements:** Look for the Input/Output Guide on the website

### Advanced Features

Once comfortable with basics, try:

```r
# Leave-one-host-out analysis (test robustness)
loo_results <- loo_host_analysis(
  results,
  Host_tree = host_tree,
  Symbiont_tree = symbiont_tree,
  Host_to_Symbiont_df = links_df
)

# Sensitivity analysis (test parameter robustness)
sensitivity_results <- sensitivity_analysis(
  Host_tree = host_tree,
  Symbiont_tree = symbiont_tree,
  Host_to_Symbiont_df = links_df,
  span_range = c(0.05, 0.1, 0.2)
)
```

### Report Issues

Found a bug or problem? Create an issue at:
https://github.com/SprockettLab/codiv/issues

---

## Quick Reference Card

### Installation
```r
install.packages(c("ape", "castor", "tidyverse", "data.table"))
remotes::install_github("SprockettLab/codiv")
library(codiv)
```

### Quick Test
```r
data(host_tree_example)
data(symbiont_tree_example)
data(host_symbiont_links_example)
results <- codiv(host_tree_example, symbiont_tree_example, 
                 host_symbiont_links_example)
head(results)
```

### Your Own Data
```r
library(ape)
host <- read.tree("host.nwk")
symbiont <- read.tree("symbiont.nwk")
links <- read.csv("links.csv")
results <- codiv(host, symbiont, links)
```

### View Results
```r
head(results)              # First few rows
summary(results)           # Statistical summary
results[results$hommola_p_value < 0.05, ]  # Significant results only
```

### Save Results
```r
write.csv(results, "results.csv")
saveRDS(results, "results.RDS")
```

---

## Getting Help

1. **In R:** `?codiv`, `?hommola`, etc.
2. **Website:** https://sprockettlab.github.io/codiv/
3. **Vignette:** `vignette("codiv")`
4. **Issues:** https://github.com/SprockettLab/codiv/issues
5. **Email:** Daniel.Sprockett@wfusm.edu

---

Happy analyzing! 🧬
