#!/bin/bash
#SBATCH --job-name=codiv_analysis
#SBATCH --output=logs/codiv_%j.log
#SBATCH --error=logs/codiv_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=default

# Optional HPC support for codiv analysis
#
# This script runs a single codiv analysis on an HPC cluster.
# Resource requirements scale with tree size:
#   - Small trees (50-200 tips): 4 cores, 8-16 GB RAM, 1-4 hours
#   - Medium trees (200-1000 tips): 8 cores, 16-32 GB RAM, 4-12 hours
#   - Large trees (1000-5000 tips): 16 cores, 64 GB RAM, 12-48 hours
#
# Usage:
#   sbatch submit_codiv_single.sh
#
# Modify parameters below for your specific analysis

# ============================================================================
# USER CONFIGURATION
# ============================================================================

# Paths (update these for your data)
HOST_TREE="data/host_tree.nwk"
SYMBIONT_TREE="data/symbiont_tree.treefile"
HOST_SYMBIONT_DATA="data/host_symbiont_links.txt"
OUTPUT_FILE="results/codiv_results.tsv"

# Parameters
PERMUTATIONS=99          # Number of permutations (99 for quick, 999+ for publication)
MIN_HOSTS=3              # Minimum hosts per node
MIN_SYMBIONT_TIPS=7      # Minimum symbiont tips per node
SPAN_FRACTION=0.1        # Span threshold (0.1 = distal nodes only)
METHODS="hommola"        # Methods: hommola, paco, parafit (can combine)
CORES=8                  # Use all allocated cores
SEED=8675309             # Random seed for reproducibility

# ============================================================================
# SCRIPT (do not modify below this line)
# ============================================================================

set -e  # Exit on error

# Create log directory if it doesn't exist
mkdir -p logs results

echo "Starting codiv analysis at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Allocated resources: $SLURM_CPUS_PER_TASK CPUs, $SLURM_MEM_PER_NODE MB"
echo ""

# Load modules (customize for your cluster)
# module load R/4.5
# module load gcc/11.2.0

# Create R script
cat > /tmp/codiv_analysis_$SLURM_JOB_ID.R << 'RSCRIPT'
#!/usr/bin/env Rscript
library(codiv)

# Read parameters from environment
host_tree <- ape::read.tree(Sys.getenv("HOST_TREE"))
symbiont_tree <- ape::read.tree(Sys.getenv("SYMBIONT_TREE"))
host_symbiont_data <- read.csv(Sys.getenv("HOST_SYMBIONT_DATA"), sep="\t")

# Run codiv analysis
message("Running codiv analysis...")
results <- codiv(
  Host_tree = host_tree,
  Symbiont_tree = symbiont_tree,
  Host_to_Symbiont_df = host_symbiont_data,
  min_hosts = as.integer(Sys.getenv("MIN_HOSTS")),
  min_symbiont_tips = as.integer(Sys.getenv("MIN_SYMBIONT_TIPS")),
  span_fraction = as.numeric(Sys.getenv("SPAN_FRACTION")),
  permutations = as.integer(Sys.getenv("PERMUTATIONS")),
  seed = as.integer(Sys.getenv("SEED")),
  verbose = TRUE,
  Save_fp = Sys.getenv("OUTPUT_FILE"),
  methods = strsplit(Sys.getenv("METHODS"), ",")[[1]],
  cores = as.integer(Sys.getenv("CORES")),
  subtree_features = TRUE
)

message("Analysis complete!")
message(paste("Results saved to:", Sys.getenv("OUTPUT_FILE")))
RSCRIPT

# Export variables for R script
export HOST_TREE
export SYMBIONT_TREE
export HOST_SYMBIONT_DATA
export OUTPUT_FILE
export PERMUTATIONS
export MIN_HOSTS
export MIN_SYMBIONT_TIPS
export SPAN_FRACTION
export METHODS
export CORES
export SEED

# Run analysis
Rscript /tmp/codiv_analysis_$SLURM_JOB_ID.R

# Cleanup
rm /tmp/codiv_analysis_$SLURM_JOB_ID.R

echo ""
echo "Job completed at $(date)"
