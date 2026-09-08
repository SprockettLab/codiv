#!/bin/bash
#SBATCH --job-name=codiv_sweep
#SBATCH --output=logs/codiv_sweep_%j.log
#SBATCH --error=logs/codiv_sweep_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=default

# Congruence sweep for the codiv methods paper (localization AUC vs. realized
# congruence). Results append to a CSV after every dataset, so this job can be
# preempted and resubmitted to resume where it left off.
#
# Resource guidance: memory scales with tree size (NHOSTS) and n_per_host;
# runtime scales with REPS x number of congruence values x PERMS. For the
# defaults below (50 reps, 20 hosts, 999 permutations) budget several hours.
#
# Usage:
#   sbatch slurm/submit_congruence_sweep.sh
#
# To resume an interrupted run, just resubmit with the same OUT_CSV.

set -e
mkdir -p logs simulations

# --- configuration (edit for your run) ------------------------------------
export REPS=50
export NHOSTS=20
export CORES=8            # match --cpus-per-task above
export PERMS=999
export NPERHOST=1,4
export OUT_CSV=simulations/congruence_sweep.csv
export OUT_PNG=simulations/congruence_sweep.png

# --- environment (customize for your cluster) -----------------------------
# Do not use `conda activate` in SLURM scripts; prepend the env to PATH:
#   export PATH=/path/to/r-env/bin:$PATH
# or load a module:
#   module load R/4.5

echo "Starting congruence sweep at $(date) (job $SLURM_JOB_ID)"
echo "REPS=$REPS NHOSTS=$NHOSTS CORES=$CORES PERMS=$PERMS"

# codiv must be installed in the R library used here, e.g.:
#   R -e 'remotes::install_github("SprockettLab/codiv")'
Rscript simulations/run_congruence_sweep.R

echo "Finished at $(date)"
