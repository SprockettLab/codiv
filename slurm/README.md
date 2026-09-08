# HPC Support for codiv

Optional SLURM scripts for running codiv on high-performance computing clusters.

**Important:** HPC support is **optional** — codiv works perfectly well on a laptop for small to medium trees. Only use these scripts if you have very large trees (1000+ tips) or want to run multiple analyses in parallel.

---

## Quick Start

1. **Edit the configuration section** in `submit_codiv_single.sh`
   - Set paths to your trees and data
   - Adjust parameters (permutations, methods, etc.)
   - Set resource requirements based on tree size (see table below)

2. **Submit the job:**
   ```bash
   sbatch submit_codiv_single.sh
   ```

3. **Monitor progress:**
   ```bash
   squeue -u $USER
   tail -f logs/codiv_*.log
   ```

4. **Retrieve results:**
   ```bash
   cat results/codiv_results.tsv
   ```

---

## Resource Requirements

Scale resources based on your tree size:

| Tree Size | Cores | RAM | Time | Example Data |
|-----------|-------|-----|------|--------------|
| Small (50-200 tips) | 4 | 8-16 GB | 1-4 hours | Pilot studies |
| Medium (200-1000 tips) | 8 | 16-32 GB | 4-12 hours | Typical gut microbiomes |
| Large (1000-5000 tips) | 16 | 32-64 GB | 12-48 hours | MAG collections |
| Very Large (5000+ tips) | 32 | 64-128 GB | 48+ hours | All MAGs from database |

**Notes:**
- Time scales non-linearly with tree size
- More permutations (999 vs 99) ~10x runtime increase
- Multiple methods (all 3 vs 1) ~2-3x runtime increase
- Parallelization gains plateau around 16 cores

---

## Available Scripts

### `submit_codiv_single.sh`

Run a single codiv analysis on a cluster.

**Features:**
- Single job submission
- Automatic logging
- Configurable parameters
- Works with standard module systems

**Usage:**
```bash
# Edit configuration
nano submit_codiv_single.sh

# Submit
sbatch submit_codiv_single.sh

# Check status
squeue -u $USER

# View output
tail -f logs/codiv_*.log
```

**Output:**
- Results: `results/codiv_results.tsv`
- Log: `logs/codiv_*.log`
- Error: `logs/codiv_*.err`

---

## Workflow: Running Multiple Analyses

### Scenario: Testing different parameters

**Option 1: Submit multiple jobs**
```bash
# Test different span fractions
for SPAN in 0.05 0.1 0.2; do
  sed "s/SPAN_FRACTION=.*/SPAN_FRACTION=$SPAN/" submit_codiv_single.sh > \
    submit_codiv_span_$SPAN.sh
  sbatch submit_codiv_span_$SPAN.sh
done

# Monitor all jobs
watch squeue -u $USER
```

**Option 2: Use job arrays** (coming in future release)
- Run the same analysis on different data subsets
- Manage many jobs through a single `sbatch` call

### Scenario: Long-running analysis with checkpointing

Use the `continue` parameter in codiv() to resume interrupted runs:

```bash
# First run (interrupted)
sbatch submit_codiv_single.sh

# Later, resubmit with continue=TRUE
# The function will pick up where it left off
sbatch submit_codiv_single.sh
```

---

## Cluster-Specific Configuration

### Module Loading

If your cluster uses modules, add to the script before the R call:

```bash
# Common examples:
module load R/4.5
module load intel/2021
module load gcc/11.2
```

### Queue Selection

Adjust SLURM partition/queue for your cluster:

```bash
#SBATCH --partition=gpu       # For GPU cluster
#SBATCH --partition=long      # For long-running jobs
#SBATCH --partition=priority  # High-priority queue
```

### Email Notifications

Get notified when job completes:

```bash
#SBATCH --mail-user=your.email@institution.edu
#SBATCH --mail-type=END,FAIL
```

### Account/Project Billing

If your cluster requires account specification:

```bash
#SBATCH --account=my_project
```

---

## Troubleshooting

### Job won't start / waiting in queue

**Check available resources:**
```bash
sinfo                  # View partitions and available nodes
squeue -p partition    # See queue for specific partition
```

**Solution:** Reduce requested resources or try different partition:
```bash
#SBATCH --cpus-per-task=4    # Try 4 instead of 8
#SBATCH --partition=gpu      # Try different queue
```

### Job times out

**Problem:** Allocated time wasn't enough

**Solution:** Increase time limit:
```bash
#SBATCH --time=24:00:00      # Increase to 24 hours
```

**Or optimize code:**
```bash
PERMUTATIONS=9               # Reduce for quick test
METHODS="hommola"            # Use only fastest method
```

### Memory errors

**Problem:** Job killed due to insufficient memory

**Solution:** Increase memory request:
```bash
#SBATCH --mem=64G            # Increase to 64 GB
```

### "Cannot find module R"

**Problem:** R module not available

**Solution:** Load R manually or check available modules:
```bash
module avail R               # See available R versions
module load R/4.5            # Load specific version
```

---

## Best Practices

1. **Test locally first:**
   ```bash
   # Run on small subset locally before submitting to cluster
   results_local <- codiv(..., span_fraction=0.2, permutations=9)
   ```

2. **Check logs before trusting results:**
   ```bash
   tail logs/codiv_*.log
   ```

3. **Use checksums for data integrity:**
   ```bash
   md5sum $HOST_TREE $SYMBIONT_TREE
   # Save checksums; verify them in output
   ```

4. **Document your run:**
   ```bash
   # Create a record of parameters used
   cat > results/run_parameters.txt << EOF
   Date: $(date)
   Job ID: $SLURM_JOB_ID
   Host tree: $HOST_TREE
   Symbiont tree: $SYMBIONT_TREE
   Parameters: permutations=$PERMUTATIONS, span=$SPAN_FRACTION
   EOF
   ```

---

## For HPC Administrators

If deploying codiv for a research group:

1. **Install as module:**
   - Install codiv system-wide or in a shared environment
   - Create a module file pointing to the installation
   - Users load with `module load codiv`

2. **Set resource defaults:**
   - Create a template script with your cluster's standard settings
   - Document typical resource requests per tree size

3. **Enable job monitoring:**
   - Set up Slurm accounting to track codiv resource usage
   - Helps users optimize requests for future runs

---

## When to Use HPC vs. Local

### Use **Local Machine:**
- Tree size < 1000 tips
- Quick exploratory analysis
- Parameter testing
- Learning codiv

### Use **HPC Cluster:**
- Tree size > 1000 tips  
- High-resolution permutation testing (999+)
- Running multiple analyses in parallel
- Production-level results for publication

---

## Additional Resources

- **SLURM Documentation:** https://slurm.schedmd.com/sbatch.html
- **R on HPC:** https://researchcomputing.org.uk/software/r/
- **Cluster-specific help:** Check your institution's HPC documentation

### `submit_congruence_sweep.sh`

Runs the methods-paper congruence sweep (`simulations/run_congruence_sweep.R`):
node-level localization AUC vs. realized congruence, for the uncollapsed and
collapsed statistics. Results append to a CSV after every dataset, so the job is
resumable — resubmit with the same `OUT_CSV` to continue.

Configure with environment variables (`REPS`, `NHOSTS`, `CORES`, `PERMS`,
`NPERHOST`). Runtime scales with `REPS x congruence-values x PERMS`; memory
scales with `NHOSTS` and `n_per_host`.
