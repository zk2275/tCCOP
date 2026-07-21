#!/bin/bash
#SBATCH --job-name=tcopp_grid
#SBATCH --array=1-6                 # 6 tasks = 3 scenarios x {B=1, B=20}
#SBATCH --cpus-per-task=32          # reps run in parallel via mclapply
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/tcopp_grid_%A_%a.out
#SBATCH --error=logs/tcopp_grid_%A_%a.err
# ---------------------------------------------------------------------
# Pitt CRC HTC submission for the 60-cell tCOPP simulation grid.
#
# Each array task computes one (scenario, B) block and writes
# results/grid_<scenario>_B<BB>.csv. After all six finish, run:
#     Rscript make_boxplots.R
#
# Adjust module/conda lines to your environment. TCOPP_NREP sets the
# number of Monte-Carlo replicates (default 100).
# ---------------------------------------------------------------------

mkdir -p logs results

# --- load R (edit to match your CRC setup) ---
module purge
module load gcc/12.2.0 r/4.4.0        # example; use `module spider r` to find yours
# If you use conda instead:
#   source activate r-tcopp

export TCOPP_NREP=100                 # replicates per cell (paper used 50)
export OMP_NUM_THREADS=1              # keep BLAS single-threaded; parallelism is over reps

echo "host=$(hostname)  array_task=${SLURM_ARRAY_TASK_ID}  cpus=${SLURM_CPUS_PER_TASK}"
Rscript run_grid_array.R              # reads SLURM_ARRAY_TASK_ID
