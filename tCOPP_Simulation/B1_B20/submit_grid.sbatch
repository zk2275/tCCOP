#!/bin/bash
#SBATCH --job-name=tcopp_grid
#SBATCH --partition=htc              # HTC partition (edit if your allocation differs)
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
# IMPORTANT (do NOT `module load gcc/...` here): the conda R packages
# (osqp, grf, ...) are compiled against conda's newer libstdc++. Loading
# a gcc module puts an OLDER libstdc++ on the path and osqp.so then fails
# with "GLIBCXX_3.4.32 not found". We therefore use ONLY the conda env
# that passed the interactive smoke test -- no `module` commands.
# TCOPP_NREP sets the number of Monte-Carlo replicates (default 100).
# ---------------------------------------------------------------------

mkdir -p logs results

# --- activate the SAME conda env used in the smoke test (R + all packages) ---
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate base
# If 'conda' is somehow not on PATH in the batch job, replace the two lines
# above with the literal path, e.g.:
#   source /ihome/wchen/zhk22/miniconda3/etc/profile.d/conda.sh
#   conda activate base

export TCOPP_NREP=100                 # replicates per cell (paper used 50)
export OMP_NUM_THREADS=1              # keep BLAS single-threaded; parallelism is over reps

echo "host=$(hostname)  array_task=${SLURM_ARRAY_TASK_ID}  cpus=${SLURM_CPUS_PER_TASK}"
Rscript run_grid_array.R              # reads SLURM_ARRAY_TASK_ID
