#!/usr/bin/env Rscript
## =====================================================================
## run_grid_array.R  --  entry point for the 60-cell grid on Pitt HPC
## ---------------------------------------------------------------------
## One SLURM array task computes one (scenario, B) block = 10 cells x nrep reps
## and writes a single tidy CSV. Six tasks cover the whole grid:
##
##   task 1 -> (a, B=1)     task 2 -> (a, B=20)
##   task 3 -> (b, B=1)     task 4 -> (b, B=20)
##   task 5 -> (c, B=1)     task 6 -> (c, B=20)
##
## Replicates within a task run in parallel with mclapply across the cores
## SLURM allocated to the task.
##
## Usage on the cluster (see submit_grid.sbatch):
##     Rscript run_grid_array.R          # reads $SLURM_ARRAY_TASK_ID
## Local / interactive:
##     Rscript run_grid_array.R 4        # run task 4 only
##     Rscript run_grid_array.R all 20   # run ALL tasks, nrep = 20 (quick test)
## =====================================================================

suppressMessages({
  library(MASS); library(numDeriv); library(rootSolve); library(osqp)
  library(densratio); library(dplyr); library(glmnet); library(quantreg)
  library(randomForest); library(grf); library(SuperLearner); library(parallel)
})

## ---- source the library, the pooled patch, then the driver (order matters)
source("tCOPP_functions.R")
source("tCOPP_MuSCI_pooled_patch.R")
source("tCOPP_grid_driver.R")

## ---- how many Monte-Carlo replicates (paper uses 50; grid files used 100) ----
NREP <- as.integer(Sys.getenv("TCOPP_NREP", unset = "100"))

## ---- the 6-task map -------------------------------------------------
TASKS <- data.frame(
  scenario = c("a","a","b","b","c","c"),
  B        = c( 1, 20,  1, 20,  1, 20),
  stringsAsFactors = FALSE
)

## ---- resolve which task(s) to run -----------------------------------
args <- commandArgs(trailingOnly = TRUE)
slurm_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))

if (length(args) >= 1 && tolower(args[1]) == "all") {
  task_ids <- seq_len(nrow(TASKS))
  if (length(args) >= 2) NREP <- as.integer(args[2])
} else if (length(args) >= 1 && !is.na(suppressWarnings(as.integer(args[1])))) {
  task_ids <- as.integer(args[1])
} else if (!is.na(slurm_id)) {
  task_ids <- slurm_id
} else {
  stop("No task id: set SLURM_ARRAY_TASK_ID, or pass a task number / 'all' on the command line.")
}

## ---- cores ----------------------------------------------------------
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                unset = Sys.getenv("SLURM_CPUS_ON_NODE", unset = NA)))
if (is.na(ncores)) ncores <- max(1, detectCores() - 1)

out_dir <- "results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ---- run each requested task ----------------------------------------
for (tid in task_ids) {
  sc <- TASKS$scenario[tid]; B <- TASKS$B[tid]
  message(sprintf("\n==== task %d : scenario (%s), B=%d, nrep=%d, cores=%d ====",
                  tid, sc, B, NREP, ncores))
  t0 <- Sys.time()

  reps <- mclapply(seq_len(NREP),
                   function(r) tryCatch(run_one_rep_grid(r, scenario = sc, B = B),
                                        error = function(e) {
                                          message(sprintf("  rep %d FAILED: %s", r,
                                                          conditionMessage(e))); NULL }),
                   mc.cores = ncores, mc.preschedule = FALSE)

  ok  <- !vapply(reps, is.null, logical(1))
  message(sprintf("  %d/%d reps ok", sum(ok), NREP))
  df  <- do.call(rbind, reps[ok])

  fout <- file.path(out_dir, sprintf("grid_%s_B%02d.csv", sc, B))
  write.csv(df, fout, row.names = FALSE)
  message(sprintf("  wrote %s  (%d rows) in %.1f min",
                  fout, nrow(df), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

message("\nAll requested tasks done.")
