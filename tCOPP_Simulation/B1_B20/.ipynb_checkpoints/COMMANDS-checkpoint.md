# tCOPP simulation grid — terminal command reference

A quick reference for the commands used to run the 60-box simulation grid on the
Pitt CRC HTC cluster and to (re)generate the figures. Commands are grouped by
what you are trying to do; the two you run most often are marked ⭐.

Files this refers to (keep them in one directory, next to `tCOPP_functions.R`):
`submit_grid.sbatch`, `run_grid_array.R`, `tCOPP_grid_driver.R`,
`tCOPP_MuSCI_pooled_patch.R`, `tCOPP_functions.R`, `make_boxplots.R`.

---

## 0. Log in and go to the project folder

```bash
ssh yourusername@htc.crc.pitt.edu        # VPN first if off-campus; confirm host on CRC site
cd "tCOPP_Simulation_June/B=1&B=20"       # QUOTES are required: the '&' in the name
                                          # otherwise runs 'cd .../B=1' in the background
pwd                                       # confirm you are where you think you are
ls                                        # confirm the scripts are here
```

Optional, less error-prone: rename the folder so the `&` stops biting you.
```bash
cd /ihome/wchen/zhk22
mv "tCOPP_Simulation_June/B=1&B=20" tCOPP_Simulation_June/B1_B20
cd tCOPP_Simulation_June/B1_B20
```

---

## 1. Smoke test BEFORE the full run (recommended)

Runs one block (scenario b, B=20) with only 4 replicates, straight through
`Rscript` (no scheduler). Doubles as a syntax/parse check.
```bash
TCOPP_NREP=4 Rscript run_grid_array.R 4
```
Then verify the influence-function cells produced real numbers (not NA):
```bash
awk -F, 'NR>1{t[$4]++; if($6=="NA")n[$4]++} END{for(m in t)printf "%-14s %d/%d NA\n",m,(n[m]+0),t[m]}' results/grid_b_B20.csv
```
You want `0/…` NA for all five methods.

---

## 2. ⭐ Submit the full simulation (SLURM job array)

Six array tasks = 3 CCOD scenarios × {B=1, B=20}, each writing one CSV.
```bash
sbatch submit_grid.sbatch                 # prints: Submitted batch job <JOBID>
```

---

## 3. Monitor the jobs

```bash
squeue -u $USER                           # list your running/pending tasks
watch -n 30 squeue -u $USER               # auto-refresh every 30s (Ctrl-C to stop)
tail -f logs/tcopp_grid_<JOBID>_1.err     # live log of task 1 (Ctrl-C to stop watching)
```
The job is finished when `squeue -u $USER` shows no more `tcopp_gr` rows.
(`sbatch` returns immediately — it does NOT wait for the run to finish, so do not
chain the plotting step onto the same line.)

Cancel a job (e.g. a leftover buggy submission):
```bash
scancel <JOBID>                           # whole job
scancel <JOBID>_4                         # a single array task
```

---

## 4. Check the jobs succeeded (after the queue clears)

```bash
ls -l results/                            # expect 6 files: grid_{a,b,c}_B{01,20}.csv
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS,ExitCode   # want COMPLETED everywhere
cat logs/tcopp_grid_<JOBID>_1.err         # read errors for one task if something failed
```
Look for `FAILED`, `TIMEOUT`, or `OUT_OF_MEMORY`; each task also prints `N/N reps ok`.

---

## 5. ⭐ Make the figures

Reads the six CSVs and writes `figures/fig_scenario_{a,b,c}.(png|pdf)`, `fig_MASTER`,
and the summary tables. NOTE: this reads existing results, so if you only changed
the plotting code (e.g. the y-axis ranges), run ONLY this — no need to resubmit.
```bash
Rscript make_boxplots.R
```

---

## Common gotchas we hit

- **`&` in the folder name** splits `cd` into a background job → "No such file or
  directory". Fix: quote the path, or rename the folder (Section 0).
- **`sbatch: error: Unable to open file`** → you are not in the folder with the
  `.sbatch` file (the `cd` silently failed). Re-check with `pwd` / `ls`.
- **`GLIBCXX_3.4.32 not found` (osqp fails, tasks die in seconds)** → the batch
  script must NOT `module load gcc/...` (an older libstdc++ shadows conda's). It
  should only activate the conda env that passed the smoke test. Already fixed in
  `submit_grid.sbatch`.
- **"100/100 reps ok" but EIF cells are NA** → the rep count is not the success
  test; always run the `awk` NA check in Section 1. (Fixed: `alpha` is now defined
  for the influence-function estimators.)
- **`No grid_*.csv found`** when plotting → the jobs have not finished yet; wait
  for `squeue` to clear and `ls results/` to show all six CSVs.

## Environment knobs

```bash
export TCOPP_NREP=100      # Monte Carlo replicates per cell (set inside submit_grid.sbatch)
export OMP_NUM_THREADS=1   # keep BLAS single-threaded; parallelism is over replicates
```
`sbatch` resources live in the `#SBATCH` header of `submit_grid.sbatch`
(`--array=1-6`, `--cpus-per-task`, `--mem`, `--time`, `--partition=htc`).
Figure y-axis ranges live at the top of `make_boxplots.R` (`CP_YLIM`, `AL_YLIM`).
