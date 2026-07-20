# tCOPP simulation grid — 60 boxplots (COPP / tCOPP-W / EIF(I) / EIF(M) / EIF(P) × {B=1,B=20} × {internal, multi} × {a,b,c})

This bundle answers the three questions and builds the full 60-cell simulation.
It is designed to drop into your existing repo alongside `tCOPP_functions.R`
and run on Pitt CRC HTC. **Nothing here has been executed** (no R in the
authoring environment) — the code was written against your existing files and
bracket-checked, but please smoke-test with a tiny `TCOPP_NREP` first.

---

## 1. Which methods in the code are "ours"? (method-name map)

Your paper's methods are **SM, COPP, tCOPP-W, tCOPP-EIF(I), tCOPP-EIF(M), tCOPP-EIF(P)**.
The three we *propose* are **tCOPP-W, tCOPP-EIF(I), tCOPP-EIF(M)** (Sec. "Estimation and
Algorithms"); tCOPP-EIF(P) is the **pooled oracle baseline**, and SM/COPP are prior-work
baselines. The `.Rmd` files use legacy labels, mapped below.

| Paper name        | Role                 | Code function (`tCOPP_functions.R`)         | Key args                                          | Legacy plot label / result var                          |
|-------------------|----------------------|---------------------------------------------|---------------------------------------------------|---------------------------------------------------------|
| SM                | baseline             | `Naive()`                                   | —                                                 | `"SM"` / `CP_low_Naive`                                 |
| COPP              | baseline (Zhang '23) | `ConformalCI()`                             | `method="BootSweight"`, `weighted=TRUE`, wttl via `w.method=1` (unit weights) | `"COPP"` / `CP_low_BootSweight`     |
| **tCOPP-W**       | **proposed**         | `ConformalCI()`                             | same, but wttl via `w.method=2` (real covariate-shift weights, with intercept) | `"COPP-TL"` / `CP_low_BootSweight_tl` |
| **tCOPP-EIF(I)**  | **proposed**         | `MuSCI_ConformalCI()` (single source)       | `conf.score="localASR"`, uses `rhat_tgSt` → `PI.tgSt` | `"COPP-EIF-(local)ASR"` / `CP_low_EIF_localASR_BootSweight` |
| **tCOPP-EIF(M)**  | **proposed (flagship)** | `MuSCI_ConformalCI_MuS()` (multi source)  | `conf.score="localASR"`, uses `rhat_fed2` → `PI.fed2` | `"COPP-EIF-MSource-(local)ASR"` / `CP_low_EIF_MSource_localASR_BootSweight` |
| tCOPP-EIF(P)      | oracle (pooled)      | `MuSCI_ConformalCI_MuS()` **`rhat_pool`→`PI.pool`** | same fit, pooled interval instead of federated | *(not surfaced by the original code — see §2)*          |
| COPP-IS           | extra baseline       | `ConformalCI()`                             | `method="BootDweight"` (double weighting)         | `"COPP-IS"` / `CP_low_BootDweight`                      |
| COPP-IS-TL        | extra baseline       | `ConformalCI()`                             | `method="BootDweight"` + `w.method=2`             | `"COPP-IS-TL"` / `CP_low_BootDweight_tl`                |

So **your reading is correct**: `COPP-TL` = **tCOPP-W**, `COPP-EIF-…-ASR` (target only) =
**tCOPP-EIF(I)**, and `COPP-EIF-MSource-localASR` = **tCOPP-EIF(M)**. Two clarifications:

- COPP vs tCOPP-W differ **only** in the covariate-shift weight: COPP calls
  `learn_weights(X, Xtest, w.method=1, …)` (returns all-ones → subsampling correction only),
  while tCOPP-W calls `learn_weights(cbind(1,X), cbind(1,Xtest), w.method=2, …)` (real
  transfer weights). Both use `weighted=TRUE` in `ConformalCI`.
- **`COPP-IS` / `COPP-IS-TL` (BootDweight / "double weighting") are baselines only** and are
  *not* among the four reported methods — don't conflate Dweight with Sweight.

## 2. Why is there no tCOPP-EIF(P) in the "Pooled" figure files?

Because **nothing in those files ever extracts the pooled interval**:

1. `MuSCI_ConformalCI_MuS()` **computes** `PI.pool` (from `rhat_pool`) internally but
   `return`s only the federated `PI.fed2` — see the two lines
   `YBleft <- rbind(YBleft, PI.fed2$lower)` / `…PI.fed2$upper`. The pooled result is
   discarded before the function returns.
2. `Fig_2a___Pooled_pc_B1.Rmd` has **no pooled experiment at all**. It runs the same seven
   experiments as the single-source/multi-source files and its final `PlotData` lists exactly
   `SM, COPP, COPP-IS, COPP-TL, COPP-IS-TL, COPP-EIF-localASR, COPP-EIF-MSource-localASR`.
   "Pooled" in the title is a leftover.
3. `tCOPP_Fig2b_pooled_result_pc_B1.Rmd` *does* contain an `Experiment 8: EIF-Pooled-localASR`
   block (the `CP_low_EIF_localASR_BootSweight_pooled` runner), but **every line of it is
   commented out**, so it never executes and never enters the plot.

The manuscript's Fig. 2(b)/(c) *do* report tCOPP-EIF(P) (filenames `…_pooled.png`), so those
panels were made by an older/manual path, not by the code as uploaded. **The patch in this
bundle fixes the root cause** by exposing the pooled interval so EIF(P) can be produced
consistently.

## 3. The 60-box grid (Figure-2 replacement)

Layout (your spec): each method appears **once, in its natural source-mode**, on a
single left-to-right axis -- the manuscript's mixed-axis "internal source only and
multiple sources" figure -- with a **B=1 row and a B=20 row**, and CP + AL side by side:

```
        |          Coverage prob.           |            Avg. length           |
  B=1   | COPP tCOPP-W EIF(I) EIF(M) EIF(P) | COPP tCOPP-W EIF(I) EIF(M) EIF(P)
  B=20  | COPP tCOPP-W EIF(I) EIF(M) EIF(P) | COPP tCOPP-W EIF(I) EIF(M) EIF(P)
```

= 5 CP + 5 AL per B row x 2 B rows = 20 boxes/scenario x 3 scenarios = **60 boxes**.
Natural modes: COPP / tCOPP-W / tCOPP-EIF(I) -> internal source D1 only;
tCOPP-EIF(M) / tCOPP-EIF(P) -> multiple sources (D1 + external). Boxes are coloured
by mode so the caption is visible. `make_boxplots.R` writes
`fig_scenario_a/b/c.(png|pdf)` for the three `\includegraphics` slots plus a stacked
`fig_MASTER`.

**Scenarios** (only the external-site outcome model changes; D1 <-> D0 always satisfy
CCOD): **(a)** external sites share the internal model (CCOD holds for all); **(b)**
external site k scales `X3^3` by k (simple violation, your Fig 2b model); **(c)** the
complex `f_k,g_k,h_k` model from the manuscript table.

**Two expected, correct artefacts to be aware of:**
- COPP / tCOPP-W / tCOPP-EIF(I) use only D1, so their boxes are (up to Monte-Carlo
  noise) **the same across (a), (b), (c)** -- only EIF(M)/EIF(P) respond to the
  external-source scenario. That invariance is itself informative.
- **tCOPP-W at B=20** is run with the native BootSweight majority-vote aggregation --
  the "configuration parity" option your multi-sampling section already sanctions.
  Because tCOPP-W's threshold is a *stable weighted quantile*, its B=1 and B=20 boxes
  should look nearly identical, which is exactly the contrast that proves aggregation
  is **necessary for the EIF thresholds but not for the weighted method**. So the box
  is not "meaningless" -- it's the control that carries the B=20 argument. (If you'd
  rather flag it, `make_boxplots.R` can grey or annotate it -- ask.)

Efficiency & integrity: each replicate generates ONE dataset and fits each model ONCE;
EIF(M) and EIF(P) are read off the SAME multi-source fit (federated vs pooled
interval), so those two cells are paired and consistent. Cost is ~4 model fits/rep.

### Files
- `tCOPP_MuSCI_pooled_patch.R` -- overrides the two MuSCI functions to also return the
  pooled interval (`out$*_pool`), unlocking tCOPP-EIF(P); federated/target stay the
  default so existing `.Rmd`s keep working. Auto-generated from your `tCOPP_functions.R`
  by `_gen_patch.py` (rerun that if you edit the library).
- `tCOPP_grid_driver.R` -- scenario DGPs + the 5 natural-mode cells + faithful CP/AL.
- `run_grid_array.R` / `submit_grid.sbatch` -- one SLURM task per (scenario, B), parallel over reps.
- `make_boxplots.R` -- the three (a)/(b)/(c) figures + master + summary CSVs.

### Run
```bash
# 0. put these next to your existing tCOPP_functions.R (uses localASR for the EIF methods)
# 1. smoke test first (this doubles as the parse-check I could not do here):
TCOPP_NREP=4 Rscript run_grid_array.R 4        # scenario (b), B=20
# 2. full grid on CRC:
sbatch submit_grid.sbatch                       # 6 array tasks -> results/grid_*.csv
# 3. figures once all six CSVs exist:
Rscript make_boxplots.R                         # -> figures/fig_scenario_{a,b,c}.png
```

### Settled per your reply
- Grid **replaces** the old (a)/(b)/(c) (Figure 2); narrative updated once results land.
- Both **B=1 and B=20** shown for all methods (to demonstrate B=20 is necessary).
- **localASR** for the EIF methods (consistent with the paper); CQR for COPP/tCOPP-W.
- `ex_site_size = 2000`, `nrep = 100`.
- tCOPP-EIF(P) recomputed from scratch via the patch -- report whatever it gives.

### Not covered here
- **Figure 1** (internal-source-only, which also includes **SM**) is a *separate* figure
  and is not regenerated by this grid (the grid's 5 methods exclude SM). Say if you want
  Figure 1 rebuilt too.
