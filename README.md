# Conformal Inference of Individualized Treatment Rules under Distributional Shift

Code accompanying the paper:

> **Conformal Inference of Individualized Treatment Rules under Distributional Shift**
> Zhiyu Sui\*, Zhuodiao Kuang\*, Lu Tang, Ying Ding
> Department of Biostatistics and Health Data Science, University of Pittsburgh

We propose **tCOPP** (transferred Conformal Off-Policy Prediction), a conformal
inference framework for individualized off-policy prediction under distributional
shift. It produces **context-based outcome prediction intervals** under a given
evaluation policy, correcting for covariate shift between source and target and
(optionally) integrating multiple external sources via federated aggregation.

## Repository structure

This tree reflects the files actually present in the repository.

```
.
├── README.md
├── LICENSE
├── tCOPP_Simulation/                  # Simulation study (paper Section 3)
│   ├── tCOPP_functions.R              # Core method implementations (all estimators)
│   ├── tCOPP_MuSCI_pooled_patch.R     # Surfaces the pooled-oracle EIF(P) interval
│   ├── tCOPP_grid_driver.R            # Single-cell driver of the simulation grid
│   ├── run_grid_array.R               # Array runner over scenarios × B × sources
│   ├── submit_grid.sbatch             # SLURM submission (Pitt CRC HTC)
│   ├── make_boxplots.R                # Builds scenario figures from results/
│   ├── _gen_patch.py                  # Helper that generates the pooled patch
│   ├── COMMANDS.md                    # Exact run commands
│   ├── README.md                      # Method-name map (paper <-> code) + notes
│   ├── results/                       # Raw grid output
│   │   └── grid_{a,b,c}_B{01,20}.csv  # scenario a/b/c, B = 1 and B = 20
│   └── figures/                       # Figures + summaries used in the paper
│       ├── fig_scenario_{a,b,c}.{png,pdf}
│       ├── grid_summary_medians.csv   # Median CP / AL per method x scenario x B
│       └── results_combined.csv
└── tCOPP_RWD/                         # Real-data analysis (paper Section 4)
    ├── COPP_realdata_NoStrip_new.Rmd      # Single-source: tCOPP-W, tCOPP-EIF(I)
    └── multisource_fixed_patch8_new.Rmd   # Multi-source: tCOPP-EIF(M) + pooled EIF(P)
```

## Methods: paper name vs. code

All estimators live in `tCOPP_Simulation/tCOPP_functions.R`. The `.Rmd`/driver
files use legacy labels; the mapping (also in `tCOPP_Simulation/README.md`) is:

| Paper name       | Role                         | Code entry point                                | Legacy label                  |
|------------------|------------------------------|-------------------------------------------------|-------------------------------|
| **tCOPP-W**      | proposed                     | `ConformalCI(..., w.method=2)`                  | `COPP-TL`                     |
| **tCOPP-EIF(I)** | proposed                     | `MuSCI_ConformalCI()` (`conf.score="localASR"`) | `COPP-EIF-(local)ASR`         |
| **tCOPP-EIF(M)** | proposed (flagship)          | `MuSCI_ConformalCI_MuS()` -> `rhat_fed2`        | `COPP-EIF-MSource-(local)ASR` |
| COPP             | baseline (Zhang et al. 2023) | `ConformalCI(..., w.method=1)`                  | `COPP`                        |
| tCOPP-EIF(P)     | pooled-oracle baseline       | `MuSCI_ConformalCI_MuS()` -> `rhat_pool`        | via `tCOPP_MuSCI_pooled_patch.R` |

Conformal scores (see paper Appendix B): the simulation uses **CQR** for COPP /
tCOPP-W and **locally weighted ASR** for the EIF estimators; the real-data
analysis uses the **absolute-residual (AR)** score. In the code, the keyword
`conf.score = "ASR"` is the *unweighted* AR score and `"localASR"` is the
locally weighted version.

## Data availability

The real-data analysis uses the **ProCESS randomized trial** and **UPMC EHR**
data, which are **restricted** (proprietary / patient privacy) and are **not**
distributed in this repository. The `tCOPP_RWD/*.Rmd` files document the analysis
pipeline but require access to those datasets to run. See the manuscript's Data
Availability statement for details. Do **not** commit patient-level data; the
included `.gitignore` blocks the relevant file patterns.

## Reproducing the simulation study

```bash
cd tCOPP_Simulation
# see COMMANDS.md for the exact grid; on an HPC cluster:
sbatch submit_grid.sbatch          # runs run_grid_array.R across the grid
Rscript make_boxplots.R            # rebuilds figures/fig_scenario_{a,b,c}
```

Outputs land in `results/` (raw CP/AL per replicate) and `figures/`
(`grid_summary_medians.csv` gives the median CP/AL per method x scenario x B
reported in the paper). Data-generating processes and settings (K = 5 sources,
~2,000 per source, ~10,000 target, nominal 1 - alpha = 0.9) are described in
paper Section 3.1.

## Requirements

- R (>= 4.0)
- R packages: `quantregForest`, `SuperLearner`, `glmnet`, `ranger`, `parallel`

## Citation

```bibtex
@article{sui_tcopp,
  title   = {Conformal Inference of Individualized Treatment Rules under Distributional Shift},
  author  = {Sui, Zhiyu and Kuang, Zhuodiao and Tang, Lu and Ding, Ying},
  journal = {Biostatistics (under review)},
  year    = {2026}
}
```

## Contact

- Zhuodiao Kuang — zhk22@pitt.edu
- Zhiyu Sui — ZHS66@pitt.edu
- Lu Tang — LUTANG@pitt.edu
- Ying Ding — YINGDING@pitt.edu

## License

MIT License (see `LICENSE`).
