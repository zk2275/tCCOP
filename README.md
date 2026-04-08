# Conformal Inference of Individualized Treatment Rules under Distributional Shift

This repository contains code for the paper:

> **Conformal Inference of Individualized Treatment Rules under Distributional Shift**  
> Zhiyu Sui\*, Zhuodiao Kuang\*, Lu Tang, Ying Ding  
> Department of Biostatistics and Health Data Science, University of Pittsburgh

## Overview

We propose **tCOPP** (transferred Conformal Off-Policy Prediction), a conformal inference framework for individualized off-policy prediction under distributional shift. Unlike population-level policy evaluation, our methods provide **context-based outcome prediction intervals** under any given decision rule, with **finite-sample coverage guarantees**.

Key features:

- **Individualized prediction**: constructs prediction intervals tailored to patient characteristics, rather than a single population-level estimate.
- **Distributional shift correction**: handles covariate shift between source (e.g., RCT) and target (e.g., EHR) populations via importance weighting.
- **Multi-source integration**: borrows strength from multiple external data sources using federated transfer learning with data-adaptive aggregation.
- **Distribution-free**: compatible with any predictive algorithm; no distributional assumptions on outcomes.

## Methods

We develop three methods under the tCOPP framework:

| Method | Description | Data Used |
|---|---|---|
| **tCOPP-W** | Weighted conformal off-policy prediction | Internal source + Target |
| **tCOPP-EIF(I)** | Efficient influence function-based, internal source only | Internal source + Target |
| **tCOPP-EIF(M)** | Efficient influence function-based, multiple sources | Internal + External sources + Target |

### Baselines

| Method | Description |
|---|---|
| **SM** | Subsampling method without covariate shift adjustment |
| **COPP** | Conformal off-policy prediction (Zhang et al., 2023) without cross-population covariate shift correction |

### Key Assumptions

1. **Standard causal assumptions**: consistency, overlap, and unconfoundedness within each source dataset.
2. **Cross-population assumptions**: population overlap and common conditional outcome distribution (CCOD) between the internal source and the target.

## Repository Structure

```
├── README.md
├── R/                        # Core method implementations
│   ├── tCOPP_W.R             # tCOPP-W algorithm
│   ├── tCOPP_EIF_I.R         # tCOPP-EIF(I) algorithm
│   ├── tCOPP_EIF_M.R         # tCOPP-EIF(M) algorithm
│   ├── baselines.R           # SM and COPP baselines
│   ├── utils.R               # Shared utility functions
│   └── weight_estimation.R   # Density ratio and propensity score estimation
├── simulation/               # Simulation studies
│   ├── sim_internal_only.R   # Simulations using internal source only (Figure 1)
│   ├── sim_multi_source.R    # Simulations with multiple sources (Figures 2–3)
│   ├── dgp.R                 # Data generating processes
│   └── run_sim.sh            # SLURM submission script for HPC
├── application/              # Real data analysis
│   ├── sepsis_analysis.R     # ProCESS trial + UPMC EHR analysis
│   └── data_preprocessing.R  # Data cleaning and covariate construction
├── figures/                  # Scripts and outputs for reproducing figures
└── results/                  # Simulation and analysis outputs
```

> **Note**: The structure above reflects the intended organization. Some files may be updated as the manuscript undergoes revision.

## Simulation Studies

### Data Generating Process

Covariates $X_1, X_2, X_3, X_4 \sim U(0,1)$ are independently generated. Potential outcomes are:

$$Y(t) = 1 + X_1 - X_2 + [S + \mathbf{1}(S=0)]X_3^2 + \exp(X_4) + t(3 - 5X_1 + 2X_2 - 3X_3 + X_4) + (1+t)(1 + \textstyle\sum_{j=1}^4 X_j)\epsilon$$

where $\epsilon \sim N(0,1)$. Under this design, CCOD holds between the internal source $D_1$ and target $D_0$.

A more challenging scenario with complex CCOD violations for external sources is also evaluated (see Table 1 in the paper).

### Settings

- **Number of sources**: $K = 5$ (1 internal + 4 external)
- **Source sample size**: ~2,000 per source (75% training, 25% calibration)
- **Target sample size**: ~10,000
- **Replications**: 50
- **Nominal coverage**: $1 - \alpha = 0.9$

### Conformal Scores

- **CQR** (Conformalized Quantile Regression): used for SM, COPP, and tCOPP-W
- **Locally weighted ASR** (Absolute Standardized Residual): used for tCOPP-EIF(I), tCOPP-EIF(M), and tCOPP-EIF(P)

### Evaluation Metrics

- **Coverage Probability (CP)**: proportion of test samples whose true outcome falls within the prediction interval
- **Average Length (AL)**: mean length of prediction intervals across test samples

## Application: Sepsis Treatment Outcome Prediction

We apply the proposed methods to predict individualized 60-day in-hospital mortality under Early Goal-Directed Therapy (EGDT) for septic shock patients.

- **Source data**: ProCESS randomized controlled trial (EGDT arm: $n = 439$; usual care arm: $n = 456$), divided into 3 geographic regions
- **Target data**: UPMC electronic health records ($n = 6{,}546$ septic patients)
- **Covariates**: 21 baseline variables including demographics, vital signs, and laboratory results
- **Analysis configurations**:
  - All ProCESS data as internal source → tCOPP-W, tCOPP-EIF(I)
  - Region 1 as internal source; Regions 2 & 3 as external → tCOPP-EIF(I), tCOPP-EIF(M)

## Requirements

- R (≥ 4.0)
- Key packages: `quantregForest`, `SuperLearner`, `glmnet`, `ranger`, `parallel`

## References

- Zhang, Y., Shi, C., and Luo, S. (2023). Conformal off-policy prediction. *AISTATS*.
- Liu, Y., Levis, A.W., Normand, S.L., et al. (2024). Multi-source conformal inference under distribution shift. *arXiv:2405.09331*.
- Tibshirani, R.J., Barber, R.F., Candes, E.J., et al. (2019). Conformal prediction under covariate shift. *NeurIPS*.
- Romano, Y., Patterson, E., and Candes, E. (2019). Conformalized quantile regression. *NeurIPS*.
- ProCESS Investigators (2014). A randomized trial of protocol-based care for early septic shock. *NEJM*, 370(18), 1683–1693.

## Contact

- Zhiyu Sui — [ZHS66@pitt.edu](mailto:ZHS66@pitt.edu)
- Zhuodiao Kuang — [zhk22@pitt.edu](mailto:zhk22@pitt.edu)
- Lu Tang — [LUTANG@pitt.edu](mailto:LUTANG@pitt.edu)
- Ying Ding — [YINGDING@pitt.edu](mailto:YINGDING@pitt.edu)

## Citation

```bibtex
@article{sui2025tcopp,
  title={Conformal Inference of Individualized Treatment Rules under Distributional Shift},
  author={Sui, Zhiyu and Kuang, Zhuodiao and Tang, Lu and Ding, Ying},
  journal={},
  year={2025}
}
```

## License

This project is licensed under the MIT License.
