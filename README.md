# Time Point Specific Fit Indices for DSEM: Simulation Study

**Master Thesis** | University of Tübingen | Prof. Holger Brandt

---

## Overview

This repository contains a comprehensive simulation study evaluating the performance of **time-specific fit indices** for **Dynamic Structural Equation Models (DSEM)** applied to Intensive Longitudinal Data (ILD).

### Research Objectives

1. **Evaluate timepoint-specific fit indices** (BRMSEA, BMc, B̂, $B̂_{adj}$, BCFI, BTLI) in detecting misspecifications in DSEM
2. **Assess detection performance** across varying conditions:
   - Sample sizes (N = 30, 60, 90, 120, 150, 300)
   - Number of time points (Nt = 5, 10, 15, 25, 50, 75)
   - Misspecification levels

---

## Key Innovation

Unlike traditional fit indices that aggregate over all time points, this study examines **time-point-specific fit indices** that can:
- Identify **when** (at which time point) model misspecification occurs
- Detect **localized misfit** in longitudinal structures
- Handle **time-varying dynamics** inherent in DSEM

---

## Repository Structure

```
Time-Point-Specific-Bayesian-Fit-Indices-for-DSEM/
├── README.md                          # This file
├── LICENSE                            # License
├── .gitignore                         # Git exclusions
├── main.R
├── run_cluster.sh                     # SLURM job submission
│
├── src/
│   ├── fitindices_utils.R         # Fit indices computation
│   └── dgp_sem.R                  # Data generation functions
├── models/
│   ├── DSEM_VAR1_RI.stan          # AR(1) latent growth model
│   └── DSEM_NULL.stan             # Null baseline model
└── results/
    └── 

```

---

## Installation & Setup

### Prerequisites

- **R** (≥ 4.0)
- **Stan** (RStan package)
- **Packages**: `rstan`, `tidyverse`, `loo`

### Quick Start

```r
# Install dependencies
install.packages()

# Clone repository
git clone https://github.com/[your-username]/dsem-fit-indices.git
cd dsem-fit-indices

# Run a simple example

```

---

## Simulation Design

### Study Conditions (Factorial Design)

| Parameter | Values |
|-----------|--------|
| **Sample size (N)** | 30, 60, 90, 120, 150, 300 |
| **Time points (Nt)** | 5, 10, 15, 25, 50, 75 |
| **Cross Loading Misspecification (M)** | 0, 0.3, 0.6 |
| **Ocurrence for event T (T_event)** | 0.25, 0.5, 0.75 |
| **Replications** | 250 per condition |


---

## Key Files Description

### Data Generation

**`scripts/dgp_sem.R`**
- `gendata01()`: Generates AR(1) latent growth model with within-time cross-loadings
- `gendata02()`: Generates model with cross-time cross-loadings
- Handles missing data (MCAR mechanism)
- Returns list with true parameters and synthetic observations

### Model Fitting

**`main.R`**
- Orchestrates full simulation pipeline
- Loops over simulation conditions
- Fits Stan models (DSEM_VAR1_RI and DSEM_NULL)
- Saves posterior samples and convergence diagnostics
- Integrates with HPC job arrays

### Fit Indices

**`scripts/fitindices_utils.R`**

Key functions:

---

## Stan Models

### `models/DSEM_VAR1_RI.stan` (True Model)

**Structure**:
- N subjects, Nt time points
- 2 latent factors (η) with AR(1) dynamics
- 6 observed indicators (y) per time point
- Random intercepts (u₀)
- Multivariate normal errors

**Parameters**:
- λ (factor loadings)
- Φ (AR coefficient matrix)
- Ψ (error covariances)

### `models/DSEM_NULL.stan` (Null Model)

**Structure**:
- Time-varying means
- Independent observations across time
- Used for computing fit indices via likelihood ratio

---

## Expected Results

Based on prior work (Saprinskiy, 2024; Fălcuşan et al.):



---

## Relevant Literature

### Foundation Studies

1. **Saprinskiy, N.** (2024). *Development and evaluation of dynamic model fit indices for dynamic structural equation models*. University of Tübingen.
   - Evaluated 4 Bayesian fit indices (BRMSEA, BMc, B̂, $B̂_{adj}$)
   - Demonstrated excellent power in detecting misspecifications
   - Introduced time-point-specific methods

2. **Fălcuşan, M., Zierahn, K. I., & Saprinskiy, N.** (2024). *Model fit assessment for dynamic structural equation models: Frequentist vs. Bayesian alternatives*. 
   - Compared traditional (RMSEA, CFI, TLI) vs. Bayesian indices
   - Found TLI and χ²/df most effective for complex DSEM
   - Highlighted need for improved time-dependent methods

### Key DSEM References

- **Asparouhov, T., Hamaker, E. L., & Muthén, B.** (2018). Dynamic structural equation models. *Structural Equation Modeling*, 25(3), 359-388.
- **Brandt, H.** (2023). Time-specific model fit in dynamic SEM. *Journal of Structural Equation Modeling*.

---


## Contributing

This is a Master thesis project. For questions or suggestions:

📧 **Contact**: josephliwanwang@gmail.com 
👨‍🎓 **Supervisor**: Prof. Dr. Holger Brandt (University of Tübingen)

---

## License

MIT License - See LICENSE file for details

---

## Citation

If you use this code in your research, please cite:

```bibtex
@mastersthesis{Wan2026dsem,
  author = {Joseph Li Wan Wang},
  title = {Time-Specific Fit Indices for Dynamic Structural Equation Models: A Simulation Study},
  school = {University of Tübingen},
  year = {2026},
  supervisor = {Prof. Dr. Holger Brandt}
}
```

---

## Status & Timeline

- ✅ Data generation: Complete
- ✅ Stan models: Complete  
- ✅ Fit indices functions: Complete
- ⏳ **Simulation runs**: In progress (Jan-Feb 2026)
- ⏳ Results analysis: Planned (Feb-Mar 2026)
- ⏳ Thesis writing: Planned (Mar-Apr 2026)

---

**Last Updated**: January 2026  
**Repository**: [github.com/your-username/dsem-fit-indices](https://github.com)

