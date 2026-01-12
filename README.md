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

## Key Files

### Data Generation

**`scripts/dgp_sem.R`**

### Model Fitting

**`main.R`**

### Fit Indices

**`scripts/fitindices_utils.R`**

---

## Relevant Literature

### Foundation Studies

1. **Saprinskiy, N.** (2024). *Development and evaluation of dynamic model fit indices for dynamic structural equation models*. University of Tübingen.
   - Evaluated 4 Bayesian fit indices (BRMSEA, BMc, B̂, $B̂_{adj}$)
   - Demonstrated excellent power in detecting misspecifications
   - Introduced time-point-specific methods

2. **Fălcuşan, M., Zierahn, K. I., & Saprinskiy, N.** (2025). *Model fit assessment for dynamic structural equation models: Frequentist vs. Bayesian alternatives*. 
   - Compared traditional (RMSEA, CFI, TLI) vs. Bayesian indices
   - Found TLI and χ²/df most effective for complex DSEM
   - Highlighted need for improved time-dependent methods
---


## Contributing

This is a Master thesis project. For questions or suggestions:

📧 **Contact**: josephliwanwang@gmail.com 
👨‍🎓 **Supervisor**: Prof. Dr. Holger Brandt (University of Tübingen)

---

## License

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

**Last Updated**: January 2026  
**Repository**: [https://github.com/JosephLWW/Time-Point-Specific-Fit-Indices-for-DSEM/](https://github.com)

