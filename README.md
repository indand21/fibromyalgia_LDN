# QSP Model of Low-Dose Naltrexone in Fibromyalgia

A quantitative systems pharmacology (QSP) model that mechanistically characterizes how
low-dose naltrexone (LDN, 4.5 mg/day) exerts therapeutic effects in fibromyalgia through
the balance of picomolar OGF-OGFr anti-inflammatory signaling, TRIF-IRF3-biased TLR4
antagonism, OPRM1 hormesis, neuroinflammatory feedback, and an expectation-mediated
placebo response. The model comprises 37 ordinary differential equations across eleven
biological modules and reproduces the clinically observed U-shaped dose-response.

This repository contains the model code, input data, analysis scripts, and simulation
outputs needed to reproduce the results. **Manuscript files are intentionally not included.**

## Repository structure

```
model/                 QSP model (mrgsolve specification)
  qsp_model.R
scripts/               Analysis and simulation pipeline (R)
  run_desolve.R          Main v4 model (deSolve) -> steady_state_v4 / full_simulation_v4
  pk_calibration.R       PK calibration to clinical datasets
  pk_calibrated_model.R  Calibrated PK model
  dose_comparison.R      Dose-response (LDN vs standard)
  dose_titration.R       Dose-titration simulation
  sensitivity_analysis.R Local (OAT) + global (LHS-PRCC) sensitivity
  monte_carlo_virtual_patients.R  Virtual population / Monte Carlo
  validate_clinical.R    Validation against published clinical data
  qsp_model_functions.R  Shared helper functions
  prepare_data.R         Network-pharmacology data preparation
  run_simulation.R       Simulation driver
parameters.R           Model parameters with literature sources
initial_conditions.R   Baseline fibromyalgia disease state
data/                  Network-pharmacology inputs (hub genes, PPI network, scores)
output/
  results/             Simulation result tables (CSV); v4 files are the current model
  figures/             Generated figures (PNG)
PROJECT_PLAN.md        Model design and module specification
SCIENTIFIC_IMPROVEMENTS.md  Notes on model refinements
```

## Requirements

- R (>= 4.5)
- R packages: `deSolve`, `dplyr`, `ggplot2`, `tidyr` (and `mrgsolve` for `model/qsp_model.R`)

```r
install.packages(c("deSolve", "dplyr", "ggplot2", "tidyr"))
```

## Reproducing the main results

The current model is v4 (with placebo response). From the `scripts/` directory:

```r
source("run_desolve.R")   # main QSP simulation -> ../output/results/*_v4.csv, ../output/figures/*_v4.png
```

Other analyses (PK calibration, dose-response, sensitivity, virtual patients, validation)
are in the correspondingly named scripts.

## Source of truth

The headline results reported for the model derive from the v4 outputs:
`output/results/steady_state_v4.csv`, `prcc_pain.csv`, and the calibrated PK files.

## License

No license file is included; please contact the authors regarding reuse.
