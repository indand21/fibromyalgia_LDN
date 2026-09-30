# QSP Fibromyalgia Low-Dose Naltrexone Model

This repository contains a hypothesis-testing quantitative systems pharmacology (QSP) model of low-dose naltrexone (LDN) in fibromyalgia. It uses extracted summary statistics from published studies, implements a Bayesian inference workflow with a frozen test-data protocol, and validating predictions against held-out PET, RCT and cytokine data.

## Running the Model

The complete model is executed with:

```r
Rscript scripts/run_all.R
```

This runs all 11 scripts in sequence from data fitting through sensitivity analysis and figure generation.

## Data Organization and Freeze Stages

Data is split into training and test sets, protected by a hash-based freeze manifest. Stage A (pharmacokinetics and bottom-up brain MOR binding) is frozen before testing against PET occupancy data. Stage B (clinical hypotheses) is frozen before testing against FINAL, INNOVA and cytokine outcomes. See `data/README.md` for detailed provenance of every dataset and `output/results/freeze_manifest.csv` for freeze timestamps.

## File Structure

- `R/`: All model code (utilities, PK, binding, disease dynamics, MCMC, prediction)
- `data/`: Training, parameter, prior and test datasets (all with `source`, `location`, `verified` provenance columns)
- `scripts/`: Numbered analysis scripts (01–11) plus `run_all.R`
- `tests/`: testthat unit tests for every module
- `output/results/`: Fitted parameters, predictions and numerical results
- `output/figures/`: Publication figures
