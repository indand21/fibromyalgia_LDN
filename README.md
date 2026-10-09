# QSP Fibromyalgia Low-Dose Naltrexone Model

This repository contains a hypothesis-testing quantitative systems pharmacology (QSP) model of low-dose naltrexone (LDN) in fibromyalgia. It uses extracted summary statistics from published studies, implements a Bayesian inference workflow with a frozen test-data protocol, and validating predictions against held-out PET, RCT and cytokine data.

## Running the Model

The complete model is executed with:

```r
Rscript scripts/run_all.R
```

This runs, in order:

1. the primary pass, scripts 01-12 (with 05b after 05) under the adopted binding stage in `output/results/binding_stage.txt`;
2. sensitivity passes for binding stages A2 and A (scripts 04, 05, 06, 07 in summary mode, and 08), each in a work copy with its own freeze manifest, written to `output/results/sensitivity_stage_<S>/`;
3. a summary across binding stages, `output/results/sensitivity_summary.csv`.

Each script runs in its own R process; logs are in `output/logs/`, per-script times in `output/logs/run_timing.txt` and the session information in `output/results/sessionInfo.txt`. The driver refuses to run over an existing freeze manifest or sensitivity output: archive the previous `output/results`, `output/figures` and `output/logs` (move them to `../_archive/<name>_<YYYYMMDD_HHMMSS>/`) before a clean re-run.

## Data Organization and Freeze Stages

Data is split into training and test sets, protected by a hash-based freeze manifest. Stage A (pharmacokinetics and bottom-up brain MOR binding) is frozen before testing against PET occupancy data. Stage B (clinical hypotheses) is frozen before testing against FINAL, INNOVA and cytokine outcomes. See `data/README.md` for detailed provenance of every dataset and `output/results/freeze_manifest.csv` for freeze timestamps.

## File Structure

- `R/`: All model code (utilities, PK, binding, disease dynamics, MCMC, prediction)
- `data/`: Training, parameter, prior and test datasets (all with `source`, `location`, `verified` provenance columns)
- `scripts/`: Numbered analysis scripts (01–11) plus `run_all.R`
- `tests/`: testthat unit tests for every module
- `output/results/`: Fitted parameters, predictions and numerical results
- `output/figures/`: Publication figures
