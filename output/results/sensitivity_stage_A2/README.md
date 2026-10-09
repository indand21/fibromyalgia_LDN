# Sensitivity pass: binding stage A2 (not the adopted stage)

Written by scripts/run_all.R on 2026-10-08 19:43.
Scripts run, in a work copy of R/, scripts/ and data/: 04_fit_hypotheses.R, 05_calibrate_variability.R, 06_test_trials.R, 07_virtual_trials.R, 08_test_cytokines.R.
Script 07 ran in SUMMARY mode (V5_SCRIPT07_MODE=summary): virtual-trial summary and responder scores only;
no power curves and no SD-sensitivity calibration runs. Script 09 (Sobol) is not run in sensitivity passes.
Inputs taken from the primary pass: exposure_cache_A2.rds, binding_params_A2.rds, pet_metrics.csv; binding_stage.txt set to A2.
The work copy's own freeze manifest is kept here as freeze_manifest_sensitivity_pass.csv; the primary manifest was never touched.
virtual_trials.csv (full) is not kept; virtual_trials_sample.csv is.
