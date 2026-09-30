# Sensitivity analysis: binding stage A2 (NOT adopted)

Run 5 adopted stage A (bottom-up) because stage A2 (Kd scaled to the Rabiner PET rows) did not improve RMSE on the
held-out PET rows (see ../stage_decision.csv). This directory holds scripts 04, 05, 06, 07, 08 and 09 rerun with
binding_stage.txt forced to A2, i.e. with hypothesis fits, held-out predictions, cytokine predictions, virtual trials
and Sobol indices all computed from the A2 exposure cache. It exists so a reader can see how much each headline number
depends on the stage choice. Nothing in the parent directory was overwritten.

How it was made: the pipeline was copied to a scratch working tree after script 03, binding_stage.txt was set to A2, and
scripts 04, 05, 06, 08, 09, 07 were run unmodified there (n_iter = 60000, n_burn = 30000, convergence gate and freeze guards
intact; the stage-B freeze happened in the copy). Script 05 (variability calibration) was rerun because script 07 consumes
omegas.csv. The full virtual_trials.csv is not included (about 500 MB); virtual_trials_sample.csv is the 1 per cent sample.
Logs and per-script runtimes are in logs/. The DOR/KOR, table and figure scripts (10, 11) were not run on this pass.
