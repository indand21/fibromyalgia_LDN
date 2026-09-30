# CHANGELOG

## 2026-09-29 — analysis started
Hypothesis-testing QSP analysis of low-dose naltrexone in fibromyalgia. Pharmacokinetics are fitted only to extracted published summary statistics with page/table locations (data/training/pk_summary.csv).

## 2026-09-29 — bottom-up binding failed PET test (RMSE 31.4). Calibrating Kd scale on Rabiner 2011 only (stage A2).

## 2026-09-29 - sampling effort raised for script 04
Script 04 first run failed the ess_bulk gate on log_tau (H1 280, H3 335; all rhat 1.00-1.02). Sampling effort raised to n_iter = 60000, n_burn = 30000 (4 chains) at the call site in scripts/04_fit_hypotheses.R. Computational only: no change to priors, likelihood, model structure, data or the gate.

## 2026-09-30 - six fixes following independent verification of the first complete run (run 3), then full rerun (run 4)
An independent verification of the first complete pipeline run confirmed the central conclusion (H1/H3 predict held-out null trials far worse than H0/H2) and found six defects. Fixed here, nothing else changed (no change to model structure, priors, likelihood, decision rules or any file in data/):
1. F1 TLR4 unit bug: pool_binding now converts concentration constants (Ki, Kd, IC50) to nM (uM x 1000, nM x 1, pM / 1000) and stops on an unrecognised unit. The curated 105.5 uM TLR4 IC50 was previously consumed as 105.5 nM-equivalent (1000x too potent). Data unchanged; consumer-side fix. Tests added.
2. F2 INNOVA estimand: scripts/06 now scores INNOVA against the unadjusted difference in change (arm change_nrs rows, day 90 = +0.31), with the SE borrowed from the adjusted row's 95% CI (approximation, recorded per row). The adjusted level difference is scored as a labelled SENSITIVITY row. The estimand is a column of test_predictions.csv.
3. F3 pseudo-replication: the PRIMARY held-out elpd uses one endpoint per trial (FINAL d84, INNOVA d90) (elpd_primary in test_elpd.csv); the all-days sum is kept only as a secondary, non-independent column. Scripts 07, 09 and 11 use the primary column.
4. F4 power curve: replaced by empirical rejection rates from the simulated diff distribution, two labelled variants (rej_rate_with_model_unc, rej_rate_sampling_only), with a definition column and a no-drug null check (power_null_check.csv).
5. F5 stage-consistent occupancy: figures and tables read the occupancy of the active binding stage (binding_stage.txt), computed from that stage's exposure cache; steady_state_occupancy_active.csv reports mean/peak/trough at 1.5, 4.5, 6 and 50 mg.
6. F6 freeze manifest: read_test_data now hashes each test file on first read (stage "T") and stops if it later changes. Tests added.
Because the TLR4 fix changes the H2 fit and stage B was already frozen, the previous outputs were archived to ../_archive/results_run3_<timestamp> and the pipeline was rerun cleanly (01-11). Full virtual_trials.csv is not committed (about 500 MB); virtual_trials_sample.csv (fixed-seed 1% sample) is committed instead.

## 2026-09-30 — bottom-up binding failed PET test (RMSE 31.4). Calibrating Kd scale on Rabiner 2011 only (stage A2).

## 2026-09-30 - run 4 follow-up (scripts 07, 09, 10, 11 rerun; 01-06 not rerun)
1. virtual_summary.csv and Fig 6A now use the same unadjusted change estimand as script 06 (INNOVA d90 = +0.31); the LMM-adjusted +0.49 is a labelled sensitivity column/dashed line. New columns obs_estimand and obs_diff_adjusted_level_sensitivity.
2. p_sig_two_sided (and its note column) removed from virtual_summary.csv; nothing reads it.
3. Script 03 no longer appends a duplicate CHANGELOG line when stage A2 triggers.
4. New output/results/test_elpd_contrasts.csv (pairwise primary elpd differences, SE, n, ratio, note that n = 2 is indicative only), also included in table3_heldout.csv.
5. No remaining reads of elpd_test_diff or power in scripts or tables.

## 2026-09-30 - run 5: correction of the binding-stage decision rule (F1) and three presentation fixes (F2-F4), then clean refit
An independent manuscript review (C1, C2, I2, I5, I8, I10) found that the pipeline carried the WORSE held-out predictor into every downstream result. Run 4 switched from stage A (bottom-up) to stage A2 (Kd scale fitted to the Rabiner rows) whenever stage-A RMSE exceeded 15 points, although A2 is worse on the rows it was not fitted to (5 held-out rows: A 22.69 vs A2 46.16; all 11 rows: A 31.44 vs A2 31.57). Changes (no change to model structure, priors, likelihood, convergence gate, freeze guards or thresholds):
1. F1 stage selection (scripts/03_test_pet.R, R/stage_select.R): A2 is still fitted on the Rabiner rows when stage A fails (RMSE > 15), but is ADOPTED for downstream use only if it improves RMSE on the non-calibration (held-out) rows relative to stage A on the same rows; otherwise stage A is kept. Decision and reason are in output/results/stage_decision.csv. pet_metrics.csv now also reports the matched held-out comparison (A and A2 on the same five rows) plus A2 all-rows and calibration-row metrics; every existing output is kept. Test added (tests/testthat/test-stage-select.R). NOTE: F1 is a POST-HOC correction made after the held-out results were seen, and is therefore exploratory in the sense the design spec defines. Because the stage changes the exposure cache that the hypothesis fits consume, everything downstream was refitted from a clean state (previous outputs archived to ../_archive/{results,figures}_run5pre_<timestamp>).
2. F2 data edit to a provenance file: data/module_screen.csv, TLR4 / MD2 row: criterion C1 now reads "No (not met at clinically reachable concentrations; ...)" and the reason states that the only available potency is a functional IC50 of about 105.5 uM in mouse microglia, some four orders of magnitude above unbound brain concentrations at 1.5 to 6 mg, and that H2 was nevertheless carried forward as a falsifiable hypothesis so that the prediction could be tested rather than assumed. Citations column and all other rows unchanged.
3. F3 Figure 2B (scripts/10_figures.R): now plots the predicted occupancy for the actual regimen and scan time of each Rabiner observation (same predict_pet() as script 03, now in R/pet.R) for both the adopted and the alternative stage, with the scored row-level predictions as hollow markers; panel C shows both stages likewise.
4. F4 (scripts/07_virtual_trials.R): the doubled-variability run is labelled as effectively doubling the placebo-response variability term (omega_p is about 2, omega_e only 1e-5 to 3e-3); virtual_summary.csv now records omega_p_used, omega_e_used and omega_label.
Also folded in (reporting only, no numerical effect): script 04 writes mcmc_diagnostics.csv (per-chain acceptance, per-hypothesis min ess_bulk and max rhat); script 02 records the provenance of each koff (koff_provenance.csv; koff_BN is an assumed 2 per h, now labelled "assumed" in table2_parameters.csv, which also gains a koff_BN row); the stale keq sensitivity note in Table 2 is replaced by "fixed at 2 per hour ... not varied"; sessionInfo.txt is written at the end of the run.
