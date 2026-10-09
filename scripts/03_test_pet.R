options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
pk <- readRDS("output/results/pk_fit.rds"); bpA <- readRDS("output/results/binding_params_A.rds")
pet <- read_test_data("data/test/pet_occupancy.csv", "A")
ed50_tab <- read_test_data("data/test/pet_ed50.csv", "A")   # Rabiner 2011 published ED50; hashed under stage T on first read

NORM_REF <- PET_NORM_REF   # predict_pet() and the reference conventions live in R/pet.R (shared with script 10)
if (!"reference" %in% names(pet)) stop("pet_occupancy.csv lacks the `reference` column")
if (!all(pet$reference %in% c("absolute", NORM_REF))) stop("unknown `reference` value in pet_occupancy.csv")
ed50_row <- ed50_tab[ed50_tab$endpoint == "MOR_occupancy" & ed50_tab$quantity == "ED50", , drop = FALSE]
if (nrow(ed50_row) != 1 || ed50_row$unit != "mg") stop("Rabiner ED50 test file: expected one MOR_occupancy ED50 row in mg")
ed50_obs <- c(value = ed50_row$value, lo95 = ed50_row$lo95, hi95 = ed50_row$hi95)
ED50_TIME_H <- 4   # single-dose scan time used for the Rabiner rows (pet_occupancy.csv time_basis: scans < 8 h post dose)
if (!all(pet$time_h[grepl("Rabiner", pet$study)] == ED50_TIME_H)) stop("Rabiner rows no longer share the 4 h scan time")
STD_DOSES <- c(1.5, 4.5, 6, 50); UNC_DOSES <- c(1.5, 4.5, 6)
KEQ_GRID <- c(0.5, 2, 8); KEQ_DEFAULT <- 2
N_DRAWS <- 500; UNC_SEED <- 20261001
N_WORKERS <- default_n_workers()

predict_pet_rows <- function(bp, rows) predict_pet(bp, rows, pk$par)
# RMSE/MAE/n_within_15 as before; wSSE (precision-weighted SSE, sd = reported sd or 10 points) added.
metrics <- function(rows, pred, stage, label = "all") data.frame(stage = stage, n = nrow(rows), RMSE = sqrt(mean((pred - rows$occupancy_pct)^2)),
  MAE = mean(abs(pred - rows$occupancy_pct)), n_within_15 = sum(abs(pred - rows$occupancy_pct) <= 15), rows = label,
  wSSE = pet_wsse_from_pred(pred, rows))
tag <- function(rows, pred, stage, role) { rows$predicted <- pred; rows$stage <- stage; rows$role <- role
  rows$sd_used <- pet_precision_sd(rows$sd); rows$weighted_residual <- (pred - rows$occupancy_pct) / rows$sd_used; rows }
is_rab <- grepl("Rabiner", pet$study); rab <- pet[is_rab, ]; rest <- pet[!is_rab, ]

# ---- (a) Stage A (bottom-up): a GENUINE held-out test against all 11 rows (no PET row is fitted) ----
predA <- predict_pet_rows(bpA, pet)
outA <- tag(pet, predA, "A (bottom-up)", "held-out (stage A fits no PET row)")
mA <- metrics(pet, predA, "A", "all rows (held-out: stage A fits no PET row)")
if (mA$RMSE > 15) {
  msg <- sprintf("bottom-up binding failed PET test (RMSE %.1f). Joint calibration of Kd and koff scales on all PET rows (stage J, adopted; post-hoc) with stage A2 (Kd scale on Rabiner 2011 only) as sensitivity.", mA$RMSE)
  cl <- if (file.exists("CHANGELOG.md")) readLines("CHANGELOG.md", warn = FALSE) else character(0)
  if (!any(grepl(msg, cl, fixed = TRUE)))
    cat(paste0("\n## ", Sys.Date(), " - ", msg, "\n"), file = "CHANGELOG.md", append = TRUE)
}

# ---- (b) Stage A2 (sensitivity): one Kd scale fitted to the six Rabiner rows only (unweighted, as before) ----
sseA2 <- function(ls) sum((predict_pet_rows(scale_binding(bpA, exp(ls), 1), rab) - rab$occupancy_pct)^2)
lsA2 <- stats::optimize(sseA2, c(-6, 6))$minimum
bpA2 <- modifyList(bpA, scale_binding(bpA, exp(lsA2), 1))
saveRDS(bpA2, "output/results/binding_params_A2.rds")
freeze_stage("A2", c("output/results/binding_params_A2.rds", "data/training/pk_summary.csv"))
predA2 <- predict_pet_rows(bpA2, pet)
outA2 <- rbind(tag(rab, predA2[is_rab], "A2 (calibration data)", "calibration (A2 Kd scale fitted to these rows)"),
               tag(rest, predA2[!is_rab], "A2 (Kd calibrated on Rabiner; held-out)", "held-out for A2"))

# ---- (c) Stage J (PRIMARY): shared Kd scale AND shared koff scale fitted jointly to ALL 11 rows, precision-weighted, ----
# ---- Lee rows on their normalised estimand; refitted at keq = 0.5, 2 and 8 per h (two starts each, in parallel)      ----
jobs <- expand.grid(keq = KEQ_GRID, koff0 = c(0.01, 0.1))
fit_job <- function(i, jobs, bpA, pet, pk_par, kd0) fit_binding_joint(bpA, pet, pk_par, keq = jobs$keq[i],
                                                                     start = c(log(kd0), log(jobs$koff0[i])))
environment(fit_job) <- globalenv()
all_fits <- with_cluster(min(nrow(jobs), N_WORKERS), function(cl)
  parallel::parLapply(cl, seq_len(nrow(jobs)), fit_job, jobs = jobs, bpA = bpA, pet = pet, pk_par = pk$par, kd0 = bpA2$kd_scale),
  seed = UNC_SEED)
fits <- lapply(KEQ_GRID, function(k) { f <- all_fits[jobs$keq == k]; f[[which.min(vapply(f, `[[`, 0, "wsse"))]] })
names(fits) <- as.character(KEQ_GRID)
keq_tab <- do.call(rbind, lapply(fits, function(f) {
  occ <- ss_mean_occupancy(UNC_DOSES, pk$par, f$bp)
  data.frame(keq = f$keq, wsse = f$wsse, rmse = f$rmse, kd_scale = f$kd_scale, log_kd_scale_se = f$se_log[[1]],
             koff_scale = f$koff_scale, log_koff_scale_se = f$se_log[[2]], hessian_pd = f$hessian_pd,
             Kd_NTX_nM = f$bp$Kd_NTX, Kd_BN_nM = f$bp$Kd_BN, koff_NTX_per_h = f$bp$koff_NTX, koff_BN_per_h = f$bp$koff_BN,
             mean_ss_occupancy_1.5mg = occ[["dose_1.5"]], mean_ss_occupancy_4.5mg = occ[["dose_4.5"]], mean_ss_occupancy_6mg = occ[["dose_6"]],
             n_starts = sum(jobs$keq == f$keq),
             start_wsse_spread = diff(range(vapply(all_fits[jobs$keq == f$keq], `[[`, 0, "wsse")))) }))
keq_choice <- select_keq(keq_tab, default = KEQ_DEFAULT)
keq_tab$adopted <- keq_tab$keq == keq_choice$keq
write.csv(keq_tab, "output/results/binding_J_keq_sensitivity.csv", row.names = FALSE)
fitJ <- fits[[as.character(keq_choice$keq)]]
bpJ <- modifyList(bpA, fitJ$bp)
bpJ$stage <- "J"; bpJ$keq_reason <- keq_choice$reason
bpJ$log_scales <- fitJ$log_scales; bpJ$log_scales_cov <- fitJ$cov; bpJ$log_scales_se <- fitJ$se_log
bpJ$Kd_NTX_source <- paste0("stage J: in vitro pooled Ki x kd_scale ", signif(fitJ$kd_scale, 4), " (fitted to all 11 PET rows); in vitro: ", bpA$Kd_NTX_source)
bpJ$Kd_BN_source <- paste0("stage J: in vitro Ki x kd_scale ", signif(fitJ$kd_scale, 4), " (shared with NTX); in vitro: ", bpA$Kd_BN_source)
bpJ$koff_NTX_source <- paste0("stage J: ", bpA$koff_NTX, " per h x koff_scale ", signif(fitJ$koff_scale, 4), " (fitted to all 11 PET rows); base: ", bpA$koff_NTX_source)
bpJ$koff_BN_source <- paste0("stage J: ", bpA$koff_BN, " per h x koff_scale ", signif(fitJ$koff_scale, 4), " (shared with NTX); base: ", bpA$koff_BN_source)
bpJ$koff_source <- "stage J: koff scale fitted jointly with the Kd scale to all 11 PET rows (post-hoc, exploratory)"
saveRDS(bpJ, "output/results/binding_params_J.rds")
saveRDS(list(fits = fits, all_fits = all_fits, jobs = jobs, keq_tab = keq_tab, keq_choice = keq_choice,
             pet_default_sd = PET_DEFAULT_SD, lee_tail_note = LEE_TAIL_NOTE, base = bpA[c("Kd_NTX", "Kd_BN", "koff_NTX", "koff_BN")]),
        "output/results/binding_fit_J.rds")
outJ <- tag(pet, fitJ$pred, "J (joint Kd + koff calibration; all rows are calibration data)", "calibration (stage J fitted to all 11 rows)")

# ---- metrics (existing rows kept; J rows and the wSSE column added) ----
m <- rbind(mA,
           metrics(rest, predA2[!is_rab], "A2 held-out", "held-out rows only (n = 5)"),
           metrics(rest, predA[!is_rab], "A held-out (matched)", "same held-out rows as A2 held-out"),
           metrics(pet, predA2, "A2 all", "all rows (6 calibration + 5 held-out)"),
           metrics(rab, predA2[is_rab], "A2 calibration", "calibration rows (Rabiner; fitted)"),
           metrics(rab, predA[is_rab], "A calibration", "calibration rows (Rabiner; not fitted under A)"),
           metrics(pet, fitJ$pred, "J all", "all rows (all are calibration data for J)"),
           metrics(rab, fitJ$pred[is_rab], "J Rabiner", "Rabiner rows (calibration)"),
           metrics(rest, fitJ$pred[!is_rab], "J Lee + Weerts", "Lee 1988 and Weerts 2008 rows (calibration)"))

# ---- (g) implied single-dose occupancy ED50 at 4 h per stage vs the Rabiner 2011 published ED50 ----
bps <- list(A = bpA, A2 = bpA2, J = bpJ)
ed50_mod <- vapply(bps, occupancy_ed50, 0, pk_par = pk$par, time_h = ED50_TIME_H)
ed50_out <- data.frame(stage = names(bps), model_ed50_mg = unname(ed50_mod), time_h = ED50_TIME_H,
  rabiner_ed50_mg = ed50_obs[["value"]], rabiner_lo95 = ed50_obs[["lo95"]], rabiner_hi95 = ed50_obs[["hi95"]],
  ratio_model_to_rabiner = unname(ed50_mod) / ed50_obs[["value"]],
  within_rabiner_ci = unname(ed50_mod) >= ed50_obs[["lo95"]] & unname(ed50_mod) <= ed50_obs[["hi95"]],
  comparison_independent = names(bps) == "A",
  note = c("independent: stage A fits no PET row",
           "NOT independent: the Rabiner ED50 summarises the Rabiner scans used to calibrate A2's Kd scale",
           "NOT independent: the Rabiner ED50 summarises the Rabiner scans used to calibrate stage J"))
write.csv(ed50_out, "output/results/occupancy_ed50_by_stage.csv", row.names = FALSE)

# ---- (e) exposure caches for A2 and J on the same dose grid as A (the A cache comes from script 02) ----
dose_grid <- as.numeric(names(readRDS("output/results/exposure_cache_A.rds")))
saveRDS(exposure_cache(dose_grid, pk$par, bpA2, max_doses = EXPOSURE_MAX_DOSES), "output/results/exposure_cache_A2.rds")
saveRDS(exposure_cache(dose_grid, pk$par, bpJ, max_doses = EXPOSURE_MAX_DOSES), "output/results/exposure_cache_J.rds")

# ---- (d) decision and adopted stage ----
dec <- binding_stage_decision(
  A = list(rmse = mA$RMSE, wsse = mA$wSSE, n = nrow(pet), ed50 = ed50_mod[["A"]]),
  A2 = list(kd_scale = bpA2$kd_scale, rmse_calib = m$RMSE[m$stage == "A2 calibration"], rmse_all = m$RMSE[m$stage == "A2 all"],
            wsse = m$wSSE[m$stage == "A2 all"], ed50 = ed50_mod[["A2"]]),
  J = list(kd_scale = fitJ$kd_scale, kd_scale_se_log = fitJ$se_log[[1]], koff_scale = fitJ$koff_scale,
           koff_scale_se_log = fitJ$se_log[[2]], wsse = fitJ$wsse, rmse = fitJ$rmse, ed50 = ed50_mod[["J"]], hessian_pd = fitJ$hessian_pd),
  keq_tab = keq_tab, keq_choice = keq_choice, ed50_obs = ed50_obs)
stage <- dec$stage_adopted
write.csv(dec, "output/results/stage_decision.csv", row.names = FALSE)
write.csv(rbind(outA, outA2, outJ), "output/results/pet_test.csv", row.names = FALSE)
write.csv(m, "output/results/pet_metrics.csv", row.names = FALSE)
writeLines(stage, "output/results/binding_stage.txt")

# ---- (f) steady-state occupancy per stage, and for the adopted stage ----
by_stage <- do.call(rbind, lapply(c("A", "A2", "J"), function(sg) {
  o <- steady_state_occupancy(readRDS(sprintf("output/results/exposure_cache_%s.rds", sg))); o$stage <- sg
  o[o$dose_mg %in% STD_DOSES, c("stage", "dose_mg", "mean_RO", "peak_RO", "trough_RO")] }))
by_stage$adopted <- by_stage$stage == stage
write.csv(by_stage, "output/results/steady_state_occupancy_by_stage.csv", row.names = FALSE)
write.csv(by_stage[by_stage$stage %in% c("A", "A2"), c("stage", "dose_mg", "mean_RO", "peak_RO", "trough_RO")],
          "output/results/steady_state_occupancy_both_stages.csv", row.names = FALSE)   # kept for existing readers (A vs A2)
for (sg in c("A2", "J")) write.csv(steady_state_occupancy(readRDS(sprintf("output/results/exposure_cache_%s.rds", sg))),
                                   sprintf("output/results/steady_state_occupancy_%s.csv", sg), row.names = FALSE)
occ <- active_occupancy("output/results")
write.csv(occ[occ$dose_mg %in% STD_DOSES, c("stage", "dose_mg", "mean_RO", "peak_RO", "trough_RO")],
          "output/results/steady_state_occupancy_active.csv", row.names = FALSE)
# DOR/KOR descriptive equilibrium occupancy for the active stage (Kd values scaled by the stage's Kd scale; 1 for stage A)
bc <- read_dataset("data/params/binding_constants.csv"); acache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
kd_sc <- if (!is.null(bps[[stage]]$kd_scale)) bps[[stage]]$kd_scale else 1
ki <- function(an, tg) kd_sc * pool_binding(bc, an, tg, "Ki")$value
dk <- do.call(rbind, lapply(c(4.5, 50), function(d) { pr <- acache[[as.character(d)]]
  do.call(rbind, lapply(c("DOR", "KOR"), function(tg) data.frame(stage = stage, dose_mg = d, target = tg,
    mean_occupancy = equilibrium_occupancy(mean(pr$Cu_NTX), mean(pr$Cu_BN), ki("NTX", tg), ki("BN", tg))))) }))
write.csv(dk, "output/results/dor_kor_occupancy_active.csv", row.names = FALSE)

# ---- occupancy uncertainty for stage J (review C4): Laplace draws on the log scales, optionally joint with PK draws ----
pk_pd <- isTRUE(pk$hessian_summary$pd)
unc_job <- function(i, dJ, dPK, pk, bp0, keq, doses) {
  p <- if (is.null(dPK)) pk$par else pk_from_theta(dPK[i, ], pk$par, pk$pars)
  v <- ss_mean_occupancy(doses, p, scale_binding(bp0, exp(dJ[i, 1]), exp(dJ[i, 2]), keq))
  c(v, n_warn = attr(v, "n_warn"))
}
environment(unc_job) <- globalenv()
unc_note <- character(0)
if (fitJ$hessian_pd) {
  dJ <- laplace_draws(fitJ$log_scales, fitJ$cov, N_DRAWS, seed = UNC_SEED)
  dPK <- if (pk_pd) laplace_draws(pk$theta, pk$hessian_summary$cov, N_DRAWS, seed = UNC_SEED + 1) else NULL
  if (!pk_pd) unc_note <- "PK Hessian not positive definite: only binding (stage J) uncertainty is propagated."
  runs <- list(binding_only = NULL, binding_and_pk = dPK)
  runs <- runs[!vapply(names(runs), function(nm) nm == "binding_and_pk" && is.null(dPK), TRUE)]
  res_unc <- with_cluster(N_WORKERS, function(cl) lapply(runs, function(dpk)
    do.call(rbind, parallel::parLapply(cl, seq_len(N_DRAWS), unc_job, dJ = dJ, dPK = dpk, pk = pk, bp0 = bpA,
                                       keq = fitJ$keq, doses = UNC_DOSES))), seed = UNC_SEED)
  point <- ss_mean_occupancy(UNC_DOSES, pk$par, bpJ)
  unc <- do.call(rbind, lapply(names(res_unc), function(nm) { M <- res_unc[[nm]]
    do.call(rbind, lapply(seq_along(UNC_DOSES), function(j) { q <- summarise_draws(M[, j])
      data.frame(stage = "J", dose_mg = UNC_DOSES[j], uncertainty_source = nm, mean_RO_point = point[[j]],
                 mean_RO_median = q[["median"]], mean_RO_lo95 = q[["lo95"]], mean_RO_hi95 = q[["hi95"]],
                 n_draws = q[["n"]], n_draws_requested = nrow(M), n_draws_not_steady = sum(M[, "n_warn"] > 0), seed = UNC_SEED,
                 method = if (nm == "binding_only") "Laplace (multivariate normal) on log kd_scale, log koff_scale; PK at its point estimate"
                          else "binding Laplace draws x independent PK Laplace draws (transformed PK parameters); APPROXIMATE: it ignores the PK-binding correlation, which widens the interval but also shifts it, so it is not a strict superset of the binding-only interval",
                 keq = fitJ$keq) })) }))
} else {
  unc_note <- "Stage-J Hessian not positive definite: Laplace draws not possible; no uncertainty band computed."
  unc <- data.frame(stage = "J", dose_mg = UNC_DOSES, uncertainty_source = NA, mean_RO_point = NA, mean_RO_median = NA,
                    mean_RO_lo95 = NA, mean_RO_hi95 = NA, n_draws = 0, n_draws_requested = N_DRAWS, n_draws_not_steady = NA,
                    seed = UNC_SEED, method = unc_note, keq = fitJ$keq)
}
unc$note <- paste(unc_note, collapse = " ")
# keq is fixed at the adopted value by the 3.84 rule, but it moves occupancy materially: its range over the keq grid is
# written as a SEPARATE, labelled structural sensitivity (not part of the Laplace band; report it alongside the band).
keq_rows <- do.call(rbind, lapply(seq_along(UNC_DOSES), function(j) {
  v <- keq_tab[[sprintf("mean_ss_occupancy_%smg", UNC_DOSES[j])]]
  data.frame(stage = "J", dose_mg = UNC_DOSES[j], uncertainty_source = "keq_sensitivity", mean_RO_point = v[keq_tab$keq == fitJ$keq],
             mean_RO_median = NA_real_, mean_RO_lo95 = min(v), mean_RO_hi95 = max(v), n_draws = length(v),
             n_draws_requested = length(v), n_draws_not_steady = NA,
             seed = NA, method = sprintf("STRUCTURAL SENSITIVITY, not a probability band: min and max of mean steady-state occupancy over stage-J refits at keq = %s per h (values %s); keq = %g kept by the 3.84 rule; report alongside the Laplace band (details in binding_J_keq_sensitivity.csv)",
                                         paste(keq_tab$keq, collapse = " / "), paste(signif(v, 3), collapse = " / "), fitJ$keq),
             keq = NA_real_, note = "keq range") }))
unc <- rbind(unc, keq_rows)
write.csv(unc, "output/results/occupancy_uncertainty_J.csv", row.names = FALSE)

# Weerts 2008 is the only daily PET row and the most heavily weighted (sd 4.9), and its scan time is an ASSUMPTION: the
# abstract gives the scan day but not the time after that day's dose (time_basis column). Prediction across plausible scan
# times at the curated dose count, so the reader can see how much of the row's influence rests on that assumption.
# R4: primary row is n_doses=3 (doses days 15/16/17, scan on day 18 at 24 h trough; R4 consistency fix).
#     The 4-dose sensitivity (doses days 15-18, scan 2-12 h after day-18 dose) is generated here as a separate block,
#     per R4-instructions.md section 3: "add 4 doses at 2-12 h as a sensitivity".
w_row <- pet[pet$study == "Weerts 2008", , drop = FALSE]
if (nrow(w_row) == 1) {
  w_grid <- do.call(rbind, lapply(c(1, 2, 4, 8, 12, 24), function(h) { r <- w_row; r$time_h <- h; r }))
  w_sens <- data.frame(study = "Weerts 2008", dose_mg = w_row$dose_mg, n_doses = w_row$n_doses, time_h = w_grid$time_h,
                       observed_pct = w_row$occupancy_pct, sd = w_row$sd,
                       predicted_J = predict_pet(bpJ, w_grid, pk$par),
                       assumed_in_data = w_grid$time_h == w_row$time_h,
                       note = sprintf("primary: n_doses=%d (R4 consistency fix); scan time is assumed (time_basis: %s)",
                                      w_row$n_doses, w_row$time_basis))
  # Sensitivity block: 4 doses (days 15-18), post-dose scan at 2-12 h on day 18. This was the pre-R4 assumption;
  # carried forward as an explicitly labelled sensitivity, not the primary.
  w_row4 <- w_row; w_row4$n_doses <- 4L
  w_grid4 <- do.call(rbind, lapply(c(2, 4, 8, 12, 24), function(h) { r <- w_row4; r$time_h <- h; r }))
  w_sens4 <- data.frame(study = "Weerts 2008", dose_mg = w_row4$dose_mg, n_doses = 4L, time_h = w_grid4$time_h,
                        observed_pct = w_row$occupancy_pct, sd = w_row$sd,
                        predicted_J = predict_pet(bpJ, w_grid4, pk$par),
                        assumed_in_data = FALSE,
                        note = ifelse(w_grid4$time_h == 24,
                          "SENSITIVITY: n_doses=4 with the 24 h scan (the pre-R4 curated row, internally inconsistent; for comparison with the primary 3-dose row)",
                          "SENSITIVITY: n_doses=4 (doses days 15-18, scan post-dose day 18 at 2-12 h)"))
  write.csv(rbind(w_sens, w_sens4), "output/results/pet_weerts_scan_time_sensitivity.csv", row.names = FALSE)
  # R4 fix round 1 (review M6): the 3- versus 4-dose comparison is computed here, at the adopted stage-J fit and PK.
  w_dose_note <- sprintf("Weerts 2008 dose count (stage J): primary 3 doses with the assumed 24 h trough scan predicts %.1f%%; the pre-R4 4 doses at 24 h predicts %.1f%%; 4 doses with a post-dose scan at 2-12 h predicts %.1f-%.1f%% (observed %.1f +/- %.1f%%; pet_weerts_scan_time_sensitivity.csv).",
                         w_sens$predicted_J[w_sens$assumed_in_data][1], w_sens4$predicted_J[w_sens4$time_h == 24],
                         min(w_sens4$predicted_J[w_sens4$time_h < 24]), max(w_sens4$predicted_J[w_sens4$time_h < 24]), w_row$occupancy_pct, w_row$sd)
} else w_dose_note <- "Weerts 2008 row not found: no dose-count comparison."

# ---- stage-J sensitivity to the PK: fp at 100 mg fixed at 0.90/0.99 (primary 0.95); terminal-phase PK; ----
# R4 fix round 1: every fp here is a presystemic CONVERSION fraction (to 6-beta-naltrexol), never an escape fraction.
# fp is its value at the 100 mg reference dose; under an adopted dose-dependent structure fp0 is the conversion fraction
# as dose -> 0 (or floss0 the parent-only presystemic loss as dose -> 0), anchored so that 100 mg is unchanged. Only when
# such a structure is ADOPTED does script 01 write its D50 = 10 and 100 mg refits (pk_fit_sens_D50_10.rds,
# pk_fit_sens_D50_100.rds), which the glob below then picks up; under the linear model they do not exist. Each file name
# is pk_fit_sens_<key>.rds with key from PK_SENS_KEYS (R/pk.R): the `key` column, not the label, identifies a row.
# The terminal fit is the one direct test of slow MOR dissociation vs long terminal PK phase. Refitted at the adopted keq.
sens_files <- sort(list.files("output/results", pattern = "^pk_fit_sens_.*[.]rds$", full.names = TRUE))
if (length(sens_files)) {
  sens_pk <- lapply(sens_files, readRDS); names(sens_pk) <- sub("^pk_fit_sens_(.*)[.]rds$", "\\1", basename(sens_files))
  sjobs <- expand.grid(k = seq_along(sens_pk), koff0 = c(0.01, 0.1))
  sens_job <- function(i, sjobs, sens_pk, bpA, pet, keq, kd0) fit_binding_joint(bpA, pet, sens_pk[[sjobs$k[i]]]$par, keq = keq,
                                                                             start = c(log(kd0), log(sjobs$koff0[i])))
  environment(sens_job) <- globalenv()
  sfits_all <- with_cluster(min(nrow(sjobs), N_WORKERS), function(cl)
    parallel::parLapply(cl, seq_len(nrow(sjobs)), sens_job, sjobs = sjobs, sens_pk = sens_pk, bpA = bpA, pet = pet,
                        keq = fitJ$keq, kd0 = fitJ$kd_scale), seed = UNC_SEED)
  # `key` (stable; review C1) identifies the row for scripts 10 and 11; `analysis` is a human-readable label only.
  # fp = conversion fraction at 100 mg; fp0 / floss0 / fp_D50_mg are NA / NA / Inf under the linear model.
  sens_row <- function(key, label, pkf, f) { occ <- ss_mean_occupancy(UNC_DOSES, pkf$par, f$bp)
    data.frame(key = key, analysis = label, pk_structure = pkf$structure, pk_dose_model = pk_dose_model_of(pkf$par), fp = pkf$par$fp,
               fp0 = if (!is.null(pkf$par$fp0)) pkf$par$fp0 else NA_real_,
               floss0 = if (!is.null(pkf$par$floss0)) pkf$par$floss0 else NA_real_,
               fp_D50_mg = if (pk_dose_model_of(pkf$par) != "linear") pkf$par$fp_D50 else Inf,
               bioavailability_4.5mg = pk_bioavailability(4.5, pkf$par),
               parent_t_half_terminal_h = pk_parent_half_lives(pkf$par)[["terminal"]], keq = f$keq,
               kd_scale = f$kd_scale, log_kd_scale_se = f$se_log[[1]], koff_scale = f$koff_scale, log_koff_scale_se = f$se_log[[2]],
               koff_NTX_per_h = f$bp$koff_NTX, wsse = f$wsse, rmse = f$rmse, hessian_pd = f$hessian_pd,
               mean_ss_occupancy_1.5mg = occ[["dose_1.5"]], mean_ss_occupancy_4.5mg = occ[["dose_4.5"]], mean_ss_occupancy_6mg = occ[["dose_6"]]) }
  occ_sens <- rbind(sens_row(OCC_SENS_KEYS[["primary"]], "primary", pk, fitJ), do.call(rbind, lapply(seq_along(sens_pk), function(k) {
    f <- sfits_all[sjobs$k == k]; f <- f[[which.min(vapply(f, `[[`, 0, "wsse"))]]
    sens_row(names(sens_pk)[k], sens_pk[[k]]$sensitivity, sens_pk[[k]], f) })))
  # STRUCTURAL sensitivity to the ASSUMED 10:1 base koff ratio. The shared-scale ruling keeps koff_BN / koff_NTX at its base
  # value (2 / 21 per h, the 2 being an assumption: no 6-BN koff is curated). Equalising the two base rates is the strongest
  # alternative; it is a single refit, reported as a structural sensitivity alongside keq.
  bp_eq <- binding_params(Kd_NTX = bpA$Kd_NTX, Kd_BN = bpA$Kd_BN, koff_NTX = bpA$koff_NTX, koff_BN = bpA$koff_NTX, keq = fitJ$keq)
  f_eq <- fit_binding_joint(bp_eq, pet, pk$par, keq = fitJ$keq, start = fitJ$log_scales)
  r_eq <- sens_row(OCC_SENS_KEYS[["koff_eq"]], sprintf("koff_BN base equal to koff_NTX (%.0f per h; primary assumes %.0f, a 10:1 ratio)",
                           bpA$koff_NTX, bpA$koff_BN), pk, f_eq)
  occ_sens <- rbind(occ_sens, r_eq)
  occ_sens$note <- ifelse(occ_sens$key == OCC_SENS_KEYS[["koff_eq"]],
    sprintf("The base 10:1 koff ratio is ASSUMED and never varied by the shared-scale calibration. Equalising it fits the PET data WORSE (weighted SSE %.1f against %.1f), so the data support the assumption, but occupancy at 4.5 mg moves from %.3f to %.3f (%+.1f pp; largest other refit shift %+.1f pp). Report alongside the keq range.",
            f_eq$wsse, fitJ$wsse, occ_sens$mean_ss_occupancy_4.5mg[1], r_eq$mean_ss_occupancy_4.5mg,
            100 * (r_eq$mean_ss_occupancy_4.5mg - occ_sens$mean_ss_occupancy_4.5mg[1]),
            { o <- 100 * (occ_sens$mean_ss_occupancy_4.5mg[!occ_sens$key %in% OCC_SENS_KEYS] - occ_sens$mean_ss_occupancy_4.5mg[1]); o[which.max(abs(o))] }), "")
  if (anyDuplicated(occ_sens$key)) stop("occupancy_sensitivity_J.csv: duplicated key")
  write.csv(occ_sens, "output/results/occupancy_sensitivity_J.csv", row.names = FALSE)
}

# Goodness of fit at the stage-J optimum, and how far the keq decision depends on the assumed default sd. The weighted SSE is
# a chi-squared statistic, so wSSE / (n - p) near 1 says the assumed sds are about right. Scaling the DEFAULT sd scales the
# contribution of every row that uses it, hence the sd at which the keq = 8 gain would reach the 3.84 threshold.
j_df <- nrow(pet) - 2L
chi2_red <- fitJ$wsse / j_df
keq_best <- keq_tab$keq[which.min(keq_tab$wsse)]
gain_keq <- keq_tab$wsse[keq_tab$keq == fitJ$keq] - min(keq_tab$wsse)
# The keq gain recomputed at other default sds, holding each keq fit's predictions fixed; the largest default sd at which it
# would still clear 3.84 is the point beyond which the adopted keq changes.
wsse_at <- function(f, sd_def) sum(((f$pred - pet$occupancy_pct) / pet_precision_sd(pet$sd, sd_def))^2)
sd_grid <- seq(2, PET_DEFAULT_SD, by = 0.05)
gain_at <- vapply(sd_grid, function(s) wsse_at(fits[[as.character(fitJ$keq)]], s) - min(vapply(fits, wsse_at, 0, sd_def = s)), 0)
sd_flip <- if (any(gain_at > stats::qchisq(0.95, 1))) max(sd_grid[gain_at > stats::qchisq(0.95, 1)]) else NA_real_
gof_note <- sprintf("Stage-J goodness of fit: weighted SSE %.2f on %d rows with 2 fitted scales, so reduced chi-squared %.2f. Being close to 1, it supports the documented default sd of %g percentage points for the rows that report none. The keq decision does, however, depend on that weight scale: the keq = %g gain over the adopted keq = %g is %.2f at a default sd of %g, and would reach the 3.84 threshold at a default sd of about %s points, below which keq = %g would be adopted instead. Report the keq range as a structural sensitivity (occupancy_uncertainty_J.csv).",
                    fitJ$wsse, nrow(pet), chi2_red, PET_DEFAULT_SD, keq_best, fitJ$keq, gain_keq, PET_DEFAULT_SD,
                    if (is.na(sd_flip)) "none in 2-10" else sprintf("%.1f", sd_flip), keq_best)

writeLines(c(sprintf("Stage-A (bottom-up) PET RMSE: %.2f percentage points (n = %d; MAE %.2f; %d within 15 points; weighted SSE %.1f). Stage A fits no PET row, so this is a genuine held-out test.",
                     mA$RMSE, mA$n, mA$MAE, mA$n_within_15, mA$wSSE),
             sprintf("Stage A2 (Kd scale %.3f fitted to the Rabiner rows only): calibration-row RMSE %.2f, all-row RMSE %.2f (sensitivity bound).",
                     bpA2$kd_scale, m$RMSE[m$stage == "A2 calibration"], m$RMSE[m$stage == "A2 all"]),
             sprintf("Stage J (adopted): kd_scale %.3f (SE of log %.3f), koff_scale %.5f (SE of log %.3f), keq %g per h; weighted SSE %.2f, RMSE %.2f on all 11 rows (all are calibration data).",
                     fitJ$kd_scale, fitJ$se_log[[1]], fitJ$koff_scale, fitJ$se_log[[2]], fitJ$keq, fitJ$wsse, fitJ$rmse),
             sprintf("Precision weights: sd = reported sd where finite, otherwise %g percentage points (Rabiner digitised single points; Lee 1988 SEM-only).", PET_DEFAULT_SD),
             sprintf("keq sensitivity (0.5, 2, 8 per h): %s", keq_choice$reason),
             sprintf("Implied single-dose occupancy ED50 at %g h: A %.2f mg, A2 %.2f mg, J %.2f mg; Rabiner 2011 published %.2f mg (%.2f-%.2f). Only the stage-A comparison is independent: for A2 and J the ED50 summarises Rabiner scans that were used for calibration.",
                     ED50_TIME_H, ed50_mod[["A"]], ed50_mod[["A2"]], ed50_mod[["J"]], ed50_obs[["value"]], ed50_obs[["lo95"]], ed50_obs[["hi95"]]),
             sprintf("Binding stage carried forward: %s. %s", stage, dec$reason),
             if (length(unc_note)) unc_note else "Occupancy uncertainty: binding-only and binding + PK Laplace bands in occupancy_uncertainty_J.csv; the keq range (0.5 to 8 per h) is a separate structural sensitivity in the same file (uncertainty_source keq_sensitivity) and must be reported alongside the band.",
             sprintf("PK sensitivity of stage J (fp at 100 mg fixed at 0.90 and 0.99 instead of 0.95; terminal-phase PK with the Verebey 24-72 h naltrexone values appended%s) and the assumed 10:1 base koff_BN/koff_NTX ratio: occupancy_sensitivity_J.csv (rows identified by the key column). PK dose model: %s.",
                     if (pk_dose_model_of(pk$par) != "linear") "; D50 of the adopted dose-dependent structure refitted at 10 and 100 mg (primary 30 mg)" else
                       "; the dose-dependent first-pass structures were not adopted, so their D50 fits are PK outputs only (pk_dose_model_comparison.csv)",
                     PK_DOSE_MODELS[[pk_dose_model_of(pk$par)]]$label),
             w_dose_note,
             gof_note,
             sprintf("Dose counts: each PET row is simulated with the number of doses in the n_doses column of pet_occupancy.csv (%s). Steady state: each exposure profile is extended until the extrapolated remaining gap in mean occupancy is at most %g (absolute); the dose counts used are %s.",
                     paste(sprintf("%s %d", substr(pet$study, 1, 12), pet$n_doses), collapse = "; "), EXPOSURE_SS_TOL,
                     paste(sprintf("%g mg: %d", by_stage$dose_mg[by_stage$stage == stage],
                                   vapply(readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))[as.character(by_stage$dose_mg[by_stage$stage == stage])],
                                          function(z) attr(z, "n_doses"), 0L)), collapse = "; ")),
             sprintf("Reference convention: %d row(s) 'absolute' (prediction = 100 * RO(t)); %d row(s) '%s' (prediction = 100 * RO(t) / RO(1 h), same dose and regimen); the convention used per row is in the `reference` column of pet_test.csv",
                     sum(pet$reference == "absolute"), sum(pet$reference == NORM_REF), NORM_REF),
             paste("Unbound fraction of 6-beta-naltrexol:", if (!is.null(pk$fu_BN_note)) pk$fu_BN_note else "see pk_fit_notes.txt"),
             paste("Lee 1988 tail, identifiability:", LEE_TAIL_NOTE)),
           "output/results/pet_test_note.txt")
freeze_stage("J", c("data/test/pet_occupancy.csv", "data/training/pk_summary.csv", "data/params/binding_constants.csv",
                    "output/results/pk_fit.rds", "output/results/binding_params_A.rds",
                    "output/results/binding_params_J.rds", "output/results/binding_fit_J.rds",
                    "data/training/pk_summary_sensitivity_terminal.csv", sens_files))
print(m); print(keq_tab); print(ed50_out[, 1:8]); print(unc[, 1:8])
