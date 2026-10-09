options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
t_start <- Sys.time()
obs <- read_dataset("data/training/pk_summary.csv")
bc <- read_dataset("data/params/binding_constants.csv")
start <- pk_default()
fu_NTX <- pool_binding(bc, "NTX", "plasma_protein", "fu")$value
fu_BN <- pool_binding(bc, "BN", "plasma_protein", "fu")$value
start$fu_NTX <- fu_NTX
start$fu_BN <- if (is.finite(fu_BN)) fu_BN else start$fu_NTX
fu_BN_note <- if (is.finite(fu_BN)) paste("fu_BN from binding_constants.csv:", pool_binding(bc, "BN", "plasma_protein", "fu")$source) else
  sprintf("fu_BN = %.2f is taken EQUAL TO THE PARENT's fu_NTX by fallback: no 6-beta-naltrexol plasma-protein-binding value was found (data/README.md 'Not found' list). It is an assumption, not a measurement, and it scales every unbound 6-beta-naltrexol concentration.", start$fu_NTX)
RES <- "output/results"

# fp (fraction of the absorbed dose converted presystemically to 6-beta-naltrexol) is FIXED (coordinator ruling,
# R1 fix round 3): with the 24 h rows the free fit ran to fp -> 1 on a ridge where only the parent's apparent values
# (CL, V1, Q, V2 divided by 1 - fp) are identified (smallest Hessian eigenvalue 0.037; WSS changes <= 0.29 for fp 0.9-0.999).
# R4 fix round 1: fp is the value AT THE REFERENCE DOSE (100 mg, Verebey); a dose-dependent structure is anchored there,
# so fp(100 mg) = 0.95 exactly whatever its D50 (R/pk.R, pk_dose_shape).
FP_FIXED <- 0.95
fp_note <- sprintf("fp fixed at %.2f at the reference dose of %g mg (the Verebey 1976 dose): the approved naltrexone labelling gives oral bioavailability of 5 to 40 per cent, i.e. a systemic parent fraction (1 - fp) of roughly 0.05 to 0.40; %.2f sits at the UPPER EDGE of that range (1 - fp = 0.05) and inside the region where the weighted SSE was flat in the R1 free fit (change <= 0.29 for fp 0.9 to 0.999). A dose-dependent structure is anchored at %g mg, so this value holds there exactly for every D50. Occupancy sensitivity at fp = 0.90 and 0.99 (at %g mg) is reported in occupancy_sensitivity_J.csv (fits pk_fit_sens_fp090.rds, pk_fit_sens_fp099.rds).",
                   FP_FIXED, FP_DREF_MG, FP_FIXED, FP_DREF_MG, FP_DREF_MG)
start$fp <- FP_FIXED
FREE_PARS <- setdiff(PK_FIT_PARS, "fp")
PAR_MAP <- function(cl) function(X, FUN) parallel::parLapply(cl, X, FUN)   # fits are deterministic: parallelism is neutral
is_ver <- grepl("^Verebey", obs$study); is_bes <- grepl("^Bested", obs$study)

# ---- 1. Compartments, under the linear dose model ----
# Review C3 (R1): Tmax is IN the objective (CV floor 0.25, because observed Tmax values sit on the sampling grid), since it
# is the only datum that identifies ka. Two-compartment only if WSS drops by > chi-square(2) 95%. Weights: R/pk.R
# pk_obs_cv / pk_obs_weight (R4 fix round 1: reported between-subject SD, else a documented assumed CV, else the floor;
# the Bested day-1 profile rows carry 1/4 each).
sel_lin <- with_cluster(2L, function(cl) pk_select_structure(obs, start, free = FREE_PARS, map_fits = PAR_MAP(cl)))
if (sel_lin$fit$convergence != 0) stop("PK fit did not converge (optim convergence code ", sel_lin$fit$convergence, ")")
# ---- 2. Identifiability gate on the linear fit (ka floor; Hessian positive definite), R/pk.R ----
f0 <- sel_lin$fit; f0$dose_model <- "linear"
hs_initial_lin <- pk_hessian_summary(f0)
g_lin <- pk_identifiability_gate(f0, obs)
fit_lin <- g_lin$fit; fit_lin$structure <- sel_lin$fit$structure; fit_lin$dose_model <- "linear"
hs_lin <- pk_hessian_summary(fit_lin)

# ---- 3. Dose dependence: linear against two one-extra-parameter structures (R4 fix round 1, review I3, I6, I7) ----
# Both alternatives are fitted on the SELECTED compartment structure with EXACTLY the linear fit's free set plus one
# parameter (fp0 for fp_dose, floss0 for parent_loss), from the linear optimum, with D50 fixed (primary 30 mg; 10 and
# 100 mg as sensitivity, also anchored at 100 mg, so each D50 changes only the shape below 100 mg). The primary fits
# (D50 30 mg) use two starts and keep the lower WSS; the sensitivity fits one. The comparison table and the rule use these
# same fits and the same objective.
# PRE-STATED ADOPTION RULE (written before any of these fits was run; R/pk.R pk_adoption_rule, applied mechanically):
#   adopt a one-extra-parameter structure only if dWSS = WSS(linear) - WSS(alternative) > 3.84 (chi-square, 1 df),
#   AND its Hessian is positive definite, AND no Verebey 1976 row is predicted worse than by the linear model by more
#   than a factor of 1.5 (ratio of fold errors); if both alternatives qualify, the larger dWSS wins; otherwise linear.
FP_D50_PRIMARY <- 30
FP_D50_ALL <- c(10, 30, 100)
ALT <- c("fp_dose", "parent_loss")
ALT_STARTS <- list(fp_dose = c(0.97, 0.995), parent_loss = c(0.3, 0.8))   # primary fits; sensitivity fits use the second
BESTED_CV_SENS <- c(0.5, 1.2)        # sensitivity of the assumed day-21 CV (primary 0.8)
obs_cv <- function(cv) { o <- obs; o$cv_assumed[is_bes & o$regimen == "daily"] <- cv; o }
alt_jobs <- function(model, d50, o, base_par, base_pars, tag, primary = d50 == FP_D50_PRIMARY && identical(o, obs))
  lapply(if (primary) ALT_STARTS[[model]] else ALT_STARTS[[model]][2], function(s0) list(id = tag, model = model, d50 = d50, obs = o, pars = c(base_pars, PK_DOSE_MODELS[[model]]$par),
         start = modifyList(base_par, stats::setNames(list(s0, d50), c(PK_DOSE_MODELS[[model]]$par, "fp_D50")))))
jobs <- c(
  unlist(lapply(ALT, function(m) unlist(lapply(FP_D50_ALL, function(d) alt_jobs(m, d, obs, fit_lin$par, fit_lin$pars, sprintf("%s|D50=%g|cv=primary", m, d))), recursive = FALSE)), recursive = FALSE),
  unlist(lapply(BESTED_CV_SENS, function(cv) c(
    list(list(id = sprintf("linear|cv=%g", cv), model = "linear", d50 = Inf, obs = obs_cv(cv), pars = fit_lin$pars, start = fit_lin$par)),
    unlist(lapply(ALT, function(m) alt_jobs(m, FP_D50_PRIMARY, obs_cv(cv), fit_lin$par, fit_lin$pars, sprintf("%s|D50=%g|cv=%g", m, FP_D50_PRIMARY, cv))), recursive = FALSE))),
    recursive = FALSE))
alt_raw <- with_cluster(min(length(jobs), default_n_workers()), function(cl) parallel::parLapply(cl, jobs, pk_fit_job))
ids <- vapply(jobs, `[[`, "", "id")
best_of <- function(id) { k <- which(ids == id); f <- alt_raw[k][[which.min(vapply(alt_raw[k], `[[`, 0, "value"))]]
  j <- jobs[[k[1]]]; f$dose_model <- j$model; f$structure <- fit_lin$structure; f$fp_D50_fixed <- j$d50; f }
alt_fit <- function(m, d50 = FP_D50_PRIMARY, cv = "primary") best_of(sprintf("%s|D50=%g|cv=%s", m, d50, cv))
# The candidate table for any set of observations: linear reference and the two alternatives (same objective, same fits).
LOGIT_BOUND <- logit(0.9999)          # an extra parameter beyond 0.9999 (or below 0.0001) is reported as a BOUNDARY fit
cand_row <- function(f, lin, o) {
  m <- f$dose_model; xp <- PK_DOSE_MODELS[[m]]$par; hs <- pk_hessian_summary(f)
  pr <- pk_predict_stats(f$par, o); pl <- pk_predict_stats(lin$par, o); rw <- pk_row_wss(pr, o)
  data.frame(dose_model = m, label = PK_DOSE_MODELS[[m]]$label, structure = f$structure,
             fp_D50_mg = if (m == "linear") NA_real_ else f$par$fp_D50, n_parameters = length(f$pars),
             extra_parameter = if (is.null(xp)) NA_character_ else xp,
             extra_estimate = if (is.null(xp)) NA_real_ else f$par[[xp]],
             extra_se_logit = if (is.null(xp)) NA_real_ else unname(hs$se[xp]),
             at_boundary = if (is.null(xp)) FALSE else abs(f$theta[[xp]]) > LOGIT_BOUND,
             wss = f$value, delta_wss = lin$value - f$value, wss_verebey = sum(rw[grepl("^Verebey", o$study)]),
             wss_bested = sum(rw[grepl("^Bested", o$study)]), hessian_pd = hs$pd, min_eigenvalue = min(hs$eigen),
             max_verebey_worsening = pk_verebey_worsening(o$mean[grepl("^Verebey", o$study)], pr[grepl("^Verebey", o$study)], pl[grepl("^Verebey", o$study)]),
             ka = f$par$ka, fp_4.5mg = fp_at_dose(4.5, f$par), parent_loss_4.5mg = ploss_at_dose(4.5, f$par),
             bioavailability_4.5mg = pk_bioavailability(4.5, f$par), bioavailability_100mg = pk_bioavailability(100, f$par),
             convergence = f$convergence, stringsAsFactors = FALSE)
}
cand_table <- function(lin, alts, o) {
  ct <- do.call(rbind, c(list(cand_row(lin, lin, o)), lapply(alts, cand_row, lin = lin, o = o)))
  rule <- pk_adoption_rule(ct); ct$qualifies <- rule$qualifies; ct$adopted <- ct$dose_model == rule$adopted
  list(table = ct, rule = rule)
}
primary <- cand_table(fit_lin, lapply(ALT, alt_fit), obs)
cand <- primary$table; rule <- primary$rule; adopted <- rule$adopted
message("Dose-model rule: ", rule$reason)
# Sensitivity of the decision to the ASSUMED day-21 CV (0.5 and 1.2; primary 0.8): same rule, refitted linear reference.
cv_sens <- do.call(rbind, lapply(BESTED_CV_SENS, function(cv) {
  lin_cv <- best_of(sprintf("linear|cv=%g", cv)); lin_cv$dose_model <- "linear"
  t <- cand_table(lin_cv, lapply(ALT, alt_fit, cv = cv), obs_cv(cv))$table
  cbind(analysis = sprintf("Bested day-21 CV %g (primary 0.8)", cv), t) }))
# D50 sensitivity of both alternatives (10, 30 = primary, 100 mg), each against the primary linear fit.
d50_tab <- do.call(rbind, lapply(ALT, function(m) do.call(rbind, lapply(FP_D50_ALL, function(d)
  cbind(analysis = sprintf("D50 %g mg%s", d, if (d == FP_D50_PRIMARY) " (primary)" else ""), cand_row(alt_fit(m, d), fit_lin, obs),
        qualifies = NA, adopted = NA)))))

# ---- 4. The adopted model goes through the same gate ----
extra <- PK_DOSE_MODELS[[adopted]]$par
if (adopted == "linear") {
  fit <- fit_lin; hs_initial <- hs_initial_lin; gate_notes <- g_lin$notes; gate_fixed <- g_lin$fixed
} else {
  fa <- alt_fit(adopted); hs_initial <- pk_hessian_summary(fa)
  g_ad <- pk_identifiability_gate(fa, obs)
  fit <- g_ad$fit; fit$structure <- fit_lin$structure; fit$dose_model <- adopted
  gate_notes <- c(g_lin$notes, g_ad$notes); gate_fixed <- c(g_lin$fixed, g_ad$fixed)
}

# ---- 5. Sensitivity PK fits (NOT used by the primary analysis; consumed by scripts/03 for occupancy sensitivity) ----
# (1) fp at 100 mg fixed at 0.90 and 0.99 instead of 0.95 (adopted structure; start = primary estimates rescaled by (1 - fp)).
# (2) Terminal-phase sensitivity: the Verebey 1976 sentence "declined slowly from 24 through 72 hr from 2.4 to 1.7 ng/ml,
#     with an apparent half-life of 96 hr" (which dose it follows is not stated) appended as NTX C24/C72 daily rows from the
#     SEPARATE file data/training/pk_summary_sensitivity_terminal.csv; structure re-selected by the same rule, with the
#     adopted dose model (its extra parameter re-estimated; review M5).
# (3) Only if a dose-dependent structure is adopted: D50 = 10 and 100 mg (refits of the adopted structure; keys D50_10,
#     D50_100). Under the linear model the D50 fits are PK outputs only (pk_alt_*.rds, pk_dose_model_comparison.csv).
# File names and the `key` field use PK_SENS_KEYS (R/pk.R); scripts 10 and 11 select by key, never by label (review C1).
# The same batch also fits the one-compartment counterpart of an adopted dose-dependent structure, for the structure files.
rescale <- function(p, fp) {         # fp at the reference dose and the apparent parent parameters, by the same factor
  k <- (1 - fp) / (1 - p$fp)
  for (nm in c("CL", "V1", "Q", "V2")) p[[nm]] <- p[[nm]] * k
  p$fp <- fp; p
}
environment(rescale) <- globalenv()
term_raw <- read_dataset("data/training/pk_summary_sensitivity_terminal.csv")
obs_term <- rbind(obs[, names(term_raw)], term_raw)
sens_jobs <- list(list(kind = "fp", key = PK_SENS_KEYS[["fp090"]], fp = 0.90), list(kind = "fp", key = PK_SENS_KEYS[["fp099"]], fp = 0.99),
                  list(kind = "terminal", key = PK_SENS_KEYS[["terminal"]]))
if (adopted != "linear")
  sens_jobs <- c(sens_jobs, list(list(kind = "D50", key = PK_SENS_KEYS[["D50_10"]], d50 = 10), list(kind = "D50", key = PK_SENS_KEYS[["D50_100"]], d50 = 100),
                                 list(kind = "one", key = "one_compartment")))
# Every input is passed as an argument (never captured from the script environment) so that the task serialises to a worker.
run_sens <- function(j, obs, obs_term, par, pars, structure, fp_fixed, free, rescale) {
  if (identical(j$kind, "fp")) c(fit_pk(obs, rescale(par, j$fp), pars = pars), list(structure = structure))
  else if (identical(j$kind, "D50")) c(fit_pk(obs, modifyList(par, list(fp_D50 = j$d50)), pars = pars), list(structure = structure))
  else if (identical(j$kind, "one")) c(fit_pk(obs, modifyList(par, list(Q = 1e-6)), pars = setdiff(pars, c("Q", "V2"))), list(structure = "one-compartment"))
  else pk_select_structure(obs_term, modifyList(par, list(fp = fp_fixed)), free = free)
}
environment(run_sens) <- globalenv()
raw <- with_cluster(length(sens_jobs), function(cl)
  parallel::parLapply(cl, sens_jobs, run_sens, obs = obs, obs_term = obs_term, par = fit$par, pars = fit$pars,
                      structure = fit$structure, fp_fixed = FP_FIXED, free = setdiff(c(FREE_PARS, extra), gate_fixed), rescale = rescale))
names(raw) <- vapply(sens_jobs, `[[`, "", "key")
# compartments under the adopted dose model (structure files; same rule as step 1): linear -> the step-1 comparison;
# dose-dependent -> its one-compartment refit against the adopted (two- or one-compartment) fit.
if (adopted == "linear") sel <- sel_lin else if (fit$structure == "two-compartment" && !("Q" %in% gate_fixed)) {
  one_ad <- raw$one_compartment; two_ad <- fit; two_ad$start_used <- "adopted fit"
  sel <- list(one = one_ad, two = two_ad, lr_threshold = stats::qchisq(0.95, 2))
  sel$fit <- if (one_ad$value - two_ad$value > sel$lr_threshold) two_ad else one_ad
} else sel <- list(one = fit, two = sel_lin$two, lr_threshold = sel_lin$lr_threshold, fit = fit)   # one-compartment adopted: kept as selected in step 1
hs <- pk_hessian_summary(fit)
if (!hs$pd) warning("PK Hessian still not positive definite after the gate: SEs unavailable")
structure_note <- if (!identical(sel$fit$structure, fit$structure))
  sprintf("NEEDS A RULING: under the adopted dose model the compartment rule prefers the %s structure, but the dose models were compared on the %s structure selected under the linear model; the adopted fit keeps the %s structure.",
          sel$fit$structure, fit$structure, fit$structure) else
  sprintf("The compartment rule selects the %s structure under the adopted dose model as well.", fit$structure)
fit$fixed <- unique(c("fp", if (adopted != "linear") "fp_D50", gate_fixed))
fit$fp_note <- fp_note
d50_note <- if (adopted != "linear") sprintf("fp_D50 fixed at %g mg (assumed; anchored at %g mg, so it changes only the shape below %g mg; D50 = 10 and 100 mg refitted as sensitivities: pk_fit_sens_D50_10.rds, pk_fit_sens_D50_100.rds).", FP_D50_PRIMARY, FP_DREF_MG, FP_DREF_MG)
fit$fixed_notes <- c(fp_note, d50_note, gate_notes)
fit$hessian_summary <- hs; fit$hessian_summary_initial <- hs_initial
fit$fp_fixed <- FP_FIXED; fit$parent_half_lives <- pk_parent_half_lives(fit$par)
fit$fu_BN_note <- fu_BN_note; fit$fu_BN_fallback <- !is.finite(fu_BN)
fit$dose_model_comparison <- cand; fit$adoption_reason <- rule$reason
fit$bioavailability <- c(`4.5` = pk_bioavailability(4.5, fit$par), `100` = pk_bioavailability(100, fit$par))

# ---- 6a. Structure files (one vs two compartments, under the adopted dose model) ----
cmp_one <- pk_hessian_summary(sel$one); cmp_two <- pk_hessian_summary(sel$two)
n_distinct <- nrow(unique(data.frame(obs$regimen, obs$analyte, obs$stat)))
pred_fit <- pk_predict_stats(fit$par, obs)
structure_cmp <- data.frame(obs[, c("study", "dose_mg", "regimen", "analyte", "stat", "unit", "mean", "sd")],
  pred_one_compartment = pk_predict_stats(sel$one$par, obs), pred_two_compartment = pk_predict_stats(sel$two$par, obs),
  pred_selected = pred_fit, dose_model = adopted)
structure_cmp$ratio_one <- structure_cmp$pred_one_compartment / structure_cmp$mean
structure_cmp$ratio_two <- structure_cmp$pred_two_compartment / structure_cmp$mean
structure_sum <- data.frame(structure = c("one-compartment", "two-compartment"), dose_model = adopted, wss = c(sel$one$value, sel$two$value),
  n_parameters = c(length(sel$one$pars), length(sel$two$pars)), residual_df = n_distinct - c(length(sel$one$pars), length(sel$two$pars)),
  hessian_pd = c(cmp_one$pd, cmp_two$pd), min_eigenvalue = c(min(cmp_one$eigen), min(cmp_two$eigen)),
  convergence = c(sel$one$convergence, sel$two$convergence), lr_threshold = sel$lr_threshold, fp_fixed = FP_FIXED,
  selected = c(fit$structure == "one-compartment", fit$structure == "two-compartment"),
  two_compartment_start = c(NA, sel$two$start_used))
c24 <- structure_cmp$stat == "C24"
c24_ok <- function(r) all(r[c24] <= 1.5 & r[c24] >= 1 / 1.5)
c24_note <- sprintf("C24 (24 h after the acute dose): observed %s; one-compartment predicts %s; two-compartment predicts %s. %s",
  paste(sprintf("%s %.1f", structure_cmp$analyte[c24], structure_cmp$mean[c24]), collapse = ", "),
  paste(sprintf("%s %.2f", structure_cmp$analyte[c24], structure_cmp$pred_one_compartment[c24]), collapse = ", "),
  paste(sprintf("%s %.2f", structure_cmp$analyte[c24], structure_cmp$pred_two_compartment[c24]), collapse = ", "),
  if (!c24_ok(structure_cmp$ratio_one) && !c24_ok(structure_cmp$ratio_two))
    "NEITHER structure reproduces both 24 h values within a factor of 1.5: residual misfit reported as a limitation (no further structure added)."
  else "At least one structure reproduces both 24 h values within a factor of 1.5.")
write.csv(structure_cmp, file.path(RES, "pk_structure_comparison.csv"), row.names = FALSE)
write.csv(structure_sum, file.path(RES, "pk_structure_summary.csv"), row.names = FALSE)

# ---- 6b. Every PK output from the adopted fit (R/pk.R pk_write_outputs), with the fit notes ----
rw <- pk_row_wss(pred_fit, obs)
n_eff_pred <- length(unique(signif(pred_fit, 6))); n_eff_pred3 <- length(unique(signif(pred_fit, 3)))
tm <- data.frame(analyte = obs$analyte, mean = obs$mean, pred = pred_fit)[obs$stat == "Tmax", ]
chronic_bn <- obs$regimen == "daily" & obs$analyte == "BN" & obs$stat == "Cmax" & is_ver
clin <- function(x) x >= 0.05 & x <= 0.40      # labelled oral bioavailability 5-40% (review M1)
notes <- c(
  sprintf("Dose model: %s. Rule (pre-stated; R/pk.R pk_adoption_rule): %s", PK_DOSE_MODELS[[adopted]]$label, rule$reason),
  sprintf("Implied systemic bioavailability of the parent (1 - fp(D)) x (1 - parent loss(D)): %.4f at 4.5 mg and %.4f at 100 mg (labelled range 0.05-0.40: %s at 4.5 mg, %s at 100 mg). Comparison of all three structures, D50 10/30/100 mg and the Bested day-21 CV 0.5/1.2 sensitivity: pk_dose_model_comparison.csv; row-level predictions under each: pk_dose_model_vs_obs.csv.",
          fit$bioavailability[["4.5"]], fit$bioavailability[["100"]], if (clin(fit$bioavailability[["4.5"]])) "inside" else "OUTSIDE",
          if (clin(fit$bioavailability[["100"]])) "inside" else "OUTSIDE"),
  sprintf("Structure: %s (WSS one-compartment %.3f [Hessian PD %s], two-compartment %.3f [Hessian PD %s] under the %s dose model; two-compartment adopted only if WSS drops by > %.2f); see pk_structure_summary.csv and pk_structure_comparison.csv. %s",
          fit$structure, sel$one$value, cmp_one$pd, sel$two$value, cmp_two$pd, adopted, qchisq(0.95, 2), structure_note),
  c24_note,
  sprintf("Objective statistics: %s. Weights obj_weight / CV^2 with CV = between-subject SD / mean where reported (Verebey 1976; Bested 2023 day-1 SDs digitised from Figure 4), else the documented assumed CV (Bested 2023 day-21 medians: %s), else the floor %.2f (Tmax %.2f); the Bested day-1 profile rows carry obj_weight 1/4 (4 repeated measures per analyte). WSS: total %.3f, Verebey rows %.3f, Bested rows %.3f.",
          paste(fit$objective_stats, collapse = ", "), paste(unique(obs$cv_assumed[is.finite(obs$cv_assumed)]), collapse = "/"),
          PK_CV_FLOOR, PK_TMAX_CV_FLOOR, sum(rw), sum(rw[is_ver]), sum(rw[is_bes])),
  sprintf("Final WSS %.4f; Hessian eigenvalues (of hessian/2): %s; positive definite: %s",
          fit$value, paste(signif(hs$eigen, 4), collapse = ", "), hs$pd),
  sprintf("ka %.3f per h (absorption half-life %.1f min)", fit$par$ka, 60 * log(2) / fit$par$ka),
  # Effective residual df counts DISTINCT PREDICTED VALUES, not distinct (regimen, analyte, stat) labels: a terminal
  # half-life is the same quantity after a single dose and at steady state in a linear model (R1 fix round 4).
  sprintf("Observations %d (%d single-dose, %d once-daily at steady state; %d at 100 mg, %d at 4.5 mg) and %d nominally distinct (dose, regimen, analyte, stat) combinations, %d distinct predicted values to 6 significant figures and %d to 3. With %d fitted parameters: nominal residual degrees of freedom %d, EFFECTIVE %d (%d at 3 significant figures). The 8 Bested day-1 rows are repeated measures and carry weight 1/4 each.",
          nrow(obs), sum(obs$regimen == "single"), sum(obs$regimen == "daily"), sum(obs$dose_mg == 100), sum(obs$dose_mg == 4.5),
          nrow(unique(data.frame(obs$dose_mg, obs$regimen, obs$analyte, obs$stat))), n_eff_pred, n_eff_pred3, length(fit$pars),
          nrow(unique(data.frame(obs$dose_mg, obs$regimen, obs$analyte, obs$stat))) - length(fit$pars),
          n_eff_pred - length(fit$pars), n_eff_pred3 - length(fit$pars)),
  sprintf("Tmax check: %s", paste(sprintf("%s observed %.2f h, predicted %.2f h", tm$analyte, tm$mean, tm$pred), collapse = "; ")),
  sprintf("Parent naltrexone disposition half-lives: distribution %.2f h, terminal %.2f h; smallest Hessian eigenvalue %.3g (largest %.3g)",
          fit$parent_half_lives[["distribution"]], fit$parent_half_lives[["terminal"]], min(hs$eigen), max(hs$eigen)),
  paste("Fixed parameters:", fit$fixed_notes),
  sprintf("Limitation (stated; no structure added): the chronic (once-daily, steady-state) 6-beta-naltrexol Cmax is predicted %.1f vs observed %.1f ng/mL (%.2f-fold). With the 24 h concentrations in the fit, linear kinetics cannot also produce the observed 1.8-fold rise of the chronic over the acute peak; Verebey 1976 itself reports no trough accumulation. The observation's CV (0.57) keeps its weight low.",
          pred_fit[chronic_bn], obs$mean[chronic_bn], pred_fit[chronic_bn] / obs$mean[chronic_bn]),
  # the chronic-BN Cmax row is already spelled out above, so it is excluded here
  pk_misfit_notes(obs, pred_fit, exclude = chronic_bn),
  paste("Unbound fraction of 6-beta-naltrexol:", fu_BN_note))
pk_write_outputs(fit, obs, RES, notes)
cat(readLines(file.path(RES, "pk_fit_notes.txt")), sep = "\n")

# Dose-model comparison outputs (PK outputs only; not read by later scripts)
dm_tab <- rbind(cbind(analysis = "primary (D50 30 mg, Bested day-21 CV 0.8)", cand), d50_tab[d50_tab$fp_D50_mg != FP_D50_PRIMARY, names(cv_sens)],
                cv_sens)
write.csv(dm_tab, file.path(RES, "pk_dose_model_comparison.csv"), row.names = FALSE)
fits_cmp <- c(list(linear = fit_lin), stats::setNames(lapply(ALT, alt_fit), ALT))
dvo <- data.frame(obs[, c("study", "dose_mg", "regimen", "analyte", "stat", "unit", "mean", "sd", "n")], cv_used = pk_obs_cv(obs),
                  obj_weight = pk_obs_weight(obs))
for (m in names(fits_cmp)) { p <- pk_predict_stats(fits_cmp[[m]]$par, obs)
  dvo[[paste0("pred_", m)]] <- p; dvo[[paste0("fold_", m)]] <- p / obs$mean; dvo[[paste0("wss_", m)]] <- pk_row_wss(p, obs) }
write.csv(dvo, file.path(RES, "pk_dose_model_vs_obs.csv"), row.names = FALSE)
for (m in ALT) for (d in FP_D50_ALL) saveRDS(alt_fit(m, d), file.path(RES, sprintf("pk_alt_%s_D50_%g.rds", m, d)))

# ---- 7. Sensitivity records and files ----
sens_record <- function(f, key, label, fp, obs_used) {
  f$key <- key; f$sensitivity <- label; f$fp_fixed <- fp; f$hessian_summary <- pk_hessian_summary(f)
  f$dose_model <- adopted; f$parent_half_lives <- pk_parent_half_lives(f$par); f$fu_BN_note <- fu_BN_note
  f$pred_table <- data.frame(obs_used[, c("study", "regimen", "analyte", "stat", "mean")], pred = pk_predict_stats(f$par, obs_used))
  f
}
sens <- list()
for (j in sens_jobs[vapply(sens_jobs, `[[`, "", "kind") != "one"]) {
  r <- raw[[j$key]]
  sens[[j$key]] <- switch(j$kind,
    fp = sens_record(r, j$key, sprintf("fp at %g mg fixed at %.2f (primary %.2f)", FP_DREF_MG, j$fp, FP_FIXED), j$fp, obs),
    D50 = sens_record(r, j$key, sprintf("D50 of the %s structure fixed at %g mg (primary %g mg; refitted)", adopted, j$d50, FP_D50_PRIMARY), FP_FIXED, obs),
    terminal = { s <- sens_record(r$fit, j$key, "terminal-phase sensitivity (Verebey 2.4 -> 1.7 ng/mL, 24-72 h, appended)", FP_FIXED, obs_term)
                 s$structure_wss <- c(one = r$one$value, two = r$two$value); s })
}
# Not retained (R1 fix round 3): a variant without the naltrexone thalf rows was fitted and is degenerate (apparent parent
# CL -> 0.014 L/h, terminal half-life unbounded, smallest Hessian eigenvalue 7e-4). See R1-report.md.
old_sens <- list.files(RES, pattern = "^pk_fit_sens_.*[.]rds$", full.names = TRUE)
if (length(old_sens)) file.remove(old_sens)        # never leave a stale sensitivity fit for script 03's glob
for (nm in names(sens)) saveRDS(sens[[nm]], file.path(RES, sprintf("pk_fit_sens_%s.rds", nm)))
sens_row <- function(key, label, s, fp) {
  data.frame(sensitivity = key, label = label, structure = s$structure, dose_model = adopted, fp = fp,
             fp0 = if (!is.null(s$par$fp0)) s$par$fp0 else NA_real_, floss0 = if (!is.null(s$par$floss0)) s$par$floss0 else NA_real_,
             fp_D50 = if (adopted != "linear") s$par$fp_D50 else Inf, wss = s$value,
             hessian_pd = s$hessian_summary$pd, min_eigenvalue = min(s$hessian_summary$eigen), ka = s$par$ka,
             # apparent = divided by the bioavailability AT 100 mg (review M4)
             CL_apparent = s$par$CL / pk_bioavailability(FP_DREF_MG, s$par), V1_apparent = s$par$V1 / pk_bioavailability(FP_DREF_MG, s$par),
             bioavailability_4.5mg = pk_bioavailability(4.5, s$par), bioavailability_100mg = pk_bioavailability(FP_DREF_MG, s$par),
             parent_t_half_distribution_h = s$parent_half_lives[["distribution"]], parent_t_half_terminal_h = s$parent_half_lives[["terminal"]],
             CLm = s$par$CLm, Vm = s$par$Vm)
}
sens_tab <- rbind(sens_row("primary", sprintf("primary (fp %.2f at %g mg; %s)", FP_FIXED, FP_DREF_MG, adopted),
                           fit, FP_FIXED),
                  do.call(rbind, lapply(names(sens), function(nm) sens_row(nm, sens[[nm]]$sensitivity, sens[[nm]], sens[[nm]]$fp_fixed))))
write.csv(sens_tab, file.path(RES, "pk_sensitivity_summary.csv"), row.names = FALSE)
write.csv(cbind(sensitivity = "terminal", sens$terminal$pred_table), file.path(RES, "pk_sensitivity_terminal_vs_obs.csv"), row.names = FALSE)
print(sens_tab)

# ---- 8. Scripted check: every PK output describes the adopted fit (review C3) ----
chk <- pk_output_check(RES, fit, obs)
if (length(chk)) stop("PK outputs are inconsistent with the adopted fit: ", paste(chk, collapse = "; "))
runtime <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
writeLines(c(sprintf("Adopted dose model: %s (%s).", adopted, PK_DOSE_MODELS[[adopted]]$label), rule$reason,
             sprintf("D50 sensitivity: %s.", paste(sprintf("%s D50 %g mg: %s %.4g, dWSS %.2f, Hessian PD %s%s", d50_tab$dose_model, d50_tab$fp_D50_mg,
                                                             d50_tab$extra_parameter, d50_tab$extra_estimate, d50_tab$delta_wss, d50_tab$hessian_pd,
                                                             ifelse(d50_tab$at_boundary, " (BOUNDARY)", "")), collapse = "; ")),
             sprintf("Bested day-21 CV sensitivity: %s.", paste(sprintf("%s: %s dWSS %.2f%s", cv_sens$analysis, cv_sens$dose_model, cv_sens$delta_wss,
                                                                       ifelse(cv_sens$adopted, " [adopted]", "")), collapse = "; ")),
             if (adopted == "linear") "The linear model is adopted: the D50 fits are PK outputs only (pk_alt_*.rds, pk_dose_model_comparison.csv) and do not enter the occupancy sensitivity." else
               "A dose-dependent structure is adopted: its D50 = 10 and 100 mg refits enter the occupancy sensitivity (pk_fit_sens_D50_10.rds, pk_fit_sens_D50_100.rds).",
             sprintf("PK output check passed (pk_output_check). Script 01 runtime %.1f min.", runtime)),
           file.path(RES, "pk_dose_model_note.txt"))
cat(readLines(file.path(RES, "pk_dose_model_note.txt")), sep = "\n")
