# Tables (R3b) from result files only: no held-out test-data file is read here (held-out results only via output/results).
#   table1_module_screen.csv
#   table2_parameters.csv (+ table2_footnotes.csv)   PK, binding (adopted stage), posterior parameters labelled by computation
#   table3_heldout.csv (+ test_elpd_contrasts.csv)   PET, held-out trials, elpd with MCSE, contrasts with both SEs,
#                                                    per-trial lpd, responder scores, cytokines, QST
#   table4_virtual.csv                               virtual-trial summary, responder medians vs observed, ppp, power
#   tableS_virtual_omega2.csv, tableS_occupancy_sensitivity.csv, tableS_identifiability.csv, tableS_sd_sensitivity.csv,
#   tableS_responder_calibrations.csv, tableS_pk_structure.csv, tableS_bruun_check.csv, tableS_sobol.csv
# Every input is required; a missing one stops the script naming it. The adopted stage is read from binding_stage.txt and
# never replaced by another stage.
options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
res <- "output/results"
dir.create(res, recursive = TRUE, showWarnings = FALSE)

need <- function(...) {
  f <- c(...); miss <- f[!file.exists(f)]
  if (length(miss)) stop("Table script cannot run; missing input file(s): ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(f)
}
rd <- function(f) { need(f); utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }
rr <- function(name) rd(file.path(res, name))
need_cols <- function(d, cols, f) {
  m <- setdiff(cols, names(d))
  if (length(m)) stop(f, " lacks required column(s): ", paste(m, collapse = ", "), call. = FALSE)
}
pick_col <- function(d, pattern, what, f) {
  hit <- grep(pattern, names(d), value = TRUE)
  if (!length(hit)) stop(f, ": no ", what, " column found among: ", paste(names(d), collapse = ", "), call. = FALSE)
  hit[1]
}
wr <- function(d, name) write.csv(d, file.path(res, name), row.names = FALSE)

need(file.path(res, "binding_stage.txt"))
stage <- readLines(file.path(res, "binding_stage.txt"))[1]
need("data/module_screen.csv", "data/params/binding_constants.csv", "data/priors/placebo_priors.csv",
     file.path(res, c("stage_decision.csv", "pk_fit.rds", "pk_fit_table.csv", "pk_hessian.csv", "pk_fit_notes.txt",
                      sprintf("binding_params_%s.rds", c("A", "A2", stage)), "koff_provenance.csv", "training_fit.csv",
                      "pet_metrics.csv", "occupancy_ed50_by_stage.csv", "test_predictions.csv", "test_elpd.csv",
                      "test_elpd_contrast_mcse.csv", "responder_scores.csv", "cytokine_test.csv", "qst_test.csv",
                      "virtual_summary.csv", "power_curve.csv", "effect_eligibility.csv",
                      "steady_state_occupancy_by_stage.csv", "p_loo_identifiability.csv", "omegas.csv",
                      "omegas_sd_sensitivity.csv", "virtual_summary_sd_sensitivity.csv", "power_curve_sd_sensitivity.csv",
                      "pk_structure_summary.csv", "pk_structure_comparison.csv", "pk_sensitivity_summary.csv",
                      "bruun_check.csv", "sobol.csv", "sobol_seed_agreement.csv")))
dec <- rr("stage_decision.csv")
if (!identical(dec$stage_adopted[1], stage))
  stop("binding_stage.txt says ", stage, " but stage_decision.csv says ", dec$stage_adopted[1], call. = FALSE)

# ---- Table 1 ----
wr(rd("data/module_screen.csv"), "table1_module_screen.csv")

# ---- Table 2: parameters ----
T2COLS <- c("section", "parameter", "estimate", "unit", "se", "se_scale", "lo95", "hi95", "source_type", "source", "note")
t2row <- function(section, parameter, estimate, unit = "", se = NA_real_, se_scale = "", lo95 = NA_real_, hi95 = NA_real_,
                  source_type, source, note = "")
  data.frame(section = section, parameter = parameter, estimate = estimate, unit = unit, se = se, se_scale = se_scale,
             lo95 = lo95, hi95 = hi95, source_type = source_type, source = source, note = note, stringsAsFactors = FALSE)
foot <- list()
addfoot <- function(id, text) foot[[length(foot) + 1]] <<- data.frame(footnote = id, text = text, stringsAsFactors = FALSE)

# PK: the adopted structure, SEs from the Hessian (gate: positive definite), fp fixed with its rationale
pk <- readRDS(file.path(res, "pk_fit.rds"))
d50_bd <- if (pk_dose_model_of(pk$par) != "linear") d50_boundary(res, pk_dose_model_of(pk$par)) else numeric(0)
pkt <- rr("pk_fit_table.csv"); need_cols(pkt, c("parameter", "estimate", "se_transformed", "lo95", "hi95", "transform", "source_type"), "pk_fit_table.csv")
hes <- rr("pk_hessian.csv"); need_cols(hes, c("fit", "eigenvalue", "positive_definite"), "pk_hessian.csv")
hf <- hes[hes$fit == "final", ]
if (!nrow(hf)) stop("pk_hessian.csv has no 'final' rows", call. = FALSE)
pd <- isTRUE(as.logical(hf$positive_definite[1]))
pk_unit <- c(ka = "1/h", CL = "L/h", V1 = "L", Q = "L/h", V2 = "L", CLm = "L/h", Vm = "L", fp = "fraction (at 100 mg)",
             fp0 = "fraction (dose -> 0)", floss0 = "fraction (dose -> 0)", fp_D50 = "mg")
fixed_pk <- grepl("^fixed", pkt$source_type)
fp_note <- if (!is.null(pk$fp_note)) pk$fp_note else pk$fixed_notes[grepl("^fp fixed", pk$fixed_notes)]
if (any(pkt$parameter == "fp") && !length(fp_note)) stop("pk_fit.rds carries no rationale for the fixed fp", call. = FALSE)
pk_rows <- t2row("PK", pkt$parameter, pkt$estimate, ifelse(pkt$parameter %in% names(pk_unit), pk_unit[pkt$parameter], ""),
  se = if (pd) pkt$se_transformed else NA_real_, se_scale = ifelse(fixed_pk, "", paste(pkt$transform, "scale")),
  lo95 = if (pd) pkt$lo95 else NA_real_, hi95 = if (pd) pkt$hi95 else NA_real_,
  source_type = ifelse(fixed_pk, sprintf("fixed at %g", pkt$estimate), sprintf("fitted (training PK summaries; %s; dose model %s)", pk$structure,
                                                                            pk_dose_model_of(pk$par))),
  source = ifelse(fixed_pk, "pk_fit_notes.txt", "This study: fitted to data/training/pk_summary.csv"),
  note = ifelse(pkt$parameter == "fp", if (length(fp_note)) fp_note[1] else "",
                if (pd) sprintf("SE from the positive-definite Hessian (smallest eigenvalue %.3g)", min(hf$eigenvalue))
                else "Hessian NOT positive definite: no SE reported"))
pk_rows$note[fixed_pk & pkt$parameter != "fp"] <- vapply(pkt$parameter[fixed_pk & pkt$parameter != "fp"], function(p) {
  n <- pk$fixed_notes[grepl(p, pk$fixed_notes, fixed = TRUE)]; if (length(n)) n[1] else "see pk_fit_notes.txt" }, "")
addfoot("PK fit notes", paste(readLines(file.path(res, "pk_fit_notes.txt"), warn = FALSE), collapse = " | "))

# Binding: the adopted stage's parameters
bc <- read_dataset("data/params/binding_constants.csv")
src_of <- function(cmp, target, quantity) {
  r <- bc$analyte == cmp & bc$target == target & bc$quantity == quantity
  if (!any(r)) stop("binding_constants.csv has no row for ", cmp, " ", target, " ", quantity, call. = FALSE)
  hum <- r & bc$species == "human"; use <- if (any(hum)) hum else r
  txt <- paste(unique(paste0(bc$source[use], " [organism: ", bc$species[use], "]")), collapse = "; ")
  if (sum(use) > 1) txt <- paste0("geometric mean of ", sum(use), " rows: ", txt)
  if (!any(hum)) txt <- paste0("non-human fallback (no human row): ", txt)
  txt
}
bp <- readRDS(file.path(res, sprintf("binding_params_%s.rds", stage)))
bpA <- readRDS(file.path(res, "binding_params_A.rds")); bpA2 <- readRDS(file.path(res, "binding_params_A2.rds"))
kp <- rr("koff_provenance.csv"); need_cols(kp, c("analyte", "koff_per_h", "koff_assumed", "source"), "koff_provenance.csv")
kpa <- function(an) { r <- kp[kp$analyte == an, ]; if (nrow(r) != 1) stop("koff_provenance.csv needs one ", an, " row", call. = FALSE); r }
base_kd <- c(NTX = bpA$Kd_NTX, BN = bpA$Kd_BN)
sec_b <- sprintf("Binding (stage %s, adopted)", stage)
if (stage == "J") {
  for (x in c("log_scales", "log_scales_se", "kd_scale", "koff_scale")) if (is.null(bp[[x]])) stop("binding_params_J.rds lacks ", x, call. = FALSE)
  keq_tab <- rr("binding_J_keq_sensitivity.csv"); need_cols(keq_tab, c("keq", "wsse", "mean_ss_occupancy_4.5mg", "adopted"), "binding_J_keq_sensitivity.csv")
  unc <- rr("occupancy_uncertainty_J.csv")
  ci <- function(i) exp(bp$log_scales[[i]] + c(-1, 1) * 1.96 * bp$log_scales_se[[i]])
  kci <- ci(1); fci <- ci(2)
  calib <- "This study: fitted to all 11 PET rows (Lee 1988, Rabiner 2011, Weerts 2008), precision-weighted; post hoc, exploratory"
  b_rows <- rbind(
    t2row(sec_b, "kd_scale (shared by NTX and 6-beta-naltrexol)", bp$kd_scale, "multiplier", bp$log_scales_se[[1]], "log scale", kci[1], kci[2],
          "calibrated on PET (stage J)", calib, sprintf("correlation of the two log scales %.2f", stats::cov2cor(bp$log_scales_cov)[1, 2])),
    t2row(sec_b, "koff_scale (shared by NTX and 6-beta-naltrexol)", bp$koff_scale, "multiplier", bp$log_scales_se[[2]], "log scale", fci[1], fci[2],
          "calibrated on PET (stage J)", calib, ""),
    t2row(sec_b, "Kd_NTX", bp$Kd_NTX, "nM", bp$log_scales_se[[1]], "log scale", base_kd[["NTX"]] * kci[1], base_kd[["NTX"]] * kci[2],
          "derived: in vitro Ki x kd_scale", src_of("NTX", "MOR", "Ki"), sprintf("in vitro base %.3g nM", base_kd[["NTX"]])),
    t2row(sec_b, "Kd_BN", bp$Kd_BN, "nM", bp$log_scales_se[[1]], "log scale", base_kd[["BN"]] * kci[1], base_kd[["BN"]] * kci[2],
          "derived: in vitro Ki x kd_scale", src_of("BN", "MOR", "Ki"), sprintf("in vitro base %.3g nM", base_kd[["BN"]])),
    t2row(sec_b, "koff_NTX", bp$koff_NTX, "1/h", bp$log_scales_se[[2]], "log scale", kpa("NTX")$koff_per_h * fci[1], kpa("NTX")$koff_per_h * fci[2],
          "derived: base koff x koff_scale", kpa("NTX")$source,
          sprintf("base %g per h (%s); residence half-life %.1f h", kpa("NTX")$koff_per_h, if (kpa("NTX")$koff_assumed) "ASSUMED" else "in vitro", log(2) / bp$koff_NTX)),
    t2row(sec_b, "koff_BN", bp$koff_BN, "1/h", bp$log_scales_se[[2]], "log scale", kpa("BN")$koff_per_h * fci[1], kpa("BN")$koff_per_h * fci[2],
          "derived: base koff x koff_scale", kpa("BN")$source,
          sprintf("base value %g per h is %s; the 10:1 base koff ratio is a structural sensitivity (occupancy_sensitivity_J.csv)", kpa("BN")$koff_per_h,
                  if (kpa("BN")$koff_assumed) "ASSUMED (no curated 6-beta-naltrexol koff)" else "in vitro")),
    t2row(sec_b, "keq (brain equilibration)", bp$keq, "1/h", source_type = sprintf("fixed at %g per h", bp$keq),
          source = "assumed (rapid passive blood-brain barrier permeation)",
          note = sprintf("Sensitivity: refitted at keq %s per h (weighted SSE %s; mean steady-state occupancy at 4.5 mg %s); %s",
                         paste(keq_tab$keq, collapse = " / "), paste(sprintf("%.2f", keq_tab$wsse), collapse = " / "),
                         paste(sprintf("%.3f", keq_tab$`mean_ss_occupancy_4.5mg`), collapse = " / "), bp$keq_reason)))
} else {
  b_rows <- rbind(
    t2row(sec_b, "Kd_NTX", bp$Kd_NTX, "nM", source_type = if (is.null(bp$kd_scale)) "in vitro" else sprintf("in vitro Ki x kd_scale %.3g", bp$kd_scale), source = src_of("NTX", "MOR", "Ki")),
    t2row(sec_b, "Kd_BN", bp$Kd_BN, "nM", source_type = if (is.null(bp$kd_scale)) "in vitro" else sprintf("in vitro Ki x kd_scale %.3g", bp$kd_scale), source = src_of("BN", "MOR", "Ki")),
    t2row(sec_b, "koff_NTX", bp$koff_NTX, "1/h", source_type = if (kpa("NTX")$koff_assumed) "assumed" else "in vitro", source = kpa("NTX")$source),
    t2row(sec_b, "koff_BN", bp$koff_BN, "1/h", source_type = if (kpa("BN")$koff_assumed) "assumed" else "in vitro", source = kpa("BN")$source),
    t2row(sec_b, "keq (brain equilibration)", bp$keq, "1/h", source_type = sprintf("fixed at %g per h", bp$keq), source = "assumed (rapid passive blood-brain barrier permeation)"))
}
if (stage != "A2") b_rows <- rbind(b_rows, t2row("Binding (sensitivity stages, not adopted)", "kd_scale, stage A2", bpA2$kd_scale, "multiplier",
  source_type = "calibrated on the Rabiner rows only (stage A2)", source = "This study: Kd scale fitted to Rabiner 2011 occupancy; koff in vitro",
  note = "sensitivity bound; downstream results in output/results/sensitivity_stage_A2/"))

# Posterior parameters: source_type by computation (posterior / prior SD; >= 0.9 is prior-determined)
tf <- rr("training_fit.csv"); need_cols(tf, c("hypothesis", "variable", "sd"), "training_fit.csv")
qmed <- pick_col(tf, "^(X?50\\.?%?|q50|median)$", "posterior median", "training_fit.csv")
qlo  <- pick_col(tf, "^(X?2\\.5\\.?%?|q2\\.5|q025)$", "posterior 2.5% quantile", "training_fit.csv")
qhi  <- pick_col(tf, "^(X?97\\.5\\.?%?|q97\\.5|q975)$", "posterior 97.5% quantile", "training_fit.csv")
priors <- read_dataset("data/priors/placebo_priors.csv")
prior_sd <- prior_sd_on_sampler_scale(tf$variable, priors)
post <- t2row(paste("Posterior,", tf$hypothesis), paste(tf$hypothesis, tf$variable), tf[[qmed]], "sampler scale", tf$sd, "posterior SD",
              tf[[qlo]], tf[[qhi]], posterior_source_type(tf$sd, prior_sd),
              "Posterior median and 95% CrI; priors in data/priors/placebo_priors.csv and R/hypotheses.R",
              sprintf("prior SD %.4g on the sampler scale; posterior SD / prior SD = %.2f", prior_sd, tf$sd / prior_sd))
t2 <- rbind(pk_rows, b_rows, post)[, T2COLS]
nosrc <- is.na(t2$source) | trimws(t2$source) == ""
if (any(nosrc)) stop("Table 2 has parameters without a source: ", paste(t2$parameter[nosrc], collapse = ", "), call. = FALSE)
wr(t2, "table2_parameters.csv")
addfoot("prior-determined", sprintf("Posterior parameters are labelled prior-determined when posterior SD / prior SD >= %.1f (computed in scripts/11_tables.R from training_fit.csv and the priors; R/reporting.R), otherwise fitted. %d of %d posterior rows are prior-determined.",
                                    PRIOR_DETERMINED_RATIO, sum(grepl("^prior-determined", post$source_type)), nrow(post)))
addfoot("pbo_emax_logit", "pbo_emax_logit (Byon 2010) is curated in data/priors/placebo_priors.csv for reference only and is used by no code.")
addfoot("stage", sprintf("Binding stage %s adopted (stage_decision.csv): %s", stage, dec$reason[1]))
wr(do.call(rbind, foot), "table2_footnotes.csv")

# ---- Table 3: held-out results, long format ----
T3COLS <- c("block", "hypothesis", "item", "observed", "pred_median", "pred_lo90", "pred_hi90", "true_lo90", "true_hi90", "inside90",
            "metric", "value", "se_observation_level", "mcse_monte_carlo", "footnote")
long <- function(block, hypothesis, item, observed = NA, pred_median = NA, pred_lo90 = NA, pred_hi90 = NA, true_lo90 = NA,
                 true_hi90 = NA, inside90 = NA, metric = NA_character_, value = NA_real_, se_obs = NA_real_, mcse = NA_real_,
                 footnote = NA_character_)
  data.frame(block = block, hypothesis = hypothesis, item = item, observed = observed, pred_median = pred_median,
             pred_lo90 = pred_lo90, pred_hi90 = pred_hi90, true_lo90 = true_lo90, true_hi90 = true_hi90, inside90 = inside90,
             metric = metric, value = value, se_observation_level = se_obs, mcse_monte_carlo = mcse, footnote = footnote,
             stringsAsFactors = FALSE)

pm <- rr("pet_metrics.csv"); need_cols(pm, c("stage", "n", "RMSE", "MAE", "n_within_15", "rows"), "pet_metrics.csv")
pet_role <- function(s) ifelse(grepl("^A2", s), "A2 (Kd calibrated on the Rabiner rows; sensitivity)",
                        ifelse(grepl("^A", s), "A (bottom-up): HELD-OUT prediction, fits no PET row",
                        ifelse(grepl("^J", s), "J (adopted): CALIBRATION fit to all rows, not a validation", s)))
met <- c("n", "RMSE", "MAE", "n_within_15", if ("wSSE" %in% names(pm)) "wSSE")
b_pet <- do.call(rbind, lapply(seq_len(nrow(pm)), function(i)
  long("PET occupancy (Lee 1988, Rabiner 2011, Weerts 2008)", "all (hypothesis-independent)", paste0("stage ", pm$stage[i], ": ", pm$rows[i]),
       metric = met, value = as.numeric(pm[i, met]), footnote = pet_role(pm$stage[i]))))
ed <- rr("occupancy_ed50_by_stage.csv"); need_cols(ed, c("stage", "model_ed50_mg", "rabiner_ed50_mg", "rabiner_lo95", "rabiner_hi95", "note"), "occupancy_ed50_by_stage.csv")
b_ed <- long("PET: single-dose occupancy ED50 at 4 h vs Rabiner 2011", "all (hypothesis-independent)", paste("stage", ed$stage),
             observed = ed$rabiner_ed50_mg, pred_median = ed$model_ed50_mg, metric = "ED50_mg", value = ed$model_ed50_mg,
             footnote = sprintf("Rabiner %.2f mg (95%% CI %.2f-%.2f); %s", ed$rabiner_ed50_mg, ed$rabiner_lo95, ed$rabiner_hi95, ed$note))

tp <- rr("test_predictions.csv")
need_cols(tp, c("hypothesis", "trial", "endpoint", "day", "observed", "pred_median", "pred_lo90", "pred_hi90", "true_lo90", "true_hi90",
                "lpd", "lpd_mcse", "inside90", "role", "obs_se", "interval_definition"), "test_predictions.csv")
dd <- tp[tp$endpoint == "diff_nrs", , drop = FALSE]
if (!nrow(dd)) stop("test_predictions.csv has no diff_nrs rows", call. = FALSE)
dd <- dd[order(match(tolower(dd$trial), c("final", "innova")), dd$day, dd$hypothesis), ]
b_tr <- long("Held-out trials (FINAL, INNOVA): LDN minus placebo difference in NRS change", dd$hypothesis,
             sprintf("%s day %g (%s)", dd$trial, dd$day, dd$role), dd$observed, dd$pred_median, dd$pred_lo90, dd$pred_hi90,
             dd$true_lo90, dd$true_hi90, dd$inside90, metric = "lpd", value = dd$lpd, se_obs = dd$obs_se, mcse = dd$lpd_mcse,
             footnote = "pred_lo90/pred_hi90: posterior predictive interval for the OBSERVED difference (inside90 judged against it); true_lo90/true_hi90: interval of the TRUE difference; se_observation_level: the observation SE used for scoring")
el <- rr("test_elpd.csv")
need_cols(el, c("hypothesis", "elpd_primary", "mcse_elpd_primary", "primary_endpoints", "elpd_all_days_secondary_nonindependent",
                "mcse_elpd_all_days_secondary_nonindependent"), "test_elpd.csv")
b_el <- rbind(
  long(paste0("Held-out trials: PRIMARY elpd, one endpoint per trial (", el$primary_endpoints[1], ")"), el$hypothesis, "primary endpoints",
       metric = "elpd_primary", value = el$elpd_primary, mcse = el$mcse_elpd_primary,
       footnote = "mcse_monte_carlo: Monte Carlo SE of the elpd from the predictive draws; it does not cover refits or other trials"),
  long("Held-out trials: SECONDARY, non-independent sum over all scored days (INNOVA days correlated)", el$hypothesis, "all days",
       metric = "elpd_all_days_secondary_nonindependent", value = el$elpd_all_days_secondary_nonindependent,
       mcse = el$mcse_elpd_all_days_secondary_nonindependent))
pr <- tp[tp$endpoint == "diff_nrs" & tp$role == "primary", ]
pr <- pr[order(pr$hypothesis, match(tolower(pr$trial), c("final", "innova"))), ]
b_pt <- long("Held-out trials: per-trial log predictive density at the primary endpoint", pr$hypothesis, sprintf("%s day %g", pr$trial, pr$day),
             metric = "lpd", value = pr$lpd, se_obs = pr$obs_se, mcse = pr$lpd_mcse)
contr <- elpd_contrasts(tp, rr("test_elpd_contrast_mcse.csv"))
if (anyNA(contr$mcse_diff_monte_carlo)) stop("test_elpd_contrast_mcse.csv lacks a Monte Carlo SE for some contrast", call. = FALSE)
wr(contr, "test_elpd_contrasts.csv")
b_ct2 <- long("Held-out trials: pairwise PRIMARY elpd contrasts (a minus b)", contr$hypothesis_a, paste(contr$hypothesis_a, "vs", contr$hypothesis_b),
              metric = "elpd_diff_a_minus_b", value = contr$elpd_diff_a_minus_b, se_obs = contr$se_diff_observation_level,
              mcse = contr$mcse_diff_monte_carlo,
              footnote = sprintf("diff / observation-level SE %.2f; diff / Monte Carlo SE %.1f; per-trial lpd differences %s. %s%s",
                                 contr$ratio_diff_over_se, contr$ratio_diff_over_mcse, contr$per_trial_lpd_diff, contr$note,
                                 ifelse(contr$se_collapsed_flag, paste0(" FLAG: ", contr$flag_note), "")))
rs <- rr("responder_scores.csv")
need_cols(rs, c("hypothesis", "trial", "cut_pct", "obs_diff", "pred_diff_median", "pred_diff_lo90", "pred_diff_hi90", "ppp_two_sided",
                "pp_prob_counts_joint", "obs_rate_ldn", "obs_rate_pbo", "pred_rate_ldn_median", "pred_rate_pbo_median",
                "fisher_p_two_sided"), "responder_scores.csv")
rs <- rs[order(match(tolower(rs$trial), c("final", "innova")), rs$cut_pct, rs$hypothesis), ]
b_rs <- long("Held-out trials: responder scores (difference in responder proportions, LDN minus placebo)", rs$hypothesis,
             sprintf("%s, >= %d%% reduction", rs$trial, rs$cut_pct), rs$obs_diff, rs$pred_diff_median, rs$pred_diff_lo90, rs$pred_diff_hi90,
             metric = "ppp_two_sided", value = rs$ppp_two_sided,
             footnote = sprintf("observed rates LDN %.2f / placebo %.2f, predicted medians %.2f / %.2f; P(observed counts) %.3g; source-trial Fisher p %.3f (descriptive)",
                                rs$obs_rate_ldn, rs$obs_rate_pbo, rs$pred_rate_ldn_median, rs$pred_rate_pbo_median, rs$pp_prob_counts_joint, rs$fisher_p_two_sided))
ct <- rr("cytokine_test.csv")
need_cols(ct, c("cytokine", "pred_median", "pred_lo90", "pred_hi90", "hypothesis", "pct_change", "predicted"), "cytokine_test.csv")
ct$predicted <- as.logical(ct$predicted)
cp <- ct[ct$predicted %in% TRUE, , drop = FALSE]
if (!nrow(cp)) stop("cytokine_test.csv has no predicted cytokines", call. = FALSE)
unpred <- sort(setdiff(unique(ct$cytokine), unique(cp$cytokine)))
fn <- if (length(unpred)) paste("No prediction (no curated half-life):", paste(unpred, collapse = ", ")) else NA_character_
b_ct <- long("Cytokines (Parkitny 2017), % change", cp$hypothesis, cp$cytokine, cp$pct_change, cp$pred_median, cp$pred_lo90, cp$pred_hi90,
             inside90 = cp$pct_change >= cp$pred_lo90 & cp$pct_change <= cp$pred_hi90, footnote = fn)
qs <- rr("qst_test.csv"); need_cols(qs, c("measure", "diff", "lo95", "hi95", "unit", "point_favours", "ci_excludes_zero"), "qst_test.csv")
b_q <- long("QST (FINAL, source-trial differences; descriptive, no model prediction)", "none (descriptive)", sprintf("%s (%s)", qs$measure, qs$unit),
            observed = qs$diff, metric = "diff_LDN_minus_placebo", value = qs$diff,
            footnote = sprintf("95%% CI %.2f to %.2f; point favours %s; CI excludes 0: %s", qs$lo95, qs$hi95, qs$point_favours, qs$ci_excludes_zero))
t3 <- rbind(b_pet, b_ed, b_tr, b_el, b_pt, b_ct2, b_rs, b_ct, b_q)[, T3COLS]
wr(t3, "table3_heldout.csv")

# ---- Table 4: virtual trials (omega as calibrated) ----
vs <- rr("virtual_summary.csv")
need_cols(vs, c("hypothesis", "trial", "omega_mult", "obs_diff", "diff_median", "diff_lo90", "diff_hi90", "ppp_two_sided"), "virtual_summary.csv")
if (!any(vs$omega_mult == 1)) stop("virtual_summary.csv has no omega_mult == 1 rows", call. = FALSE)
v1 <- vs[vs$omega_mult == 1, ]
pw <- rr("power_curve.csv"); need_cols(pw, c("hypothesis", "role", "trial", "n_per_arm", "rej_rate_with_model_unc"), "power_curve.csv")
el4 <- rr("effect_eligibility.csv"); need_cols(el4, c("hypothesis", "net_effect_pp", "eligible", "power_role"), "effect_eligibility.csv")
POWER_N <- c(50, 100, 200)
t4 <- do.call(rbind, lapply(seq_len(nrow(v1)), function(i) {
  v <- v1[i, ]; out <- data.frame(hypothesis = v$hypothesis, trial = v$trial, calibration = "primary",
    obs_diff = v$obs_diff, sim_diff_median = v$diff_median, sim_diff_lo90 = v$diff_lo90, sim_diff_hi90 = v$diff_hi90,
    ppp_two_sided_diff = v$ppp_two_sided, stringsAsFactors = FALSE)
  for (cut in c(15, 30, 50)) {
    for (a in c("ldn", "pbo")) {
      out[[sprintf("resp%d_%s_pred_median", cut, a)]] <- v[[sprintf("resp%d_%s_median", cut, a)]]
      out[[sprintf("resp%d_%s_observed", cut, a)]] <- v[[sprintf("resp%d_%s_obs", cut, a)]]
    }
    r <- rs[rs$hypothesis == v$hypothesis & tolower(rs$trial) == tolower(v$trial) & rs$cut_pct == cut, ]
    out[[sprintf("resp%d_ppp_two_sided", cut)]] <- if (nrow(r) == 1) r$ppp_two_sided else NA_real_
  }
  e <- el4[el4$hypothesis == v$hypothesis, ]
  out$net_effect_day84_pp <- if (nrow(e)) e$net_effect_pp else NA_real_
  p <- pw[pw$hypothesis == v$hypothesis & tolower(pw$trial) == tolower(v$trial) & is.finite(pw$rej_rate_with_model_unc), ]
  out$power_role <- if (nrow(p)) p$role[1] else if (nrow(e)) e$power_role else "not run"
  for (n in POWER_N) out[[sprintf("rejection_rate_n%d", n)]] <- if (any(p$n_per_arm == n)) p$rej_rate_with_model_unc[p$n_per_arm == n][1] else NA_real_
  out$note <- paste("Responder medians are simulated rates at the actual design; resp*_ppp_two_sided from responder_scores.csv.",
                    "rejection_rate_n*: one-sided 5% test per arm size; assurance for eligible drug hypotheses, type I error for H0, blank where not run.",
                    "Responder calibration is a REPORTED MISFIT: every hypothesis over-predicts INNOVA's responder rates (see tableS_responder_calibrations.csv).")
  out }))
wr(t4, "table4_virtual.csv")
if (any(vs$omega_mult == 2)) wr(vs[vs$omega_mult == 2, ], "tableS_virtual_omega2.csv")

# ---- Supplementary: structural occupancy sensitivity at 1.5, 4.5 and 6 mg ----
DOSES <- c(1.5, 4.5, 6)
byst <- rr("steady_state_occupancy_by_stage.csv"); need_cols(byst, c("stage", "dose_mg", "mean_RO"), "steady_state_occupancy_by_stage.csv")
occ_of <- function(sg) vapply(DOSES, function(d) { v <- byst$mean_RO[byst$stage == sg & abs(byst$dose_mg - d) < 1e-9]
  if (length(v) != 1) stop("steady_state_occupancy_by_stage.csv: no ", d, " mg row for stage ", sg, call. = FALSE); v }, 0)
srow <- function(analysis, kind, point, lo = rep(NA_real_, 3), hi = rep(NA_real_, 3), wsse = NA_real_, note = "") {
  d <- data.frame(analysis = analysis, kind = kind, wsse_pet = wsse, note = note, stringsAsFactors = FALSE)
  for (j in seq_along(DOSES)) { d[[sprintf("occ_%gmg", DOSES[j])]] <- point[j]; d[[sprintf("lo_%gmg", DOSES[j])]] <- lo[j]; d[[sprintf("hi_%gmg", DOSES[j])]] <- hi[j] }
  d }
s_occ <- list(srow(sprintf("stage %s (adopted), point", stage), "point estimate", occ_of(stage)))
if (stage == "J") {
  unc <- rr("occupancy_uncertainty_J.csv"); osn <- rr("occupancy_sensitivity_J.csv")
  need_cols(osn, c("analysis", "wsse", "mean_ss_occupancy_1.5mg", "mean_ss_occupancy_4.5mg", "mean_ss_occupancy_6mg"), "occupancy_sensitivity_J.csv")
  ub <- function(src, col) vapply(DOSES, function(d) { v <- unc[[col]][unc$uncertainty_source == src & abs(unc$dose_mg - d) < 1e-9]
    if (length(v) != 1) stop("occupancy_uncertainty_J.csv: no ", src, " row at ", d, " mg", call. = FALSE); v }, 0)
  for (src in c("binding_only", "binding_and_pk"))
    s_occ[[length(s_occ) + 1]] <- srow(sprintf("Laplace 95%% band, %s", c(binding_only = "binding only", binding_and_pk = "binding + PK")[[src]]), "probability band", ub(src, "mean_RO_median"),
                                       ub(src, "mean_RO_lo95"), ub(src, "mean_RO_hi95"), note = unc$method[unc$uncertainty_source == src][1])
  s_occ[[length(s_occ) + 1]] <- srow("keq 0.5 to 8 per h (refitted; min and max)", "structural sensitivity", ub("keq_sensitivity", "mean_RO_point"),
                                     ub("keq_sensitivity", "mean_RO_lo95"), ub("keq_sensitivity", "mean_RO_hi95"))
  need_cols(osn, "key", "occupancy_sensitivity_J.csv")
  prim_occ <- osn$`mean_ss_occupancy_4.5mg`[osn$key == OCC_SENS_KEYS[["primary"]]][1]; prim_wsse <- osn$wsse[osn$key == OCC_SENS_KEYS[["primary"]]][1]
  for (i in seq_len(nrow(osn))) if (osn$key[i] != OCC_SENS_KEYS[["primary"]])
    s_occ[[length(s_occ) + 1]] <- srow(osn$analysis[i], "structural sensitivity (refit)",
      c(osn$`mean_ss_occupancy_1.5mg`[i], osn$`mean_ss_occupancy_4.5mg`[i], osn$`mean_ss_occupancy_6mg`[i]), wsse = osn$wsse[i],
      note = if ("fp_D50_mg" %in% names(osn) && osn$key[i] %in% PK_SENS_KEYS[c("D50_10", "D50_100")] && osn$fp_D50_mg[i] %in% d50_bd)
        "BOUNDARY: the dose-dependent parameter of this PK refit sits at its bound (near-singular Hessian); not an ordinary sensitivity"
        else if (osn$key[i] == OCC_SENS_KEYS[["koff_eq"]])   # computed here, so the table never carries a fixed ranking claim
        sprintf("ASSUMED 10:1 base koff ratio set to 1:1: a worse PET fit (weighted SSE %.1f against %.1f); occupancy at 4.5 mg %+.1f pp",
                osn$wsse[i], prim_wsse, 100 * (osn$`mean_ss_occupancy_4.5mg`[i] - prim_occ))
        else if ("note" %in% names(osn)) osn$note[i] else "")
}
for (sg in setdiff(c("A2", "A"), stage)) s_occ[[length(s_occ) + 1]] <- srow(sprintf("stage %s", sg), "outer bound (other binding stage)", occ_of(sg))
wr(do.call(rbind, s_occ), "tableS_occupancy_sensitivity.csv")

# ---- Supplementary: identifiability (p_loo, prior-determined fraction) ----
pl <- rr("p_loo_identifiability.csv")
need_cols(pl, c("hypothesis", "n_free_parameters", "p_loo", "prior_determined_fraction", "max_pareto_k", "p_loo_reliable", "note"), "p_loo_identifiability.csv")
pl$n_params_prior_determined_by_sd_ratio <- vapply(pl$hypothesis, function(h) sum(grepl("^prior-determined", post$source_type[tf$hypothesis == h])), 0L)
wr(pl, "tableS_identifiability.csv")

# ---- Supplementary: the SD-sensitivity calibration (consistency check, not a candidate calibration) ----
CHECK_LABEL <- "consistency check, not a candidate calibration"
osd <- rr("omegas_sd_sensitivity.csv")
need_cols(osd, c("hypothesis", "calibration", "omega_p", "omega_e", "achieved_sd_pbo", "target_sd_pbo", "achieved_sd_ldn", "target_sd_ldn", "at_bound"), "omegas_sd_sensitivity.csv")
osd$role <- ifelse(startsWith(osd$calibration, "SENSITIVITY"), CHECK_LABEL, "primary calibration")
wr(osd, "tableS_sd_sensitivity.csv")
vsd <- rr("virtual_summary_sd_sensitivity.csv"); need_cols(vsd, c("hypothesis", "trial"), "virtual_summary_sd_sensitivity.csv")
rc <- do.call(rbind, lapply(seq_len(nrow(v1)), function(i) { v <- v1[i, ]
  s <- vsd[vsd$hypothesis == v$hypothesis & vsd$trial == v$trial, ]
  if (nrow(s) != 1) stop("virtual_summary_sd_sensitivity.csv: no single row for ", v$hypothesis, " ", v$trial, call. = FALSE)
  do.call(rbind, lapply(c(15, 30, 50), function(cut) data.frame(hypothesis = v$hypothesis, trial = v$trial, cut_pct = cut,
    observed_ldn = v[[sprintf("resp%d_ldn_obs", cut)]], observed_pbo = v[[sprintf("resp%d_pbo_obs", cut)]],
    primary_pred_ldn = v[[sprintf("resp%d_ldn_median", cut)]], primary_pred_pbo = v[[sprintf("resp%d_pbo_median", cut)]],
    sd_sensitivity_pred_ldn = s[[sprintf("resp%d_ldn_median", cut)]], sd_sensitivity_pred_pbo = s[[sprintf("resp%d_pbo_median", cut)]],
    sd_sensitivity_label = CHECK_LABEL,
    note = "REPORTED MISFIT: the primary calibration over-predicts INNOVA's responder rates and the SD-sensitivity calibration under-predicts FINAL's; no calibration is chosen from test data")))
}))
wr(rc, "tableS_responder_calibrations.csv")

# ---- Supplementary: PK structure comparison (one vs two compartments, with C24), and the PK sensitivity fits ----
pss <- rr("pk_structure_summary.csv"); psc <- rr("pk_structure_comparison.csv"); pse <- rr("pk_sensitivity_summary.csv")
pdm <- rr("pk_dose_model_comparison.csv")
need_cols(pdm, c("analysis", "dose_model", "wss", "wss_verebey", "wss_bested", "delta_wss", "extra_parameter", "extra_estimate", "extra_se_logit",
                 "hessian_pd", "max_verebey_worsening", "bioavailability_4.5mg", "bioavailability_100mg", "at_boundary", "adopted"), "pk_dose_model_comparison.csv")
need_cols(pss, c("structure", "wss", "n_parameters", "selected"), "pk_structure_summary.csv")
need_cols(psc, c("study", "regimen", "analyte", "stat", "mean", "pred_one_compartment", "pred_two_compartment"), "pk_structure_comparison.csv")
tabS_pk <- rbind(
  data.frame(block = "structure summary", item = pss$structure,
             detail = sprintf("WSS %.3f; %d parameters; residual df %s; Hessian PD %s; selected %s", pss$wss, pss$n_parameters,
                              if ("residual_df" %in% names(pss)) pss$residual_df else "NA", if ("hessian_pd" %in% names(pss)) pss$hessian_pd else "NA", pss$selected),
             observed = NA_real_, one_compartment = NA_real_, two_compartment = NA_real_, stringsAsFactors = FALSE),
  data.frame(block = "row-level fit (all training rows, including C24)", item = sprintf("%s %s %s %s", psc$study, psc$regimen, psc$analyte, psc$stat),
             detail = "", observed = psc$mean, one_compartment = psc$pred_one_compartment, two_compartment = psc$pred_two_compartment),
  data.frame(block = "PK sensitivity fits", item = pse$label,
             detail = sprintf("%s; fp at 100 mg %.2f; WSS %.2f; parent terminal half-life %.1f h", pse$structure, pse$fp, pse$wss, pse$parent_t_half_terminal_h),
             observed = NA_real_, one_compartment = NA_real_, two_compartment = NA_real_),
  # R4 fix round 1: the dose-model comparison (linear, fp(D), parent-only loss), with its D50 and Bested-CV sensitivities
  data.frame(block = "dose-model comparison (pk_dose_model_comparison.csv)", item = sprintf("%s | %s", pdm$analysis, pdm$dose_model),
             detail = sprintf("WSS %.2f (Verebey %.2f, Bested %.2f); dWSS vs linear %.2f; %s; Hessian PD %s; worst Verebey row %.2f-fold of linear; bioavailability %.4f at 4.5 mg, %.4f at 100 mg%s%s",
                              pdm$wss, pdm$wss_verebey, pdm$wss_bested, pdm$delta_wss,
                              ifelse(is.na(pdm$extra_parameter), "no extra parameter", sprintf("%s %.4g (SE logit %.3g)", pdm$extra_parameter, pdm$extra_estimate, pdm$extra_se_logit)),
                              pdm$hessian_pd, pdm$max_verebey_worsening, pdm$bioavailability_4.5mg, pdm$bioavailability_100mg,
                              ifelse(as.logical(pdm$at_boundary) %in% TRUE, "; BOUNDARY fit", ""), ifelse(as.logical(pdm$adopted) %in% TRUE, "; ADOPTED", "")),
             observed = NA_real_, one_compartment = NA_real_, two_compartment = NA_real_))
wr(tabS_pk, "tableS_pk_structure.csv")

# ---- Supplementary: Bruun-Plesner check (NA where not estimable, with the notes, never a number in its place) ----
bru <- rr("bruun_check.csv"); need_cols(bru, c("hypothesis", "model_ED50_mg", "model_ED95_mg", "observed_ED50_mg", "observed_ED95_mg"), "bruun_check.csv")
bru$model_ED50_display <- ifelse(is.finite(bru$model_ED50_mg), sprintf("%.2f", bru$model_ED50_mg), "NA (not estimable)")
bru$model_ED95_display <- ifelse(is.finite(bru$model_ED95_mg), sprintf("%.2f", bru$model_ED95_mg), "NA (not estimable)")
wr(bru, "tableS_bruun_check.csv")

# ---- Supplementary: Sobol indices with the two-seed agreement (S1 by Sobol 2007, ST by Jansen; Jansen's S1 never shown) ----
sb <- rr("sobol.csv"); need_cols(sb, c("hypothesis", "day", "parameter", "S1", "S1_lo", "S1_hi", "ST", "ST_lo", "ST_hi", "eligible"), "sobol.csv")
sa <- rr("sobol_seed_agreement.csv"); need_cols(sa, c("hypothesis", "day", "parameter", "index", "abs_diff", "max_ci_half_width", "disagree"), "sobol_seed_agreement.csv")
sb <- sb[as.logical(sb$eligible) %in% TRUE & !is.na(sb$parameter), c("hypothesis", "day", "parameter", "S1", "S1_lo", "S1_hi", "ST", "ST_lo", "ST_hi")]
for (ix in c("S1", "ST")) {
  a <- sa[sa$index == ix, ]; est2 <- grep("^estimate_seed", names(a), value = TRUE)
  k <- match(paste(sb$hypothesis, sb$day, sb$parameter), paste(a$hypothesis, a$day, a$parameter))
  if (anyNA(k)) stop("sobol_seed_agreement.csv lacks ", ix, " rows for some indices in sobol.csv", call. = FALSE)
  if (length(est2) == 2) sb[[paste0(ix, "_second_seed")]] <- a[[est2[2]]][k]
  sb[[paste0(ix, "_seed_abs_diff")]] <- a$abs_diff[k]; sb[[paste0(ix, "_seeds_disagree")]] <- as.logical(a$disagree[k])
}
sb$note <- "S1: Sobol 2007 estimator; ST: Jansen; same model runs (Saltelli et al. 2010). seeds_disagree: the two seeds differ by more than the larger bootstrap 95% CI half-width."
wr(sb, "tableS_sobol.csv")
cat("Tables written to", res, "\n")
