# Binding-stage choice (R1, 2026-10-01; replaces the run-5 A-versus-A2 held-out RMSE rule, which the expert review found
# wrong: an UNWEIGHTED RMSE on five non-calibration rows, four of them Lee 1988 duration points with SEM-only dispersion,
# selected for kinetics and discarded potency, where every LDN conclusion depends on potency at 1.5-6 mg; review C1, M8).
#
# The documented choice is now:
#   * stage A (bottom-up, in vitro constants) is scored against all 11 PET rows as a genuine held-out test;
#   * stage A2 (one Kd scale fitted to the six Rabiner 2011 rows) is kept as a sensitivity bound;
#   * stage J (joint Kd scale AND koff scale fitted to all 11 rows, precision-weighted, Lee rows on their normalised
#     estimand) is ADOPTED. It uses the PET data for calibration after bottom-up prediction failed: post-hoc, exploratory.
rmse <- function(obs, pred) sqrt(mean((pred - obs)^2))

# keq (brain equilibration rate, 1/h) is not identifiable from the available data; it is fixed at `default` unless another
# value lowers the weighted SSE by more than `threshold` (3.84 = chi-square(1) 95% point; a likelihood-ratio heuristic for
# one parameter, valid because the WSSE is -2 log-likelihood up to a constant).
select_keq <- function(tab, default = 2, threshold = stats::qchisq(0.95, 1)) {
  stopifnot(all(c("keq", "wsse") %in% names(tab)), default %in% tab$keq)
  w_def <- tab$wsse[tab$keq == default]
  i <- which.min(tab$wsse); gain <- w_def - tab$wsse[i]
  if (tab$keq[i] != default && gain > threshold)
    list(keq = tab$keq[i], gain = gain,
         reason = sprintf("keq = %g per h lowers the weighted SSE by %.2f (> %.2f, the chi-square(1) 95%% point) relative to the default %g per h: adopted.",
                          tab$keq[i], gain, threshold, default))
  else
    list(keq = default, gain = gain,
         reason = sprintf("no keq value lowers the weighted SSE by more than %.2f (the chi-square(1) 95%% point) relative to the default %g per h (largest gain %.2f, at keq = %g): default kept.",
                          threshold, default, gain, tab$keq[i]))
}

# One-row record of the binding-stage choice. A, A2, J are lists of summaries; keq_tab has columns keq and wsse.
binding_stage_decision <- function(A, A2, J, keq_tab, keq_choice, ed50_obs) {
  k <- stats::setNames(as.list(keq_tab$wsse), paste0("keq_wsse_", keq_tab$keq))
  reason <- paste0(
    sprintf("Stage A (bottom-up, no PET row fitted) is the genuine held-out PET test and fails it: RMSE %.1f points on all %d rows (weighted SSE %.1f), implied single-dose occupancy ED50 %.2f mg vs Rabiner 2011 %.2f mg (%.2f-%.2f). ",
            A$rmse, A$n, A$wsse, A$ed50, ed50_obs[["value"]], ed50_obs[["lo95"]], ed50_obs[["hi95"]]),
    sprintf("Stage A2 (one Kd scale %.2f fitted to the six Rabiner rows only) fits those rows (RMSE %.1f) but not the Lee 1988 duration rows (all-row RMSE %.1f, weighted SSE %.1f), because its koff is unchanged; kept as a sensitivity bound. ",
            A2$kd_scale, A2$rmse_calib, A2$rmse_all, A2$wsse),
    sprintf("Stage J fits a shared Kd scale (%.2f, SE of log %.3f) and a shared koff scale (%.4f, SE of log %.3f) jointly to all 11 rows, precision-weighted (weighted SSE %.1f, RMSE %.1f; Hessian %s). ",
            J$kd_scale, J$kd_scale_se_log, J$koff_scale, J$koff_scale_se_log, J$wsse, J$rmse,
            if (isTRUE(J$hessian_pd)) "positive definite" else "NOT positive definite"),
    sprintf("keq sensitivity: %s ", keq_choice$reason),
    "Stage J is adopted because it is the only stage that fits the potency (Rabiner) and duration (Lee) data jointly; per-row residuals are in pet_test.csv. ",
    "It uses the PET data for calibration after bottom-up prediction failed, so it is post-hoc and exploratory, and PET is no longer a held-out test for the adopted stage; ",
    "the Rabiner ED50 comparison is not independent for A2 or J (it summarises scans used to calibrate them).")
  data.frame(stage_A_rmse_heldout = A$rmse, stage_A_wsse_heldout = A$wsse, n_heldout_A = A$n, stage_A_ed50_mg = A$ed50,
             stage_A2_kd_scale = A2$kd_scale, stage_A2_rmse_calibration = A2$rmse_calib, stage_A2_rmse_all = A2$rmse_all,
             stage_A2_wsse_all = A2$wsse, stage_A2_ed50_mg = A2$ed50,
             stage_J_kd_scale = J$kd_scale, stage_J_log_kd_scale_se = J$kd_scale_se_log,
             stage_J_koff_scale = J$koff_scale, stage_J_log_koff_scale_se = J$koff_scale_se_log,
             stage_J_wsse = J$wsse, stage_J_rmse_all = J$rmse, stage_J_hessian_pd = isTRUE(J$hessian_pd), stage_J_ed50_mg = J$ed50,
             as.data.frame(k, check.names = FALSE),
             keq_adopted = keq_choice$keq, keq_reason = keq_choice$reason, lee_tail_identifiability = LEE_TAIL_NOTE,
             rabiner_ed50_mg = ed50_obs[["value"]], rabiner_ed50_lo95 = ed50_obs[["lo95"]], rabiner_ed50_hi95 = ed50_obs[["hi95"]],
             stage_adopted = "J", reason = reason, stringsAsFactors = FALSE, check.names = FALSE)
}
