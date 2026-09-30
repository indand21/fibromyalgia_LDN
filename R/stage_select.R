# Binding-stage selection (correction of the run-4 rule, see CHANGELOG 2026-09-30 run 5).
# Stage A2 (a Kd scale fitted to the calibration rows) is adopted for downstream use ONLY if it lowers RMSE on the
# NON-calibration (held-out) rows relative to stage A evaluated on those same rows. Otherwise stage A is kept.
rmse <- function(obs, pred) sqrt(mean((pred - obs)^2))

select_binding_stage <- function(obs_heldout, pred_A_heldout, pred_A2_heldout, obs_all, pred_A_all, pred_A2_all) {
  stopifnot(length(obs_heldout) == length(pred_A_heldout), length(obs_heldout) == length(pred_A2_heldout),
            length(obs_all) == length(pred_A_all), length(obs_all) == length(pred_A2_all), length(obs_heldout) > 0)
  rA <- rmse(obs_heldout, pred_A_heldout); rA2 <- rmse(obs_heldout, pred_A2_heldout)
  adopt <- if (rA2 < rA) "A2" else "A"
  reason <- if (adopt == "A2")
    sprintf("A2 lowers RMSE on the %d held-out rows (%.2f vs %.2f for stage A on the same rows): adopted.", length(obs_heldout), rA2, rA)
  else
    sprintf("A2 does not improve RMSE on the %d held-out (non-calibration) rows (%.2f vs %.2f for stage A on the same rows): stage A kept. A2 wins only on the rows it was fitted to.",
            length(obs_heldout), rA2, rA)
  data.frame(stage_A_rmse_heldout = rA, stage_A2_rmse_heldout = rA2,
             stage_A_rmse_all = rmse(obs_all, pred_A_all), stage_A2_rmse_all = rmse(obs_all, pred_A2_all),
             n_heldout = length(obs_heldout), n_all = length(obs_all), stage_adopted = adopt, reason = reason,
             stringsAsFactors = FALSE)
}
