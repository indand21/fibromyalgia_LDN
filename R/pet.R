# PET occupancy predictions for the actual regimen and observation time of each row (shared by scripts 03 and 10).
# Rows with reference == "normalised to 1 h scan" (Lee 1988) are percent blockade relative to that study's own 1 h scan,
# so the comparable prediction is 100 * RO(t) / RO(1 h) for the same dose and regimen (1 h after the last dose).
# "absolute" rows are 100 * RO(t).
# The number of doses given before the scan comes from the `n_doses` column of pet_occupancy.csv (R1 fix round 4: the only
# daily row, Weerts 2008, is 50 mg/day on days 15-19 with the scan on day 18, i.e. 4 doses, not the 7 assumed before).
# When that column is absent -- figure grids and occupancy_ed50() build rows themselves -- PET_DEFAULT_N_DOSES applies.
# Rows sharing a dose, regimen and dose count share one simulation (to the latest time needed); the result per row is the
# same as simulating each row on its own (tested), it is only faster, which matters inside the stage-J optimiser.
PET_NORM_REF <- "normalised to 1 h scan"
PET_DEFAULT_N_DOSES <- c(single = 1, daily = 7)
predict_pet <- function(bp, rows, pk_par) {
  out <- rep(NA_real_, nrow(rows))
  ref <- if (is.null(rows$reference)) rep("absolute", nrow(rows)) else rows$reference
  nd <- if (is.null(rows$n_doses)) unname(PET_DEFAULT_N_DOSES[rows$regimen]) else as.integer(rows$n_doses)
  if (anyNA(nd)) stop("predict_pet: cannot determine the number of doses (regimen must be ",
                      paste(names(PET_DEFAULT_N_DOSES), collapse = "/"), " or n_doses given)")
  key <- paste(rows$dose_mg, rows$regimen, nd, sep = "|")
  for (k in unique(key)) {
    idx <- which(key == k); r <- rows[idx[1], ]; n <- nd[idx[1]]
    t_obs <- 24 * (n - 1) + rows$time_h[idx]
    t_ref <- 24 * (n - 1) + 1
    s <- simulate_pk(r$dose_mg, seq(0, max(t_obs, t_ref) + 1, 0.25), pk_par, 24, n)
    o <- simulate_occupancy(s, bp)
    ro_t <- stats::approx(o$time, o$RO, t_obs)$y
    ro_ref <- stats::approx(o$time, o$RO, t_ref)$y
    out[idx] <- ifelse(!is.na(ref[idx]) & ref[idx] == PET_NORM_REF, 100 * ro_t / ro_ref, 100 * ro_t)
  }
  out
}
