# PET occupancy predictions for the actual regimen and observation time of each row (shared by scripts 03 and 10).
# Rows with reference == "normalised to 1 h scan" (Lee 1988) are percent blockade relative to that study's own 1 h scan,
# so the comparable prediction is 100 * RO(t) / RO(1 h) for the same dose and regimen (1 h after the last dose).
# "absolute" rows are 100 * RO(t).
PET_NORM_REF <- "normalised to 1 h scan"
predict_pet <- function(bp, rows, pk_par) vapply(seq_len(nrow(rows)), function(i) {
  r <- rows[i, ]; n <- if (r$regimen == "daily") 7 else 1
  t_obs <- 24 * (n - 1) + r$time_h
  t_ref <- 24 * (n - 1) + 1
  s <- simulate_pk(r$dose_mg, seq(0, max(t_obs, t_ref) + 1, 0.25), pk_par, 24, n)
  o <- simulate_occupancy(s, bp)
  ro_t <- stats::approx(o$time, o$RO, t_obs)$y
  if (identical(r$reference, PET_NORM_REF)) 100 * ro_t / stats::approx(o$time, o$RO, t_ref)$y else 100 * ro_t
}, 0)
