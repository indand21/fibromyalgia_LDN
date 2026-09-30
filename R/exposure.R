exposure_profile <- function(dose_mg, pkp, bp, tau = 24, n_doses = 14, step = 0.25) {
  hours <- seq(0, tau, by = step)
  if (dose_mg == 0) return(data.frame(hour = hours, Cu_NTX = 0, Cu_BN = 0, RO = 0))
  pk <- simulate_pk(dose_mg, seq(0, tau * n_doses, by = step), pkp, tau, n_doses)
  occ <- simulate_occupancy(pk, bp)
  t_last <- tau * (n_doses - 1); t_prev <- tau * (n_doses - 2)
  last <- pk$time >= t_last & pk$time <= t_last + tau
  prev <- pk$time >= t_prev & pk$time <= t_prev + tau
  if (abs(mean(occ$RO[last]) - mean(occ$RO[prev])) > 0.01 * max(mean(occ$RO[last]), 1e-6))
    warning("exposure_profile: occupancy not at steady state after ", n_doses, " doses (dose ", dose_mg, " mg)")
  data.frame(hour = pk$time[last] - t_last, Cu_NTX = pk$Cu_NTX[last], Cu_BN = pk$Cu_BN[last], RO = occ$RO[last])
}

exposure_cache <- function(doses, pkp, bp, ...) {
  doses <- sort(unique(doses))
  stats::setNames(lapply(doses, exposure_profile, pkp = pkp, bp = bp, ...), as.character(doses))
}

# Means over one dosing interval [0, 24): the endpoint duplicates hour 0 at steady state.
profile_means <- function(prof, gamma, IC50_TLR4, Kd_OGFr) {
  w <- prof$hour < max(prof$hour)
  list(mRO = mean(prof$RO[w]),
       free_term = mean((1 - prof$RO[w])^gamma),
       tlr4 = mean(prof$Cu_NTX[w] / (IC50_TLR4 + prof$Cu_NTX[w])),
       ogfr = mean(prof$Cu_NTX[w] / (Kd_OGFr + prof$Cu_NTX[w])))
}

# Steady-state occupancy summary (mean, peak, trough over one dosing interval) per cached dose.
steady_state_occupancy <- function(cache)
  do.call(rbind, lapply(names(cache), function(d) data.frame(dose_mg = as.numeric(d),
    mean_RO = mean(cache[[d]]$RO), peak_RO = max(cache[[d]]$RO), trough_RO = min(cache[[d]]$RO))))

# Occupancy table for the ACTIVE binding stage (binding_stage.txt), always computed from that stage's exposure cache.
# Never falls back to another stage: a missing stage file or cache is an error.
active_occupancy <- function(res = file.path(v5_root(), "output", "results")) {
  sf <- file.path(res, "binding_stage.txt")
  if (!file.exists(sf)) stop("binding_stage.txt missing: cannot determine the active binding stage")
  stage <- readLines(sf)[1]
  cf <- file.path(res, sprintf("exposure_cache_%s.rds", stage))
  if (!file.exists(cf)) stop("exposure cache for active stage ", stage, " missing: ", cf)
  out <- steady_state_occupancy(readRDS(cf)); out$stage <- stage
  out[order(out$dose_mg), ]
}
