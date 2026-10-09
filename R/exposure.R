# R1 fix round 4 item 1. The old rule stopped when the mean occupancy of the last two intervals differed by <= 1%. Under a
# slow koff the interval means approach their limit geometrically with ratio r, so the REMAINING gap is delta * r / (1 - r),
# which for r ~ 0.85-0.9 is 6-10 times the last step: the rule stopped several percentage points short (at stage J, 0.729
# against a converged 0.746 at 4.5 mg). The criterion is now that extrapolated remaining gap, in ABSOLUTE occupancy.
#
# exposure_ss_gap(m) takes the means of three consecutive intervals, oldest first, and returns the extrapolated remaining
# gap: r = (m3 - m2) / (m2 - m1), gap = |m3 - m2| * r / (1 - r). Degenerate cases are guarded: r >= 1 (not contracting)
# returns Inf, so the run is extended or a warning is issued; a non-monotone or flat sequence (numerical noise near the
# limit) falls back to the last step itself, which is then compared with the tolerance in the usual way.
EXPOSURE_SS_TOL <- 0.001    # absolute occupancy (0.1 percentage points)
EXPOSURE_MAX_DOSES <- 512L  # generous cap: a slow-koff stage-J profile needs 30-60 doses, a 10x slower one more

exposure_ss_gap <- function(m) {
  d1 <- m[2] - m[1]; d2 <- m[3] - m[2]
  if (!all(is.finite(c(d1, d2)))) return(Inf)
  if (d2 == 0) return(0)
  if (d1 == 0 || sign(d1) != sign(d2)) return(abs(d2))
  r <- d2 / d1
  if (!is.finite(r) || r >= 1) return(Inf)
  abs(d2) * r / (1 - r)
}

# Steady-state profile over one dosing interval. The dosing run is extended until exposure_ss_gap() certifies that the mean
# occupancy is within `ss_tol` of its limit, up to max_doses (default: no extension, so existing callers are unchanged).
# The step policy uses the estimated ratio to jump straight to the dose count the geometric tail implies (with a 20% margin)
# instead of doubling blindly, which costs 2-4x less; if the ratio is degenerate it falls back to doubling. A warning is
# issued only when the gap is still not certified at max_doses. Attributes: "n_doses" used, and the certified "ss_gap".
exposure_profile <- function(dose_mg, pkp, bp, tau = 24, n_doses = 14, step = 0.25, max_doses = n_doses,
                             ss_tol = EXPOSURE_SS_TOL) {
  hours <- seq(0, tau, by = step)
  if (dose_mg == 0)
    return(structure(data.frame(hour = hours, Cu_NTX = 0, Cu_BN = 0, RO = 0), n_doses = 0L, ss_gap = 0))
  n_doses <- max(as.integer(n_doses), 4L)
  repeat {
    pk <- simulate_pk(dose_mg, seq(0, tau * n_doses, by = step), pkp, tau, n_doses)
    occ <- simulate_occupancy(pk, bp)
    im <- function(k) { w <- pk$time >= tau * (n_doses - k) & pk$time <= tau * (n_doses - k + 1); mean(occ$RO[w]) }
    m <- c(im(3), im(2), im(1))
    gap <- exposure_ss_gap(m)
    if (gap <= ss_tol || n_doses >= max_doses) break
    r <- (m[3] - m[2]) / (m[2] - m[1])
    extra <- if (is.finite(r) && r > 0 && r < 1) ceiling(1.2 * log(ss_tol / gap) / log(r)) else n_doses
    n_doses <- as.integer(min(as.integer(max(n_doses + 2L, n_doses + extra)), max_doses))
  }
  if (gap > ss_tol)
    warning("exposure_profile: occupancy not at steady state after ", n_doses, " doses (dose ", dose_mg,
            " mg; extrapolated remaining gap ", signif(gap, 3), ")")
  last <- pk$time >= tau * (n_doses - 1) & pk$time <= tau * n_doses
  structure(data.frame(hour = pk$time[last] - tau * (n_doses - 1), Cu_NTX = pk$Cu_NTX[last], Cu_BN = pk$Cu_BN[last],
                       RO = occ$RO[last]),
            n_doses = as.integer(n_doses), ss_gap = as.numeric(gap))
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
