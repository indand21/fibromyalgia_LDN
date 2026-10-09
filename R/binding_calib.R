# Joint binding calibration (stage J; review C1, C4, M8; controller ruling 2026-10-01).
# Two scales are fitted jointly against ALL PET rows, on the log scale:
#   kd_scale   multiplies BOTH dissociation constants (Kd_NTX, Kd_BN)  -> potency (in vitro -> in vivo transfer)
#   koff_scale multiplies BOTH dissociation rates (koff_NTX, koff_BN)   -> duration (the Lee 1988 168 h tail)
# kon = koff / Kd is recomputed, so a Kd scale alone leaves koff unchanged and slows kon.
# PET thereby becomes CALIBRATION data for binding; this is post-hoc and exploratory (the bottom-up stage A, which was
# not fitted to any PET row, is the only genuine held-out PET test).

# Precision weights: the reported SD where finite, otherwise a documented default of 10 percentage points. The default
# applies to the Rabiner 2011 digitised single-subject points (no dispersion) and to Lee 1988 (dispersion reported as SEM
# with n = 3-5 only; SEM is not an SD of a single observation and is not converted here).
PET_DEFAULT_SD <- 10
pet_precision_sd <- function(sd, default = PET_DEFAULT_SD) ifelse(is.finite(sd), sd, default)

pet_wsse_from_pred <- function(pred, rows) sum(((pred - rows$occupancy_pct) / pet_precision_sd(rows$sd))^2)
pet_wsse <- function(bp, rows, pk_par) pet_wsse_from_pred(predict_pet(bp, rows, pk_par), rows)

scale_binding <- function(bp0, kd_scale = 1, koff_scale = 1, keq = bp0$keq) {
  b <- binding_params(Kd_NTX = bp0$Kd_NTX * kd_scale, Kd_BN = bp0$Kd_BN * kd_scale,
                      koff_NTX = bp0$koff_NTX * koff_scale, koff_BN = bp0$koff_BN * koff_scale, keq = keq)
  b$kd_scale <- kd_scale; b$koff_scale <- koff_scale
  b
}

J_PAR_NAMES <- c("log_kd_scale", "log_koff_scale")

# Minimise the precision-weighted SSE over (log kd_scale, log koff_scale) at a fixed keq. The WSSE is -2 log-likelihood for
# Gaussian errors with known SD, so the Laplace covariance of the log scales is solve(hessian / 2).
fit_binding_joint <- function(bp0, rows, pk_par, keq = 2, start = c(log(10), log(0.1))) {
  stopifnot(all(c("dose_mg", "regimen", "time_h", "reference", "occupancy_pct", "sd") %in% names(rows)))
  obj <- function(ls) {
    v <- pet_wsse(scale_binding(bp0, exp(ls[1]), exp(ls[2]), keq), rows, pk_par)
    if (!is.finite(v)) 1e10 else v
  }
  o1 <- stats::optim(start, obj, method = "Nelder-Mead", control = list(reltol = 1e-10, maxit = 500))
  o2 <- stats::optim(o1$par, obj, method = "BFGS", control = list(reltol = 1e-12, maxit = 200))
  best <- if (o2$value <= o1$value) o2 else o1
  ls <- stats::setNames(best$par, J_PAR_NAMES)
  H <- stats::optimHess(ls, obj); dimnames(H) <- list(J_PAR_NAMES, J_PAR_NAMES)
  ev <- eigen(H / 2, symmetric = TRUE)$values
  pd <- all(is.finite(ev)) && min(ev) > 0
  cov <- if (pd) solve(H / 2) else matrix(NA_real_, 2, 2, dimnames = dimnames(H))
  bp <- scale_binding(bp0, exp(ls[[1]]), exp(ls[[2]]), keq)
  pred <- predict_pet(bp, rows, pk_par)
  list(log_scales = ls, kd_scale = exp(ls[[1]]), koff_scale = exp(ls[[2]]), keq = keq,
       se_log = stats::setNames(sqrt(diag(cov)), J_PAR_NAMES), cov = cov, hessian = H, hessian_eigen = ev, hessian_pd = pd,
       wsse = pet_wsse_from_pred(pred, rows), rmse = sqrt(mean((pred - rows$occupancy_pct)^2)), pred = pred,
       sd_used = pet_precision_sd(rows$sd), convergence = best$convergence, bp = bp)
}

# Multivariate normal draws (Laplace approximation). A seed, if given, is applied locally (the caller's RNG state is kept).
# A seed pins the GENERATOR as well as the seed (R1 fix round 4): otherwise the draws depend on whatever RNGkind happens to
# be in force, which parallel::clusterSetRNGStream() can change (see rng_state() in R/parallel.R). The caller's RNG state,
# kind included, is restored on exit.
LAPLACE_RNG_KIND <- "Mersenne-Twister"
laplace_draws <- function(mode, cov, n, seed = NULL, kind = LAPLACE_RNG_KIND) {
  R <- tryCatch(chol(cov), error = function(e) stop("laplace_draws: covariance is not positive definite"))
  if (!is.null(seed)) {
    st <- rng_state()                                  # captured BEFORE set.seed(); restored on exit
    on.exit(rng_restore(st), add = TRUE)
    set.seed(seed, kind = kind)
  }
  z <- matrix(stats::rnorm(n * length(mode)), n, length(mode))
  d <- sweep(z %*% R, 2, mode, "+")
  colnames(d) <- names(mode)
  d
}

# Median and 95% interval of a set of draws. `n` is the number of NON-missing draws the quantiles are actually based on.
summarise_draws <- function(values, probs = c(0.5, 0.025, 0.975)) {
  ok <- is.finite(values)
  q <- stats::quantile(values[ok], probs, names = FALSE)
  c(median = q[1], lo95 = q[2], hi95 = q[3], n = sum(ok))
}

# Single-dose dose giving 50% absolute occupancy at time_h (the Rabiner 2011 ED50 estimand: scans under 8 h post dose).
occupancy_ed50 <- function(bp, pk_par, time_h = 4, lower = 1e-3, upper = 1e3) {
  f <- function(ld) predict_pet(bp, data.frame(dose_mg = exp(ld), regimen = "single", time_h = time_h, reference = "absolute"),
                                pk_par) - 50
  exp(stats::uniroot(f, log(c(lower, upper)), tol = 1e-8)$root)
}

# Mean steady-state occupancy (same definition as steady_state_occupancy(): mean over one 24 h interval of the exposure
# profile, dosing extended up to EXPOSURE_MAX_DOSES until steady state) at each dose. Used per Laplace draw; the count of
# doses still not at steady state (warnings) is returned as attribute "n_warn".
ss_mean_occupancy <- function(doses, pk_par, bp) {
  nw <- 0L
  v <- withCallingHandlers(
    vapply(doses, function(d) mean(exposure_profile(d, pk_par, bp, max_doses = EXPOSURE_MAX_DOSES)$RO), 0),
    warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
  structure(stats::setNames(v, paste0("dose_", doses)), n_warn = nw)
}

# Identifiability statement written into the stage-J outputs (fit output, stage_decision.csv, pet_test_note.txt).
LEE_TAIL_NOTE <- paste(
  "The Lee 1988 blockade tail (30% of the 1 h value at 168 h after 50 mg) can be explained either by slow MOR dissociation",
  "(the stage-J explanation: a shared koff scale) or by a long terminal plasma phase of naltrexone/6-beta-naltrexol. The",
  "adopted two-compartment PK has a terminal phase of about 10 h, well short of the apparent 96 h half-life Verebey 1976 reports",
  "from 24 to 72 h in its abstract (held out of the primary fit as pk_summary_sensitivity_terminal.csv, since the abstract",
  "does not say which dose it follows).",
  "Slow dissociation pushes steady-state occupancy at 1.5-6 mg UP: receptor-bound drug carries over between daily doses, so",
  "occupancy accumulates and the peak-trough swing flattens. A long terminal PK phase pushes it DOWN, towards the",
  "rapid-equilibrium (A2-like) values: the tail is then made by low persisting concentrations that occupy receptors after",
  "50 mg but are 10-30 times lower after 1.5-6 mg, so they add little occupancy at low doses.",
  "The terminal-phase sensitivity fit (R1 fix round 3) is the one direct test of this in these data, and it PARTLY separates",
  "them: given the 24-72 h values the PK reaches a 23 h terminal phase, the fitted koff scale becomes less",
  "extreme, the PET fit improves, and occupancy at 4.5 mg falls (occupancy_sensitivity_J.csv) -- the predicted direction.",
  "It is not decisive, because that PK costs a large increase in weighted SSE against the reported 10 h half-lives (which are",
  "a 0-24 h biexponential estimand, not the model's terminal slope), so the data are internally inconsistent on this point.",
  "Slow dissociation still carries most of the tail under every analysis: koff_NTX stays far below the in vitro",
  "value (koff_provenance.csv). Occupancy at 4.5 mg should therefore be read across both analyses, not from the primary alone.")
