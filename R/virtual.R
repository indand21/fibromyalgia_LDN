# Between-subject variability, virtual trial replicates and responder rates (Task 9).

# Top of the 0-10 NRS. Each trial's LOWER baseline bound is its own documented minimum-pain eligibility criterion
# (0 when none applies), curated with sources in data/params/baseline_bounds.csv (R2 fix round 1).
NRS_MAX <- 10
DEFAULT_BASE_BOUNDS <- c(0, NRS_MAX)

read_baseline_bounds <- function(file = file.path(v5_root(), "data", "params", "baseline_bounds.csv")) {
  b <- read_dataset(file)
  need <- c("trial", "lower_nrs", "upper_nrs")
  miss <- setdiff(need, names(b))
  if (length(miss)) stop(basename(file), " lacks column(s): ", paste(miss, collapse = ", "))
  bad <- !is.finite(b$lower_nrs) | !is.finite(b$upper_nrs) | b$lower_nrs < 0 | b$upper_nrs > NRS_MAX |
         b$lower_nrs >= b$upper_nrs
  if (any(bad)) stop(basename(file), ": lower_nrs/upper_nrs must satisfy 0 <= lower_nrs < upper_nrs <= ", NRS_MAX,
                     " (rows: ", paste(b$trial[bad], collapse = ", "), ")")
  if (anyDuplicated(tolower(b$trial))) stop(basename(file), ": duplicated trial rows")
  b
}

baseline_bounds_for <- function(bounds, trial) {
  i <- which(tolower(bounds$trial) == tolower(trial))
  if (length(i) != 1L) stop("no baseline bounds curated for trial '", trial, "' (data/params/baseline_bounds.csv)")
  c(bounds$lower_nrs[i], bounds$upper_nrs[i])
}

# The parent of a baseline distribution is solved once per (mean, sd, bounds) and remembered: calibrate_omegas()
# and virtual_trials() call simulate_individuals() many times with the same baseline.
.trunc_parent_memo <- new.env(parent = emptyenv())
baseline_parent <- function(mean, sd, bounds) {
  key <- sprintf("%.12g|%.12g|%.12g|%.12g", mean, sd, bounds[1], bounds[2])
  if (is.null(.trunc_parent_memo[[key]])) .trunc_parent_memo[[key]] <- match_trunc_parent(mean, sd, bounds[1], bounds[2])
  .trunc_parent_memo[[key]]
}

# Simulate n individuals in one arm. pmax_trial may be a scalar or a length-n vector (one placebo ceiling per subject).
# arm_shift (scalar or length n) is the model-structure residual added to relative pain before scaling by pain0.
# It is supplied by the caller and is a TRIAL-level term (review I7): virtual_trials() passes the SAME value to
# both arms of a replicate trial, so it inflates each arm's change but cancels in their difference.
# base_mean/base_sd are the trial's REPORTED baseline moments and base_bounds its enrolment bounds.
simulate_individuals <- function(hyp, th, schedule, day, cache, pmax_trial, n, omega_p, omega_e,
                                 base_mean, base_sd, sd_res, seed = NULL, arm_shift = 0,
                                 base_bounds = DEFAULT_BASE_BOUNDS) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(pmax_trial) %in% c(1, n), length(arm_shift) %in% c(1, n), length(base_bounds) == 2)
  pmax_i <- inv_logit(logit(pmax_trial) + stats::rnorm(n, 0, omega_p))
  emax <- if (is.null(th$Emax)) 0 else th$Emax
  emax_i <- if (emax > 0) inv_logit(logit(emax) + stats::rnorm(n, 0, omega_e)) else rep(0, n)
  rel <- simulate_arm_vec(hyp, th, schedule, day, cache, pmax_i, emax_i)[day + 1, ]
  # Baselines are TRUNCATED to the trial's bounds, never censored onto them (review I3: censoring put 23.7% of
  # baselines at exactly 4.0). The reported moments describe the ENROLLED population, which has already passed
  # the screen, so the parent normal is moment-matched to make the TRUNCATED distribution reproduce them (R2 fix
  # round 1); truncating N(reported) instead applied the screen twice (Younger 5.08/1.50 -> 5.68/1.10).
  bp <- baseline_parent(base_mean, base_sd, base_bounds)
  pain0 <- rtnorm_trunc(n, bp$mu, bp$sigma, base_bounds[1], base_bounds[2])
  # The endpoint clip is PHYSICAL (the NRS cannot leave [0, 10]), so it stays, but it is reported rather than
  # silent: attr(, "pain_end_clip") gives the fraction of subjects it binds for at each end.
  e <- stats::rnorm(n, 0, sd_res)                   # measurement noise (same draw order as before)
  pain_raw <- pain0 * (rel + arm_shift) + e
  pain_end <- pmin(pmax(pain_raw, 0), NRS_MAX)
  out <- data.frame(pain0 = pain0, pain_end = pain_end)
  attr(out, "pain_end_clip") <- c(low = mean(pain_raw < 0), high = mean(pain_raw > NRS_MAX))
  attr(out, "noise_term") <- e                      # the measurement-noise draws, for diagnostics (R3a item 7)
  attr(out, "baseline_parent") <- bp
  out
}

responder_rates <- function(pain0, pain_end, cuts = c(0.15, 0.3, 0.5)) {
  red <- (pain0 - pain_end) / pain0
  stats::setNames(vapply(cuts, function(c) mean(red >= c), 0), paste0("resp", round(100 * cuts)))
}

# Moment-match the SD of % pain reduction in each Younger 2013 condition (training data) at the posterior median.
# `arms` is the training table; only rows of `trial_pattern` with a non-NA sd_pct are used (no test data).
# `bounds` is the curated per-trial baseline-bounds table (read_baseline_bounds()); without it the physical
# bounds [0, 10] are used; the curated table records [0.5, 10] for Younger 2013 (a floor chosen by a seed-stability rule).
# `seed` fixes the simulated cohort (default 11, as before); the Younger floor in baseline_bounds.csv was chosen so
# that omega_p does not depend on it (R2 fix round 2).
calibrate_omegas <- function(fit, arms, cache, sd_res = 0.8, n = 4000, trial_pattern = "Younger.*2013",
                             bounds = NULL, seed = 11) {
  if (!"sd_pct" %in% names(arms)) stop("training data lacks an sd_pct column; curate the SD of % reduction for Younger 2013")
  arms <- arms[grepl(trial_pattern, arms$trial, ignore.case = TRUE) & !is.na(arms$sd_pct), ]
  pbo <- arms[arms$arm == "placebo", ]; ldn <- arms[arms$arm != "placebo", ]
  if (nrow(pbo) == 0 || (fit$hyp != "H0" && nrow(ldn) == 0))
    stop("sd_pct is missing for the needed Younger 2013 arms (placebo", if (fit$hyp != "H0") " and LDN", ")")
  if (nrow(pbo) > 1) warning("more than one Younger 2013 placebo row qualifies; using the first")
  if (fit$hyp != "H0" && nrow(ldn) > 1) warning("more than one Younger 2013 LDN row qualifies; using the first")
  pbo <- pbo[1, ]; ldn <- if (nrow(ldn)) ldn[1, ] else NULL
  # Put baseline_nrs on the 0-10 NRS scale used by sd_res_nrs and the simulated NRS (0-100 VAS -> /10).
  if (!"scale" %in% names(arms)) warning("efficacy_training.csv has no `scale` column; assuming baseline_nrs is on the 0-10 scale")
  scale_of <- function(arm) if ("scale" %in% names(arm) && !is.na(arm$scale)) arm$scale else "0-10"
  rescale <- function(arm) { sc <- scale_of(arm)
    if (identical(sc, "0-100")) { arm$baseline_nrs <- arm$baseline_nrs / 10
      if ("baseline_sd" %in% names(arm)) arm$baseline_sd <- arm$baseline_sd / 10 }
    else if (!identical(sc, "0-10")) stop("unknown scale '", sc, "' in efficacy_training.csv (expected 0-100 or 0-10)")
    arm }
  scale_used <- paste0("pbo:", scale_of(pbo), if (!is.null(ldn)) paste0("; ldn:", scale_of(ldn)), " (0-100 baseline divided by 10)")
  pbo <- rescale(pbo); if (!is.null(ldn)) ldn <- rescale(ldn)
  # baseline SD: from the arm's baseline_sd when present and finite, else an assumed 1.5
  bsd <- function(arm) if ("baseline_sd" %in% names(arm) && is.finite(arm$baseline_sd)) arm$baseline_sd else 1.5
  src <- function(arm) if ("baseline_sd" %in% names(arm) && is.finite(arm$baseline_sd)) "data" else "assumed 1.5"
  med <- apply(fit$draws, 2, stats::median); th <- theta_from(med, fit$hyp, fit$fixed)
  eta <- paste0("eta_", pbo$trial)
  if (!eta %in% names(med)) stop("posterior has no ", eta)
  pm <- inv_logit(med[["mu_p"]] + med[[eta]])
  bb <- if (is.null(bounds)) DEFAULT_BASE_BOUNDS else baseline_bounds_for(bounds, pbo$trial)
  bb_src <- if (is.null(bounds)) "default physical bounds" else "data/params/baseline_bounds.csv"
  sim_one <- function(arm, op, oe)
    simulate_individuals(fit$hyp, th, parse_schedule(arm$schedule), arm$day, cache, pm, n, op, oe,
                         arm$baseline_nrs, bsd(arm), sd_res, seed = seed, base_bounds = bb)
  sd_pct <- function(arm, op, oe) { x <- sim_one(arm, op, oe); stats::sd(100 * (x$pain0 - x$pain_end) / x$pain0) }
  # How often the physical pain_end clip binds at the calibrated omegas (review I3: report, do not hide)
  clip_of <- function(arm, op, oe) attr(sim_one(arm, op, oe), "pain_end_clip")
  near_bound <- function(o) o < 0.01 || o > 2.99
  op <- stats::optimize(function(o) (sd_pct(pbo, o, 0) - pbo$sd_pct)^2, c(0, 3))$minimum
  ach_p <- sd_pct(pbo, op, 0)
  if (fit$hyp == "H0") { oe <- 0; ach_l <- NA_real_; tgt_l <- NA_real_; bnd_e <- FALSE }
  else {
    oe <- stats::optimize(function(o) (sd_pct(ldn, op, o) - ldn$sd_pct)^2, c(0, 3))$minimum
    ach_l <- sd_pct(ldn, op, oe); tgt_l <- ldn$sd_pct; bnd_e <- near_bound(oe)
  }
  bnd_p <- near_bound(op); at_bound <- bnd_p || bnd_e
  bad_p <- abs(ach_p - pbo$sd_pct) > 1; bad_l <- !is.na(ach_l) && abs(ach_l - tgt_l) > 1
  if (at_bound || bad_p || bad_l)
    warning(fit$hyp, ": omega calibration did not match the target SD of % reduction (placebo achieved ",
            round(ach_p, 2), " vs ", pbo$sd_pct, if (!is.na(ach_l)) paste0("; LDN achieved ", round(ach_l, 2), " vs ", tgt_l),
            "; omega_p=", round(op, 3), ", omega_e=", round(oe, 3), "; at_bound=", at_bound, ")")
  cl_p <- clip_of(pbo, op, 0); cl_l <- if (is.null(ldn)) c(low = NA_real_, high = NA_real_) else clip_of(ldn, op, oe)
  # The calibrated placebo cohort's own measurement-noise term 100 e / pain0 (R3a item 7: the Younger floor study
  # reads it from here instead of re-creating the cohort's RNG order by hand) and its smallest baseline.
  xp <- sim_one(pbo, op, 0)
  noise_sd <- stats::sd(100 * attr(xp, "noise_term") / xp$pain0)
  bp <- baseline_parent(pbo$baseline_nrs, bsd(pbo), bb)
  list(omega_p = op, omega_e = oe, sd_res = sd_res,
       achieved_sd_pbo = ach_p, target_sd_pbo = pbo$sd_pct, achieved_sd_ldn = ach_l, target_sd_ldn = tgt_l,
       at_bound = at_bound, scale_used = scale_used, base_sd_source = paste0("pbo:", src(pbo), if (!is.null(ldn)) paste0("; ldn:", src(ldn))),
       pain_end_clip_low_pbo = cl_p[["low"]], pain_end_clip_high_pbo = cl_p[["high"]],
       pain_end_clip_low_ldn = cl_l[["low"]], pain_end_clip_high_ldn = cl_l[["high"]],
       noise_sd_pbo = noise_sd, min_baseline_pbo = min(xp$pain0),
       baseline_bounds = sprintf("[%g, %g] (%s)", bb[1], bb[2], bb_src),
       baseline_parent_mu = bp$mu, baseline_parent_sigma = bp$sigma, baseline_parent_matched = bp$matched,
       baseline_sampling = paste0("truncated normal on the trial's bounds, parent moment-matched to the reported ",
                                  "(enrolled) mean and SD (review I3; R2 fix round 1): ", bp$note))
}

# design: list(trial, schedule_ldn, day, n_ldn, n_pbo, base_mean, base_sd[, base_bounds])
# base_mean/base_sd are the trial's REPORTED baseline moments; base_bounds its enrolment bounds from
# data/params/baseline_bounds.csv (default [0, 10]).
# The model-structure residual (exp(log_sigma_rel)) is a TRIAL-level term (review I7): each replicate trial gets
# ONE draw z ~ N(0, 1) per posterior draw, shifting the relative pain of BOTH of its arms by z * sigma_rel.
# A misspecified disease course therefore moves the whole replicate together: it inflates each arm's absolute
# change but cancels in the LDN-minus-placebo difference, instead of entering it at root 2 times sigma_rel.
# Replicates (and posterior draws) are independent of one another.
# Returned attributes: "pain_end_clip" (one row per posterior draw and arm: the fraction of simulated subjects
# whose end-of-trial NRS hit 0 or 10) and "baseline_parent" (the moment-matched baseline distribution used).
virtual_trials <- function(fit, design, cache, om, n_rep = 1000, n_draws = 200, seed = 99) {
  set.seed(seed); idx <- round(seq(1, nrow(fit$draws), length.out = n_draws))
  bb <- if (is.null(design$base_bounds)) DEFAULT_BASE_BOUNDS else design$base_bounds
  clips <- list()
  out <- do.call(rbind, lapply(idx, function(i) {
    par <- fit$draws[i, ]; th <- theta_from(par, fit$hyp, fit$fixed)
    pm <- inv_logit(par[["mu_p"]] + stats::rnorm(n_rep, 0, exp(par[["log_tau"]])))
    sig_rel <- if ("log_sigma_rel" %in% names(par)) exp(par[["log_sigma_rel"]]) else 0
    shift_rep <- stats::rnorm(n_rep) * sig_rel   # one shared residual per replicate trial, used by both arms
    one <- function(schedule, n_arm, arm) {
      # all replicates in ONE vectorised call: subject j belongs to replicate rep_id[j] with trial placebo pm[rep_id[j]]
      rep_id <- rep(seq_len(n_rep), each = n_arm)
      shift <- rep(shift_rep, each = n_arm)
      x <- simulate_individuals(fit$hyp, th, parse_schedule(schedule), design$day, cache, pm[rep_id],
             n_rep * n_arm, om$omega_p, om$omega_e, design$base_mean, design$base_sd, om$sd_res, arm_shift = shift,
             base_bounds = bb)
      cl <- attr(x, "pain_end_clip")
      clips[[length(clips) + 1L]] <<- data.frame(draw = i, arm = arm, low = cl[["low"]], high = cl[["high"]])
      red <- (x$pain0 - x$pain_end) / x$pain0
      chg <- x$pain_end - x$pain0
      by_rep <- function(v, f = mean) as.vector(tapply(v, rep_id, f))
      list(change = by_rep(chg), sd = by_rep(chg, stats::sd),
           resp = cbind(resp15 = by_rep(red >= 0.15), resp30 = by_rep(red >= 0.30), resp50 = by_rep(red >= 0.50)))
    }
    a <- one(design$schedule_ldn, design$n_ldn, "ldn"); b <- one("0:0", design$n_pbo, "placebo")
    data.frame(draw = i, rep = seq_len(n_rep), diff = a$change - b$change, sd_ldn = a$sd, sd_pbo = b$sd,
               resp30_ldn = a$resp[, "resp30"], resp30_pbo = b$resp[, "resp30"],
               resp15_ldn = a$resp[, "resp15"], resp15_pbo = b$resp[, "resp15"],
               resp50_ldn = a$resp[, "resp50"], resp50_pbo = b$resp[, "resp50"])
  }))
  attr(out, "pain_end_clip") <- do.call(rbind, clips)
  attr(out, "baseline_parent") <- baseline_parent(design$base_mean, design$base_sd, bb)
  out
}
