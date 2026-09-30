# Between-subject variability, virtual trial replicates and responder rates (Task 9).

# Simulate n individuals in one arm. pmax_trial may be a scalar or a length-n vector (one placebo ceiling per subject).
# arm_shift (scalar or length n) is an arm-level model-structure residual added to relative pain before scaling by pain0.
simulate_individuals <- function(hyp, th, schedule, day, cache, pmax_trial, n, omega_p, omega_e,
                                 base_mean, base_sd, sd_res, seed = NULL, arm_shift = 0) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(pmax_trial) %in% c(1, n), length(arm_shift) %in% c(1, n))
  pmax_i <- inv_logit(logit(pmax_trial) + stats::rnorm(n, 0, omega_p))
  emax <- if (is.null(th$Emax)) 0 else th$Emax
  emax_i <- if (emax > 0) inv_logit(logit(emax) + stats::rnorm(n, 0, omega_e)) else rep(0, n)
  rel <- simulate_arm_vec(hyp, th, schedule, day, cache, pmax_i, emax_i)[day + 1, ]
  pain0 <- pmin(pmax(stats::rnorm(n, base_mean, base_sd), 4), 10)   # FM trials enrol NRS >= 4
  pain_end <- pmin(pmax(pain0 * (rel + arm_shift) + stats::rnorm(n, 0, sd_res), 0), 10)
  data.frame(pain0 = pain0, pain_end = pain_end)
}

responder_rates <- function(pain0, pain_end, cuts = c(0.15, 0.3, 0.5)) {
  red <- (pain0 - pain_end) / pain0
  stats::setNames(vapply(cuts, function(c) mean(red >= c), 0), paste0("resp", round(100 * cuts)))
}

# Moment-match the SD of % pain reduction in each Younger 2013 condition (training data) at the posterior median.
# `arms` is the training table; only rows of `trial_pattern` with a non-NA sd_pct are used (no test data).
calibrate_omegas <- function(fit, arms, cache, sd_res = 0.8, n = 4000, trial_pattern = "Younger.*2013") {
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
  sd_pct <- function(arm, op, oe) {
    x <- simulate_individuals(fit$hyp, th, parse_schedule(arm$schedule), arm$day, cache, pm, n, op, oe,
                              arm$baseline_nrs, bsd(arm), sd_res, seed = 11)
    stats::sd(100 * (x$pain0 - x$pain_end) / x$pain0)
  }
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
  list(omega_p = op, omega_e = oe, sd_res = sd_res,
       achieved_sd_pbo = ach_p, target_sd_pbo = pbo$sd_pct, achieved_sd_ldn = ach_l, target_sd_ldn = tgt_l,
       at_bound = at_bound, scale_used = scale_used, base_sd_source = paste0("pbo:", src(pbo), if (!is.null(ldn)) paste0("; ldn:", src(ldn))))
}

# design: list(trial, schedule_ldn, day, n_ldn, n_pbo, base_mean, base_sd)
# Each replicate arm also gets one arm-level draw z ~ N(0,1) shifting its relative pain by z * sigma_rel
# (the model-structure residual, exp(log_sigma_rel)), so replicate diffs carry the same residual the fit was given.
virtual_trials <- function(fit, design, cache, om, n_rep = 1000, n_draws = 200, seed = 99) {
  set.seed(seed); idx <- round(seq(1, nrow(fit$draws), length.out = n_draws))
  do.call(rbind, lapply(idx, function(i) {
    par <- fit$draws[i, ]; th <- theta_from(par, fit$hyp, fit$fixed)
    pm <- inv_logit(par[["mu_p"]] + stats::rnorm(n_rep, 0, exp(par[["log_tau"]])))
    sig_rel <- if ("log_sigma_rel" %in% names(par)) exp(par[["log_sigma_rel"]]) else 0
    one <- function(schedule, n_arm) {
      # all replicates in ONE vectorised call: subject j belongs to replicate rep_id[j] with trial placebo pm[rep_id[j]]
      rep_id <- rep(seq_len(n_rep), each = n_arm)
      shift <- rep(stats::rnorm(n_rep) * sig_rel, each = n_arm)
      x <- simulate_individuals(fit$hyp, th, parse_schedule(schedule), design$day, cache, pm[rep_id],
             n_rep * n_arm, om$omega_p, om$omega_e, design$base_mean, design$base_sd, om$sd_res, arm_shift = shift)
      red <- (x$pain0 - x$pain_end) / x$pain0
      chg <- x$pain_end - x$pain0
      by_rep <- function(v, f = mean) as.vector(tapply(v, rep_id, f))
      list(change = by_rep(chg), sd = by_rep(chg, stats::sd),
           resp = cbind(resp15 = by_rep(red >= 0.15), resp30 = by_rep(red >= 0.30), resp50 = by_rep(red >= 0.50)))
    }
    a <- one(design$schedule_ldn, design$n_ldn); b <- one("0:0", design$n_pbo)
    data.frame(draw = i, rep = seq_len(n_rep), diff = a$change - b$change, sd_ldn = a$sd, sd_pbo = b$sd,
               resp30_ldn = a$resp[, "resp30"], resp30_pbo = b$resp[, "resp30"],
               resp15_ldn = a$resp[, "resp15"], resp15_pbo = b$resp[, "resp15"],
               resp50_ldn = a$resp[, "resp50"], resp50_pbo = b$resp[, "resp50"])
  }))
}
