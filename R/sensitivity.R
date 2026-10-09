# Effect-guarded Sobol analysis and power, and the task runners that let scripts 07 and 09 send their independent
# simulations to a PSOCK cluster (R3a item 2, review I2).
#
# Reproducibility across worker counts. Every task carries everything it needs (a posterior THINNED to the draws it
# uses, the design, the omegas, the exposure cache and its own seed) and is run by a function defined in R/, which
# with_cluster() sources on every worker. A virtual-trial task installs its own Mersenne-Twister stream from its
# seed (with_mt_seed; virtual_trials() then calls set.seed(seed) itself), so its numbers depend only on the task,
# never on the worker, the number of workers, load balancing or the caller's RNG kind. A Sobol evaluation draws
# no random numbers at all. `map` is lapply by default; a cluster map is
#   function(X, FUN) parallel::parLapplyLB(cl, X, FUN)
# and must give results identical to lapply (tested).

# The rows of fit$draws that virtual_trials()/predict_trial() use for n_draws, and nothing else, so a task does not
# ship the whole posterior to a worker. virtual_trials() on the thinned fit picks exactly the same draws, because
# round(seq(1, n, length.out = n)) is 1:n. The original row numbers are kept in fit$draw_index.
thin_fit <- function(fit, n_draws) {
  idx <- round(seq(1, nrow(fit$draws), length.out = n_draws))
  fit$draws <- fit$draws[idx, , drop = FALSE]
  fit$draw_index <- if (is.null(fit$draw_index)) idx else fit$draw_index[idx]
  fit
}

# sigma_rel = 0 (exp(-Inf)): the sampling-variability-only variant of the power curve.
no_shift <- function(f) { if ("log_sigma_rel" %in% colnames(f$draws)) f$draws[, "log_sigma_rel"] <- -Inf; f }

# One virtual-trial job. task: fit (thinned), design, cache, om, n_rep, n_draws, seed[, sampling_only][, keep].
# keep = "diff" returns only the replicate differences (power), otherwise the full virtual_trials() output with
# its `draw` column mapped back to the original posterior row numbers.
run_vt_task <- function(task) {
  f <- if (isTRUE(task$sampling_only)) no_shift(task$fit) else task$fit
  v <- with_mt_seed(task$seed, virtual_trials(f, task$design, task$cache, task$om, n_rep = task$n_rep,
                                              n_draws = task$n_draws, seed = task$seed))
  if (identical(task$keep, "diff")) return(v$diff)
  if (!is.null(f$draw_index)) {
    v$draw <- f$draw_index[v$draw]
    cl <- attr(v, "pain_end_clip"); cl$draw <- f$draw_index[cl$draw]; attr(v, "pain_end_clip") <- cl
  }
  v
}

# Net treatment effect in percentage points of pain reduction, placebo-arm relative pain minus LDN-arm relative pain
# (x 100), at the given day(s), for one parameter vector. The trial placebo ceiling is inv_logit(mu_p). Several
# days come from ONE simulation to max(days) (the recursion is the same whatever the horizon).
net_effect_par <- function(par, hyp, fixed, cache, days, schedule = "0:4.5") {
  th <- theta_from(par, hyp, fixed); pm <- inv_logit(par[["mu_p"]]); nd <- max(days)
  a <- simulate_arm(hyp, th, parse_schedule(schedule), nd, cache, pm)$rel
  b <- simulate_arm(hyp, th, parse_schedule("0:0"), nd, cache, pm)$rel
  100 * (b[days + 1] - a[days + 1])
}

# ... at the posterior median.
net_effect_at <- function(fit, cache, days, schedule = "0:4.5")
  net_effect_par(apply(fit$draws, 2, stats::median), fit$hyp, fit$fixed, cache, days, schedule)

# Magnitude guard (review I2). A hypothesis is eligible for Sobol and for power only if its net effect at the
# posterior-median parameters (not the posterior median of the net effect) at `day` EXCEEDS threshold_pp (a benefit, in percentage points of pain reduction). The previous guard
# fell back only when the best hypothesis was NAMED H0, and so decomposed H2's 0.001% effect. A hypothesis with
# no drug term (H0) is never eligible; in the power analysis it is the type-I-error check.
effect_eligibility <- function(fits, cache, day = 84, threshold_pp = 0.5) {
  do.call(rbind, lapply(names(fits), function(h) {
    f <- fits[[h]]; e <- net_effect_at(f, cache, day); null <- identical(f$hyp, "H0")
    elig <- !null && is.finite(e) && e > threshold_pp
    data.frame(hypothesis = h, day = day, net_effect_pp = e, threshold_pp = threshold_pp, eligible = elig,
               sobol_role = if (null) "not run: no drug term" else if (elig) "run" else
                 sprintf("not run: net effect %.3g pp is below the %.2g pp threshold", e, threshold_pp),
               # an eligible hypothesis's rejection rate averages over posterior draws (some with near-zero effect),
               # so it is assurance (posterior-averaged power), which is why it plateaus below 1 (R3a fix round 1)
               power_role = if (null) "type I error check (no drug term)" else if (elig)
                 "assurance (posterior-averaged power)" else
                 sprintf("not run: net effect %.3g pp is below the %.2g pp threshold", e, threshold_pp))
  }))
}

# Model evaluations for one chunk of Sobol rows (no random numbers drawn): a rows x days matrix of net effects.
sobol_eval_task <- function(task) {
  v <- vapply(seq_len(nrow(task$X)), function(i) {
    par <- task$med; par[task$free] <- as.numeric(task$X[i, ])
    net_effect_par(par, task$hyp, task$fixed, task$cache, task$days)
  }, numeric(length(task$days)))
  matrix(v, ncol = length(task$days), byrow = TRUE)
}

# Jansen Sobol indices of the net effect at each day in `days` for one hypothesis. Varied: the mechanism parameters
# (everything except mu_p, tau, sigma_rel and the trial etas), each uniform over its marginal posterior 95% interval
# (review M4: posterior correlation is ignored), all other parameters at posterior medians. The sample matrices come
# from with_mt_seed(seed); the model evaluations go through `map` in n_chunks chunks; the bootstrap CIs of day j use
# with_mt_seed(seed + j). One evaluation gives every day, so both days use the same design points.
sobol_hypothesis <- function(fit, cache, days = c(14, 84), N = 20000, nboot = 200, seed = 3, map = lapply,
                             n_chunks = 40, hypothesis = fit$hyp) {
  free <- setdiff(colnames(fit$draws), c("mu_p", "log_tau", "log_sigma_rel", grep("^eta_", colnames(fit$draws), value = TRUE)))
  # soboljansen needs d >= 2; every drug hypothesis varies kout and kpbo plus Emax, so d >= 3
  if (length(free) < 2) stop("Sobol needs at least 2 varied parameters, got: ", paste(free, collapse = ", "))
  q <- apply(fit$draws[, free, drop = FALSE], 2, stats::quantile, c(0.025, 0.975))
  med <- apply(fit$draws, 2, stats::median)
  X12 <- with_mt_seed(seed, {
    box <- function() { x <- as.data.frame(lapply(free, function(p) stats::runif(N, q[1, p], q[2, p]))); names(x) <- free; x }
    list(box(), box())
  })
  # ST from Jansen (1999); S1 from Sobol (2007) / Saltelli (2002), on the SAME design and evaluations. Jansen's
  # first-order estimator is HIGH-VARIANCE here: it takes V = var(y_A) from one sample, so for an inert parameter it
  # reduces to a common, random-signed offset set by the output's kurtosis (about 9-10 at H1 day 14) and the sampling
  # ratio var(y_B)/var(y_A) (on the H1 posterior at N = 2000, seed 3, an inert log_kpbo got S1 = 0.17 against
  # ST = 2e-5; with seed 17 the offset was 0.007). It is not a systematic upward bias. The Sobol 2007 S1 stays near 0
  # for inert parameters (it can exceed ST slightly there). Jansen's S1 is kept as S1_jansen. V cannot be pooled as
  # var(c(y_A, y_B)) without reimplementing the package's estimators (sensitivity 1.30 hard-codes var(y_A)), so V is
  # left as it is. N = 20000 (R3a fix round 1): at N = 2000 a second seed moved half the non-inert indices outside the
  # first seed's bootstrap CIs.
  sa <- sensitivity::soboljansen(model = NULL, X1 = X12[[1]], X2 = X12[[2]], nboot = nboot)
  s7 <- sensitivity::sobol2007(model = NULL, X1 = X12[[1]], X2 = X12[[2]], nboot = nboot)
  if (!identical(sa$X, s7$X)) stop("sobol_hypothesis: the Jansen and Sobol 2007 designs differ; cannot share evaluations")
  X <- as.matrix(sa$X)
  grp <- ceiling(seq_len(nrow(X)) / ceiling(nrow(X) / n_chunks))
  tasks <- lapply(split(seq_len(nrow(X)), grp), function(r)
    list(X = X[r, , drop = FALSE], med = med, free = free, hyp = fit$hyp, fixed = fit$fixed, cache = cache, days = days))
  Y <- do.call(rbind, map(tasks, sobol_eval_task))
  med_eff <- net_effect_par(med, fit$hyp, fit$fixed, cache, days)
  do.call(rbind, lapply(seq_along(days), function(j) {
    sj <- sa; s1 <- s7
    sj <- with_mt_seed(seed + j, sensitivity::tell(sj, Y[, j]))
    s1 <- with_mt_seed(seed + 100 + j, sensitivity::tell(s1, Y[, j]))
    data.frame(hypothesis = hypothesis, day = days[j], parameter = free,
               S1 = s1$S$original, S1_lo = s1$S$`min. c.i.`, S1_hi = s1$S$`max. c.i.`,
               ST = sj$T$original, ST_lo = sj$T$`min. c.i.`, ST_hi = sj$T$`max. c.i.`,
               S1_jansen = sj$S$original, seed = seed,
               S1_estimator = "sensitivity::sobol2007 (Sobol 2007 / Saltelli 2002 first-order)",
               ST_estimator = "sensitivity::soboljansen (Jansen 1999 total-order)",
               net_effect_median_pp = med_eff[j], output_sd_pp = stats::sd(Y[, j]),
               N = N, nboot = nboot, n_evaluations = nrow(X), row.names = NULL)
  }))
}

POWER_DEFINITION <- paste(
  "rej_rate_with_model_unc: share of replicate trials (n_rep replicates x n_draws posterior draws) whose simulated",
  "difference in NRS change (LDN minus placebo) is below the empirical 5% quantile of the H0 (no-drug) difference",
  "distribution simulated at the same n per arm with an independent seed, with model-structure uncertainty",
  "(TRIAL-level sigma_rel as fitted, shared by both arms of a replicate per review I7, so it largely cancels in the",
  "difference) in both distributions | rej_rate_sampling_only: the same with sigma_rel = 0 in both distributions |",
  "one-sided, nominal 0.05. role says whether a row is genuine power (a hypothesis whose day-84 net effect exceeds",
  "the threshold) or the type-I-error check (H0 against an independent H0 critical value).")

# Rejection-rate curves (review I2). fits must contain H0; hyps are the ELIGIBLE drug hypotheses (genuine power);
# H0 is always added as the type-I-error check. Seeds are as in the original script 07: critical value from H0
# with seed 1001, each candidate with seed 99, the H0 null check with seed 2002.
power_curves <- function(fits, oms, designs, cache, hyps, ns = c(50, 100, 150, 200, 300, 500), n_rep = 200,
                         n_draws = 50, map = lapply, seeds = c(crit = 1001, cand = 99, null = 2002), roles = NULL) {
  if (is.null(fits$H0) || is.null(oms$H0)) stop("power_curves needs the H0 posterior and omegas (critical values and type-I check)")
  miss <- setdiff(hyps, names(fits)); if (length(miss)) stop("power_curves: no posterior for ", paste(miss, collapse = ", "))
  hyps <- setdiff(hyps, "H0")
  thin <- lapply(fits[c("H0", hyps)], thin_fit, n_draws = n_draws)
  tasks <- list(); keys <- list()
  add <- function(key, h, d2, so, seed) {
    tasks[[length(tasks) + 1L]] <<- list(fit = thin[[h]], design = d2, cache = cache, om = oms[[h]], n_rep = n_rep,
                                         n_draws = n_draws, seed = seed, sampling_only = so, keep = "diff")
    keys[[length(keys) + 1L]] <<- key
  }
  for (dz in designs) for (n in ns) for (variant in c("unc", "samp")) {
    d2 <- utils::modifyList(dz, list(n_ldn = n, n_pbo = n)); so <- variant == "samp"
    base <- paste(dz$trial, n, variant, sep = "|")
    add(paste(base, "crit"), "H0", d2, so, seeds[["crit"]])
    add(paste(base, "null"), "H0", d2, so, seeds[["null"]])
    for (h in hyps) add(paste(base, h), h, d2, so, seeds[["cand"]])
  }
  # largest trials first, so a load-balanced cluster does not end on a queue of n = 500 jobs; results are keyed,
  # so the scheduling order cannot change them
  ord <- order(-vapply(tasks, function(t) t$design$n_ldn + t$design$n_pbo, 0))
  res <- vector("list", length(tasks)); res[ord] <- map(tasks[ord], run_vt_task)
  names(res) <- unlist(keys)
  role_of <- function(h) if (!is.null(roles) && !is.null(roles[[h]])) roles[[h]] else
    if (h == "H0") "type I error check (no drug term; H0 against an independent H0 critical value)" else
      "assurance (posterior-averaged power)"
  out <- list()
  for (dz in designs) for (n in ns) for (h in c("H0", hyps)) {
    r <- lapply(c(unc = "unc", samp = "samp"), function(variant) {
      base <- paste(dz$trial, n, variant, sep = "|")
      crit <- unname(stats::quantile(res[[paste(base, "crit")]], 0.05))
      v <- res[[paste(base, if (h == "H0") "null" else h)]]
      list(rate = mean(v < crit), crit = crit, sd = stats::sd(v), mean = mean(v), n_reps = length(v))
    })
    out[[length(out) + 1L]] <- data.frame(hypothesis = h, role = role_of(h), trial = dz$trial, n_per_arm = n,
      rej_rate_with_model_unc = r$unc$rate, rej_rate_sampling_only = r$samp$rate,
      crit_diff_with_model_unc = r$unc$crit, crit_diff_sampling_only = r$samp$crit,
      mean_diff_with_model_unc = r$unc$mean, sd_diff_with_model_unc = r$unc$sd, sd_diff_sampling_only = r$samp$sd,
      n_reps = r$unc$n_reps, n_rep = n_rep, n_draws = n_draws, definition = POWER_DEFINITION)
  }
  out <- do.call(rbind, out)
  # A hypothesis named in `roles` but not run (below the effect threshold) gets one explanatory row, so the file
  # reads on its own (R3a fix round 1).
  for (h in setdiff(names(roles), c("H0", hyps))) {
    na <- out[1, , drop = FALSE]; na[1, ] <- NA
    na$hypothesis <- h; na$n_rep <- n_rep; na$n_draws <- n_draws; na$definition <- POWER_DEFINITION
    na$role <- paste0(roles[[h]], "; no rejection rates computed (see effect_eligibility.csv)")
    out <- rbind(out, na)
  }
  rownames(out) <- NULL
  out
}

# Seed agreement of two Sobol runs that differ only in their seed (R3a fix round 1). For every hypothesis, day,
# parameter and index (S1, ST): both estimates, both bootstrap 95% CI half-widths, and a flag where the two
# estimates differ by more than the larger half-width.
sobol_seed_agreement <- function(a, b, seeds) {
  a <- a[!is.na(a$parameter), ]; b <- b[!is.na(b$parameter), ]
  ka <- paste(a$hypothesis, a$day, a$parameter); kb <- paste(b$hypothesis, b$day, b$parameter)
  if (length(ka) != length(kb) || !setequal(ka, kb))
    stop("sobol_seed_agreement: the two runs must cover the same hypothesis/day/parameter rows")
  b <- b[match(ka, kb), ]
  do.call(rbind, lapply(c("S1", "ST"), function(ix) {
    hwa <- (a[[paste0(ix, "_hi")]] - a[[paste0(ix, "_lo")]]) / 2; hwb <- (b[[paste0(ix, "_hi")]] - b[[paste0(ix, "_lo")]]) / 2
    d <- abs(a[[ix]] - b[[ix]]); mx <- pmax(hwa, hwb)
    g <- data.frame(hypothesis = a$hypothesis, day = a$day, parameter = a$parameter, index = ix,
                    ea = a[[ix]], eb = b[[ix]], ha = hwa, hb = hwb, abs_diff = d, max_ci_half_width = mx, disagree = d > mx)
    names(g)[5:8] <- c(paste0("estimate_seed", seeds[1]), paste0("estimate_seed", seeds[2]),
                       paste0("ci_half_width_seed", seeds[1]), paste0("ci_half_width_seed", seeds[2]))
    rownames(g) <- NULL
    g
  }))
}

# The row(s) written for a hypothesis whose effect is below the threshold: same columns as sobol_hypothesis(),
# indices NA, the net effect at each day recorded (review I2: say so rather than skip silently).
sobol_skipped <- function(hypothesis, days, net_effect_pp)
  data.frame(hypothesis = hypothesis, day = days, parameter = NA_character_, S1 = NA_real_, S1_lo = NA_real_,
             S1_hi = NA_real_, ST = NA_real_, ST_lo = NA_real_, ST_hi = NA_real_, S1_jansen = NA_real_, seed = NA_real_,
             S1_estimator = NA_character_, ST_estimator = NA_character_, net_effect_median_pp = net_effect_pp,
             output_sd_pp = NA_real_, N = NA_integer_, nboot = NA_integer_, n_evaluations = NA_integer_)
