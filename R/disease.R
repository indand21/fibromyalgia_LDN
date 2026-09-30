dose_by_day <- function(schedule, n_days) {
  if (n_days < 1) stop("n_days must be >= 1")
  schedule$dose_mg[findInterval(seq_len(n_days) - 1, schedule$start_day)]
}

get_prof <- function(cache, dose) {
  key <- as.character(dose)
  if (is.null(cache[[key]])) stop("exposure cache has no profile for dose ", key, " mg")
  cache[[key]]
}

# Only for parameters a hypothesis does not use: NULL or NA means "parameter not used by this hypothesis": Inf makes its term 0.
theta_or <- function(x, default) if (is.null(x) || is.na(x)) default else x

# Returns daily drug inhibition terms (vector over subjects) and updated adaptation state A.
daily_inh <- function(hyp, theta, m, A, Emax) {
  switch(hyp,
    H0 = list(inh = 0 * Emax, A = A),
    H1 = { A <- A + (1 - exp(-theta$kA)) * (m$mRO - A); list(inh = Emax * A * m$free_term, A = A) },
    H2 = list(inh = Emax * m$tlr4, A = A),
    H3 = list(inh = Emax * m$ogfr, A = A),
    stop("unknown hypothesis ", hyp))
}

simulate_arm_vec <- function(hyp, theta, schedule, n_days, cache, Pmax_vec, Emax_vec) {
  req <- list(H0 = character(), H1 = c("kA", "gamma"), H2 = "IC50_TLR4", H3 = "Kd_OGFr")[[hyp]]
  if (is.null(req)) stop("unknown hypothesis ", hyp)
  need <- c("kout", "kpbo", req)
  ok <- vapply(need, function(k) { v <- theta[[k]]; is.numeric(v) && length(v) == 1 && is.finite(v) }, logical(1))
  if (!all(ok)) stop(hyp, " requires finite theta: ", paste(need[!ok], collapse = ", "))
  n <- length(Pmax_vec); stopifnot(length(Emax_vec) == n)
  stopifnot(all(Pmax_vec >= 0 & Pmax_vec < 1))
  if (hyp != "H0" && (anyNA(Emax_vec) || any(Emax_vec < 0 | Emax_vec > 1)))
    stop(hyp, " requires Emax in [0, 1]")
  doses <- dose_by_day(schedule, n_days)
  means <- lapply(unique(doses), function(d) profile_means(get_prof(cache, d),
    gamma = theta_or(theta$gamma, 1),
    IC50_TLR4 = theta_or(theta$IC50_TLR4, Inf),
    Kd_OGFr = theta_or(theta$Kd_OGFr, Inf)))
  names(means) <- as.character(unique(doses))
  S <- matrix(1, n_days + 1, n); A <- rep(0, n); decay <- exp(-theta$kout)
  for (d in seq_len(n_days)) {
    r <- daily_inh(hyp, theta, means[[as.character(doses[d])]], A, Emax_vec)
    A <- r$A; target <- 1 - pmin(r$inh, 0.999)  # numerical guard only; Emax <= 1 and occupancy terms <= 1 keep inh < 1 in practice
    S[d + 1, ] <- target + (S[d, ] - target) * decay
  }
  t <- 0:n_days
  Pbo <- outer(1 - exp(-theta$kpbo * t), Pmax_vec)
  S * (1 - Pbo)
}

simulate_arm <- function(hyp, theta, schedule, n_days, cache, Pmax) {
  if (is.null(theta$Emax) && hyp != "H0") stop(hyp, " requires theta$Emax")
  Emax <- if (is.null(theta$Emax)) 0 else theta$Emax
  rel <- simulate_arm_vec(hyp, theta, schedule, n_days, cache, Pmax, Emax)[, 1]
  Pbo <- Pmax * (1 - exp(-theta$kpbo * (0:n_days)))
  data.frame(day = 0:n_days, S = rel / (1 - Pbo), Pbo = Pbo, rel = rel)
}

arm_outcome <- function(sim, day, outcome, baseline_nrs) {
  rel <- sim$rel[match(day, sim$day)]
  switch(outcome,
    pct_reduction = 100 * (1 - rel),
    change_nrs = baseline_nrs * (rel - 1),
    stop("unknown outcome ", outcome))
}
