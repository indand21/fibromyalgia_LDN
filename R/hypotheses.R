HYP_PARS <- list(H0 = character(0), H1 = c("logit_Emax", "log_kA", "log_gamma"),
                 H2 = "logit_Emax", H3 = c("logit_Emax", "log10_Kd_OGFr"))

init_par <- function(hyp, trials) {
  p <- c(mu_p = logit(0.12), log_tau = log(0.3), log_kpbo = log(0.1), log_kout = log(0.1),
         log_sigma_rel = log(0.03),
         stats::setNames(rep(0, length(trials)), paste0("eta_", trials)))
  extra <- c(logit_Emax = logit(0.2), log_kA = log(0.1), log_gamma = 0, log10_Kd_OGFr = 0)
  c(p, extra[HYP_PARS[[hyp]]])
}

theta_from <- function(par, hyp, fixed) {
  th <- list(kpbo = exp(par[["log_kpbo"]]), kout = exp(par[["log_kout"]]),
             IC50_TLR4 = fixed$IC50_TLR4, Kd_OGFr = fixed$Kd_OGFr)
  if (hyp != "H0") th$Emax <- inv_logit(par[["logit_Emax"]])
  if (hyp == "H1") { th$kA <- exp(par[["log_kA"]]); th$gamma <- exp(par[["log_gamma"]]) }
  if (hyp == "H3" && (is.null(fixed$Kd_OGFr) || is.na(fixed$Kd_OGFr))) th$Kd_OGFr <- 10^par[["log10_Kd_OGFr"]]
  th
}

pmax_trial <- function(par, trial) inv_logit(par[["mu_p"]] + par[[paste0("eta_", trial)]])

# Scale conversion for a prior that a curated file records on the NATURAL scale but the sampler needs on the LOG
# scale (review C2). Byon 2010 Table III reports the placebo onset rate as kpbo = 0.106 /day with a standard
# ERROR of 0.00833 /day; log_prior() needs the SD of log(kpbo), which by the delta method is sd/value =
# 0.00833/0.106 = 0.0786, about 9.4x wider than the per-day number. The CSV stays the source of record (a per-day
# value with its per-day SE); the conversion happens here, at use. The `scale` column decides and is mandatory:
# a `per_day` row is always converted, a row already on the log scale is used as it stands.
PRIOR_LOG_SCALES <- c("log", "log_per_day")
PRIOR_NATURAL_SCALES <- "per_day"
log_scale_prior <- function(priors, nm) {
  i <- which(priors$parameter == nm)
  if (length(i) != 1L) stop("priors need exactly one '", nm, "' row (found ", length(i), ")")
  v <- priors$value[i]; s <- priors$sd[i]
  if (!is.finite(v) || !is.finite(s)) stop(nm, ": prior needs a finite value and a finite sd")
  if (!"scale" %in% names(priors))
    stop(nm, ": priors have no `scale` column, so a natural-scale standard error cannot be told apart from a log-scale SD")
  sc <- priors$scale[i]
  sc <- if (is.na(sc)) "" else trimws(as.character(sc))
  if (!nzchar(sc))
    stop(nm, ": prior row has a blank `scale`; cannot tell a natural-scale standard error from a log-scale SD")
  if (sc %in% PRIOR_LOG_SCALES) return(list(mean = v, sd = s, converted = FALSE, scale = sc))
  if (!sc %in% PRIOR_NATURAL_SCALES)
    stop(nm, ": unsupported prior scale '", sc, "' (expected one of ",
         paste(c(PRIOR_NATURAL_SCALES, PRIOR_LOG_SCALES), collapse = ", "), ")")
  if (v <= 0) stop(nm, ": a ", sc, " prior value must be positive to convert to the log scale (found ", v, ")")
  list(mean = log(v), sd = s / v, converted = TRUE, scale = sc)
}

log_prior <- function(par, hyp, priors, trials) {
  need <- c("pbo_rel_mean", "pbo_logit_sd", "kpbo_per_day")
  miss <- setdiff(need, priors$parameter)
  if (length(miss)) stop("priors missing required parameter(s): ", paste(miss, collapse = ", "))
  pv <- function(nm, col = "value") priors[[col]][priors$parameter == nm]
  kp <- log_scale_prior(priors, "kpbo_per_day")
  tau <- exp(par[["log_tau"]])
  # mu_p carries the meta-analytic pooled placebo mean with its SE (a prior LOCATION, legitimately tight);
  # between-trial heterogeneity is carried by tau, and a trial's ceiling is inv_logit(mu_p + eta_trial)
  # with eta ~ N(0, tau) and tau ~ half-normal(0.5) (review I9; controller ruling 2026-10-01).
  lp <- stats::dnorm(par[["mu_p"]], logit(pv("pbo_rel_mean")), pv("pbo_logit_sd"), log = TRUE) +
    stats::dnorm(tau, 0, 0.5, log = TRUE) + par[["log_tau"]] +                 # half-normal(0.5) + Jacobian
    sum(stats::dnorm(par[paste0("eta_", trials)], 0, tau, log = TRUE)) +
    stats::dnorm(par[["log_kpbo"]], kp$mean, kp$sd, log = TRUE) +              # per-day SE converted to a log-scale SD
    stats::dnorm(par[["log_kout"]], log(0.1), 1, log = TRUE) +
    stats::dnorm(exp(par[["log_sigma_rel"]]), 0, 0.05, log = TRUE) + par[["log_sigma_rel"]]  # half-normal(0.05) + Jacobian
  if (hyp != "H0") { z <- par[["logit_Emax"]]                                     # Emax ~ U(0,1), stable Jacobian
                     lp <- lp + stats::plogis(z, log.p = TRUE) + stats::plogis(-z, log.p = TRUE) }
  if (hyp == "H1") lp <- lp + stats::dnorm(par[["log_kA"]], log(0.1), 1, log = TRUE) +
    stats::dnorm(par[["log_gamma"]], 0, 0.5, log = TRUE)
  if ("log10_Kd_OGFr" %in% names(par)) lp <- lp + stats::dnorm(par[["log10_Kd_OGFr"]], 0, 2, log = TRUE)
  lp
}

# Scale of an outcome on its own units: pct_reduction is in percent, change_nrs in NRS points.
outcome_scale <- function(outcome, baseline_nrs) ifelse(outcome == "pct_reduction", 100, baseline_nrs)

# SIGNED scale of the model-structure residual on an outcome's own units. The residual is a perturbation of the
# relative pain `rel`, and d(pct_reduction)/d(rel) = -100 while d(change_nrs)/d(rel) = baseline_nrs. The magnitude
# is outcome_scale(); the sign matters only now that one shared term enters several arms at once (review I7).
outcome_shift <- function(outcome, baseline_nrs) ifelse(outcome == "pct_reduction", -100, baseline_nrs)

# ---- Shared trial-level model-structure residual (review I7), integrated analytically ----
# For the arms i = 1..k of ONE trial:
#   y_i = m_i + z * u_i + e_i,   z ~ N(0, 1) shared by the trial,   e_i ~ N(0, se_i^2),
#   u_i = sigma_rel * outcome_shift_i.
# Hence r = y - m ~ N(0, Sigma) with Sigma = D + u u', D = diag(se^2): a rank-one (compound-symmetric when the
# u_i are equal) covariance across the arms of a trial, independent between trials. With q = sum_i u_i^2/se_i^2,
# the matrix determinant lemma gives   log|Sigma| = sum_i log se_i^2 + log(1 + q),
# and the Sherman-Morrison identity gives
#   r' Sigma^{-1} r = sum_i r_i^2/se_i^2 - (sum_i r_i u_i / se_i^2)^2 / (1 + q).
# So z is integrated out in closed form; no quadrature or data augmentation is needed. A one-arm trial reduces to
# dnorm(r, 0, sqrt(se^2 + u^2)), which is exactly the old per-arm formula, as does sigma_rel = 0.
shared_resid_loglik <- function(r, se, u) {
  k <- length(r); iv <- 1 / se^2; q <- sum(u^2 * iv)
  -0.5 * (k * log(2 * pi) + sum(log(se^2)) + log1p(q) + sum(r^2 * iv) - sum(r * u * iv)^2 / (1 + q))
}

# Pointwise terms for PSIS-LOO (R2 fix round 1, IMPORTANT 2). With a shared term the likelihood does not
# factorise over the arms of a trial, and leave-one-out then needs the CONDITIONAL p(y_i | y_-i, theta)
# (Burkner, Gabry & Vehtari 2021, Comput Stat 36:1243, "Efficient leave-one-out cross-validation for Bayesian
# non-factorized normal and Student-t models"). With P = Sigma^-1, g = P r and cbar_i = P_ii:
#   log p(y_i | y_-i, theta) = -1/2 log(2 pi) + 1/2 log(cbar_i) - 1/2 g_i^2 / cbar_i .
# P comes in closed form from Sherman-Morrison: P = D^-1 - (D^-1 u)(D^-1 u)' / (1 + q), so
#   g_i = r_i/se_i^2 - (u_i/se_i^2) * sum_j(r_j u_j/se_j^2) / (1 + q),   cbar_i = 1/se_i^2 - (u_i/se_i^2)^2 / (1 + q).
# These terms do NOT sum to the joint (only the joint, shared_resid_loglik(), enters the posterior), and unlike
# the chain-rule terms used in R2 they do not depend on the order of the rows. A one-arm trial reduces to
# dnorm(r, 0, sqrt(se^2 + u^2)). (Arm-level LOO still violates independence through eta_trial; review I5.)
shared_resid_loglik_pointwise <- function(r, se, u) {
  iv <- 1 / se^2; q <- sum(u^2 * iv); ui <- u * iv
  g <- r * iv - ui * sum(r * ui) / (1 + q)
  cbar <- iv - ui^2 / (1 + q)
  -0.5 * log(2 * pi) + 0.5 * log(cbar) - 0.5 * g^2 / cbar
}

log_lik_arms <- function(par, hyp, arms, cache, fixed, pointwise = FALSE) {
  if (any(!is.finite(arms$se)))
    stop("log_lik_arms: non-finite se in arm(s): ", paste(paste(arms$trial, arms$arm)[!is.finite(arms$se)], collapse = "; "),
         "; arms without reported dispersion must be excluded before the likelihood")
  th <- theta_from(par, hyp, fixed); sig <- exp(par[["log_sigma_rel"]])
  fit <- vapply(seq_len(nrow(arms)), function(i) {
    a <- arms[i, ]
    sim <- simulate_arm(hyp, th, parse_schedule(a$schedule), a$day, cache, pmax_trial(par, a$trial))
    arm_outcome(sim, a$day, a$outcome, a$baseline_nrs)
  }, 0)
  r <- arms$mean - fit
  u <- sig * outcome_shift(arms$outcome, arms$baseline_nrs)
  grp <- split(seq_len(nrow(arms)), arms$trial)                  # one shared residual per TRIAL
  if (!pointwise) return(sum(vapply(grp, function(ix) shared_resid_loglik(r[ix], arms$se[ix], u[ix]), 0)))
  ll <- numeric(nrow(arms))
  for (ix in grp) ll[ix] <- shared_resid_loglik_pointwise(r[ix], arms$se[ix], u[ix])
  ll
}

# PSIS-LOO on the training arms for one fitted hypothesis (script 04; R2 items 7 and fix round 1). The pointwise
# terms are the leave-one-out conditionals of log_lik_arms(pointwise = TRUE). p_loo is reported against the number
# of free parameters; p_loo is flagged unreliable when any Pareto k exceeds PARETO_K_RELIABLE (Vehtari et al. 2024).
PARETO_K_RELIABLE <- 0.7
training_loo <- function(f, hyp, arms, cache, fixed, n_keep = 2000) {
  keep <- round(seq(1, nrow(f$draws), length.out = n_keep))
  ll <- do.call(rbind, lapply(keep, function(i) log_lik_arms(f$draws[i, ], hyp, arms, cache, fixed, pointwise = TRUE)))
  reff <- loo::relative_eff(exp(ll), chain_id = f$chain[keep])
  lo <- loo::loo(ll, r_eff = reff)
  k <- lo$diagnostics$pareto_k
  p_loo <- lo$estimates["p_loo", 1]; n_par <- ncol(f$draws)
  data.frame(hypothesis = hyp, elpd_loo = lo$estimates["elpd_loo", 1], se = lo$estimates["elpd_loo", 2],
             p_loo = p_loo, n_params = n_par, prior_determined_fraction = 1 - p_loo / n_par,
             max_pareto_k = max(k), n_k_gt_0.7 = sum(k > PARETO_K_RELIABLE),
             ranking_valid = sum(k > PARETO_K_RELIABLE) == 0, p_loo_reliable = max(k) <= PARETO_K_RELIABLE)
}

# The citable p_loo table (output/results/p_loo_identifiability.csv) from training_loo() rows. Each row's `note`
# states whether its prior-determined fraction can be read as a number or is indicative only; the reliability
# flag rests on the Pareto k diagnostic of PSIS-LOO (R2 fix round 2, ruling 3).
p_loo_table <- function(loo_tab, n_obs) {
  note <- ifelse(loo_tab$p_loo_reliable,
    sprintf(paste("p_loo reliable: max Pareto k = %.2f <= %.1f. prior_determined_fraction = 1 - p_loo / n_free_parameters.",
                  "The reliability flag rests on the Pareto k diagnostic of PSIS-LOO."), loo_tab$max_pareto_k, PARETO_K_RELIABLE),
    sprintf(paste("p_loo UNRELIABLE: max Pareto k = %.2f > %.1f, so p_loo and the prior_determined_fraction are",
                  "indicative only. The reliability flag rests on the Pareto k diagnostic of PSIS-LOO."), loo_tab$max_pareto_k, PARETO_K_RELIABLE))
  data.frame(hypothesis = loo_tab$hypothesis, n_free_parameters = loo_tab$n_params,
    p_loo = loo_tab$p_loo, prior_determined_fraction = loo_tab$prior_determined_fraction,
    elpd_loo = loo_tab$elpd_loo, elpd_loo_se = loo_tab$se,
    max_pareto_k = loo_tab$max_pareto_k, p_loo_reliable = loo_tab$p_loo_reliable, ranking_valid = loo_tab$ranking_valid,
    n_training_observations = n_obs, note = note,
    method_note = paste("p_loo from the leave-one-out conditionals p(y_i | y_-i, theta) of the non-factorising likelihood",
                        "(Burkner, Gabry & Vehtari 2021). Both training arms come from ONE trial and share its eta and",
                        "model-structure residual, so the training elpd must not be used for ranking."))
}

# Training arms used in the fit, with every exclusion reason (R3a, review M2). An arm without a reported dispersion
# cannot enter a likelihood that needs an SE (never imputed). Where the arm's own `location` records that its
# endpoint is not a pain-specific score (Younger 2009: overall fibromyalgia symptom severity on a 0-100 VAS),
# that second and stronger reason is stated as well, quoting the data.
training_arms_table <- function(arms_all) {
  use <- is.finite(arms_all$se)
  loc <- if ("location" %in% names(arms_all)) arms_all$location else rep("", nrow(arms_all))
  not_pain <- grepl("not a pain-specific score", loc, ignore.case = TRUE)
  m <- regexpr("overall [^;]*not a pain-specific score", loc, ignore.case = TRUE)
  endpoint_txt <- rep("", length(loc)); endpoint_txt[m > 0] <- regmatches(loc, m)
  endpoint_txt[not_pain & endpoint_txt == ""] <- "not a pain-specific score"
  reason <- ifelse(use, "", "no dispersion reported in source")
  reason <- ifelse(!use & not_pain,
                   paste0(reason, "; endpoint mismatch: ", endpoint_txt, " (efficacy_training.csv location column)"), reason)
  data.frame(trial = arms_all$trial, arm = arms_all$arm, day = arms_all$day, mean = arms_all$mean, se = arms_all$se,
             status = ifelse(use, "included", "excluded"), reason = reason)
}

# Dose at which an effect first reaches `frac` of its maximum over the grid, by linear interpolation on the
# rising limb (doses ascending, up to the peak). NA when the effect never rises above zero.
dose_for_fraction <- function(doses, eff, frac) {
  stopifnot(length(doses) == length(eff), !is.unsorted(doses, strictly = TRUE))
  mx <- max(eff); if (!(mx > 0)) return(NA_real_)
  thr <- frac * mx; pk <- which.max(eff)
  if (eff[1] >= thr) return(doses[1])
  for (i in seq_len(pk)[-1]) if (eff[i] >= thr && eff[i - 1] < thr)
    return(doses[i - 1] + (thr - eff[i - 1]) / (eff[i] - eff[i - 1]) * (doses[i] - doses[i - 1]))
  NA_real_
}

# Peak and model EDx of an effect-versus-dose grid (R3a fix round 1). EDx is defined relative to the GRID maximum,
# so it is only meaningful when the peak is interior: if the effect is still rising at the top dose (H3, and H2's
# linear curve) the "maximum" is the grid top and EDx moves with it (H3's ED95 was 5.53 mg at a 6 mg top and 44.2 mg
# at a 50 mg top, a spurious match with the observed 5.40 mg). EDx is then NA with the reason in ed_note; likewise
# when the peak sits at the lowest non-zero dose (the true peak may lie below the grid) or there is no effect.
bruun_ed <- function(doses, eff) {
  nz <- doses[doses > 0]; pk <- doses[which.max(eff)]
  pos <- if (!(max(eff) > 0)) "no effect" else if (pk == min(nz)) "lowest non-zero dose on grid (true peak may lie below)" else
    if (pk == max(doses)) "highest dose on grid" else "interior"
  note <- switch(pos, "interior" = "",
                 "highest dose on grid" = "not estimable: effect still rising at grid top (EDx would move with the grid top)",
                 "no effect" = "not estimable: no effect on the grid",
                 "not estimable: peak at the lowest non-zero dose (true peak may lie below the grid)")
  ok <- pos == "interior"
  list(peak_dose = pk, peak_position = pos, ED50 = if (ok) dose_for_fraction(doses, eff, 0.5) else NA_real_,
       ED95 = if (ok) dose_for_fraction(doses, eff, 0.95) else NA_real_, ed_note = note)
}
BRUUN_COMPARISON_NOTE <- paste(
  "descriptive comparison only: the observed Bruun-Plesner ED50/ED95 are quantal estimands (the dose at which 50%/95%",
  "of patients respond, from an up-and-down design), while the model EDx is a graded fraction (50%/95%) of the maximum",
  "mean effect on the dose grid")

# Bruun-Plesner consistency check (R3a, review M7): effect against dose at `day` (posterior medians, no placebo),
# the dose of peak effect, where on the grid the peak sits, and the model ED50 and ED95 (dose reaching 50% / 95%
# of the peak effect on the rising limb) next to the observed up-and-down ED50 and ED95. Every grid dose must have
# an exposure-cache profile; scripts/02 builds the cache with doses down to 0.1 mg for this.
bruun_check <- function(fits, cache, doses, obs_ED50, obs_ED95, day = 365) {
  miss <- setdiff(as.character(doses), names(cache))
  if (length(miss)) stop("bruun_check: exposure cache has no profile for dose(s) ", paste(miss, collapse = ", "),
                         " mg; add them to the dose grid in scripts/02_binding_bottomup.R")
  doses <- sort(doses); nz <- doses[doses > 0]
  grid <- list(); summ <- list()
  for (h in names(fits)) {
    f <- fits[[h]]; med <- apply(f$draws, 2, stats::median); th <- theta_from(med, h, f$fixed)
    eff <- vapply(doses, function(d) 1 - utils::tail(simulate_arm(h, th, parse_schedule(paste0("0:", d)), day, cache, 0)$rel, 1), 0)
    grid[[h]] <- data.frame(hypothesis = h, dose_mg = doses, effect = eff)
    ed <- bruun_ed(doses, eff)
    summ[[h]] <- data.frame(hypothesis = h, peak_dose_mg = ed$peak_dose, peak_effect = max(eff), peak_position = ed$peak_position,
                            model_ED50_mg = ed$ED50, model_ED95_mg = ed$ED95, ed_note = ed$ed_note,
                            observed_ED50_mg = obs_ED50, observed_ED95_mg = obs_ED95, lowest_nonzero_dose_mg = min(nz),
                            day = day, comparison_note = BRUUN_COMPARISON_NOTE,
                            definition = paste("effect = 1 - relative pain at day", day, "without placebo (posterior medians);",
                                               "model EDx = dose at which the effect first reaches x% of its grid maximum (linear",
                                               "interpolation on the rising limb), reported only when the peak is interior;",
                                               "observed values are Bruun-Plesner 2020 up-and-down estimates"))
  }
  list(grid = do.call(rbind, grid), summary = do.call(rbind, summ))
}
