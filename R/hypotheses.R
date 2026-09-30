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

log_prior <- function(par, hyp, priors, trials) {
  need <- c("pbo_rel_mean", "pbo_logit_sd", "kpbo_per_day")
  miss <- setdiff(need, priors$parameter)
  if (length(miss)) stop("priors missing required parameter(s): ", paste(miss, collapse = ", "))
  pv <- function(nm, col = "value") priors[[col]][priors$parameter == nm]
  tau <- exp(par[["log_tau"]])
  lp <- stats::dnorm(par[["mu_p"]], logit(pv("pbo_rel_mean")), pv("pbo_logit_sd"), log = TRUE) +
    stats::dnorm(tau, 0, 0.5, log = TRUE) + par[["log_tau"]] +                 # half-normal(0.5) + Jacobian
    sum(stats::dnorm(par[paste0("eta_", trials)], 0, tau, log = TRUE)) +
    stats::dnorm(par[["log_kpbo"]], log(pv("kpbo_per_day")), pv("kpbo_per_day", "sd"), log = TRUE) +
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

log_lik_arms <- function(par, hyp, arms, cache, fixed, pointwise = FALSE) {
  if (any(!is.finite(arms$se)))
    stop("log_lik_arms: non-finite se in arm(s): ", paste(paste(arms$trial, arms$arm)[!is.finite(arms$se)], collapse = "; "),
         "; arms without reported dispersion must be excluded before the likelihood")
  th <- theta_from(par, hyp, fixed); sig <- exp(par[["log_sigma_rel"]])
  ll <- vapply(seq_len(nrow(arms)), function(i) {
    a <- arms[i, ]
    sim <- simulate_arm(hyp, th, parse_schedule(a$schedule), a$day, cache, pmax_trial(par, a$trial))
    sd_i <- sqrt(a$se^2 + (sig * outcome_scale(a$outcome, a$baseline_nrs))^2)
    stats::dnorm(a$mean, arm_outcome(sim, a$day, a$outcome, a$baseline_nrs), sd_i, log = TRUE)
  }, 0)
  if (pointwise) ll else sum(ll)
}
