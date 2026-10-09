# Reporting helpers shared by scripts 10 and 11 and by the run driver's sensitivity summary (R3b).
# Pure functions of result tables: no file access, no random numbers.

# ---- Prior SD of every sampled parameter, on the scale the sampler uses (Table 2 relabel, review M1) ----
# Closed forms of the priors in log_prior() (R/hypotheses.R); tested against log_prior() by numerical integration.
#   mu_p           N(logit(pbo_rel_mean), pbo_logit_sd)                       -> pbo_logit_sd
#   log_tau        tau ~ half-normal(0.5), so log tau = log 0.5 + log|Z|       -> sqrt(pi^2 / 8) = 1.111
#   log_sigma_rel  sigma_rel ~ half-normal(0.05), same argument                -> sqrt(pi^2 / 8)
#   eta_<trial>    eta | tau ~ N(0, tau), tau ~ half-normal(0.5): marginal     -> sqrt(E[tau^2]) = 0.5
#   log_kpbo       N(log 0.106, per-day SE / value), log_scale_prior()         -> 0.0786
#   log_kout, log_kA  N(log 0.1, 1)                                            -> 1
#   log_gamma      N(0, 0.5)                                                   -> 0.5
#   logit_Emax     Emax ~ U(0, 1), so logit(Emax) is standard logistic         -> pi / sqrt(3) = 1.814
#   log10_Kd_OGFr  N(0, 2)                                                     -> 2
PRIOR_SD_HALFNORMAL_LOG <- sqrt(pi^2 / 8)
prior_sd_on_sampler_scale <- function(variable, priors) {
  one <- function(v) {
    if (v == "mu_p") return(priors$value[priors$parameter == "pbo_logit_sd"])
    if (v %in% c("log_tau", "log_sigma_rel")) return(PRIOR_SD_HALFNORMAL_LOG)
    if (startsWith(v, "eta_")) return(0.5)
    if (v == "log_kpbo") return(log_scale_prior(priors, "kpbo_per_day")$sd)
    if (v %in% c("log_kout", "log_kA")) return(1)
    if (v == "log_gamma") return(0.5)
    if (v == "logit_Emax") return(pi / sqrt(3))
    if (v == "log10_Kd_OGFr") return(2)
    stop("prior_sd_on_sampler_scale: no prior SD defined for parameter '", v,
         "'; add its closed form here (Table 2 labels are computed, never typed by hand)", call. = FALSE)
  }
  vapply(variable, one, 0, USE.NAMES = FALSE)
}

# Table 2 source_type by computation (controller ruling, R3b item 3): posterior SD / prior SD >= threshold means the
# two training observations did not narrow the prior, so the parameter is labelled prior-determined.
PRIOR_DETERMINED_RATIO <- 0.9
posterior_source_type <- function(post_sd, prior_sd, threshold = PRIOR_DETERMINED_RATIO) {
  r <- post_sd / prior_sd
  ifelse(!is.finite(r), "fitted (training; posterior / prior SD ratio not computable)",
         ifelse(r >= threshold, sprintf("prior-determined (posterior SD / prior SD = %.2f)", r),
                sprintf("fitted (training; posterior SD / prior SD = %.2f)", r)))
}

# ---- Pairwise held-out elpd contrasts on the PRIMARY observations ----
# diff = sum of per-observation lpd differences; se_diff = sd(per-observation differences) * sqrt(n) (the
# loo-compare observation-level SE); mcse_diff = Monte Carlo SE of the same difference from paired predictive draws
# (script 06, test_elpd_contrast_mcse.csv). With n = 2, se_diff = |d1 - d2| exactly, so it is unstable in BOTH
# directions: tiny when the two per-trial differences happen to agree, large when they disagree in sign. A ratio
# |diff / se_diff| >= 3 can only arise when the two differences share a sign and agree within a factor of 2, i.e. from
# a collapsed SE; such rows are flagged (controller ruling: the H1-vs-H3 ratio of 11.4 is not strong evidence).
SE_COLLAPSE_RATIO <- 3
elpd_contrasts <- function(tp, cm = NULL) {
  need <- c("hypothesis", "trial", "endpoint", "day", "role", "lpd")
  miss <- setdiff(need, names(tp))
  if (length(miss)) stop("elpd_contrasts: test predictions lack column(s) ", paste(miss, collapse = ", "), call. = FALSE)
  pr <- tp[tp$endpoint == "diff_nrs" & tp$role == "primary", c("hypothesis", "trial", "day", "lpd")]
  if (!nrow(pr)) stop("elpd_contrasts: no primary diff_nrs rows", call. = FALSE)
  key <- unique(pr[c("trial", "day")]); key <- key[order(tolower(key$trial), key$day), ]
  n_obs <- nrow(key)
  lp <- function(h) { x <- pr[pr$hypothesis == h, ]
    v <- x$lpd[match(paste(key$trial, key$day), paste(x$trial, x$day))]
    if (anyNA(v)) stop("elpd_contrasts: hypothesis ", h, " lacks a primary lpd for some trial", call. = FALSE)
    v }
  mc <- function(a, b) { if (is.null(cm)) return(NA_real_)
    v <- cm$mcse_diff[cm$hypothesis_a == a & cm$hypothesis_b == b]; if (length(v) == 1) v else NA_real_ }
  hh <- sort(unique(pr$hypothesis)); out <- list()
  for (i in seq_along(hh)) for (j in seq_along(hh)) if (i < j) {
    d <- lp(hh[i]) - lp(hh[j]); est <- sum(d); se <- stats::sd(d) * sqrt(n_obs); m <- mc(hh[i], hh[j])
    ratio <- est / se
    collapsed <- n_obs <= 2 && is.finite(ratio) && abs(ratio) >= SE_COLLAPSE_RATIO
    out[[length(out) + 1]] <- data.frame(hypothesis_a = hh[i], hypothesis_b = hh[j], elpd_diff_a_minus_b = est,
      se_diff_observation_level = se, n_obs = n_obs, ratio_diff_over_se = ratio,
      mcse_diff_monte_carlo = m, ratio_diff_over_mcse = est / m,
      per_trial_lpd_diff = paste(sprintf("%s d%g %+.3f", key$trial, key$day, d), collapse = "; "),
      se_collapsed_flag = collapsed,
      flag_note = if (collapsed) sprintf("RATIO INFLATED BY A COLLAPSED SE: the %d per-trial differences agree in sign and size, so se_diff = |d1 - d2| is spuriously small; not strong evidence", n_obs) else "",
      se_kind = paste("se_diff_observation_level = SE over the n held-out primary endpoints (loo-compare style);",
                      "mcse_diff_monte_carlo = Monte Carlo (posterior-predictive draw) SE of the same difference, from paired draws"),
      note = sprintf("PRIMARY endpoints only (%s). With n = %d the observation-level SE is unstable in both directions (spuriously small when the per-trial differences agree, large when they disagree); indicative only, not decisive.",
                     paste(sprintf("%s d%g", key$trial, key$day), collapse = ", "), n_obs),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# ---- Fig 3 axis categories (R4 fix round 1, review I4): ONE function labels both the observed arms and the posterior
# predictions, from the same trial, arm and day columns, so the two layers always share their x categories.
# "younger2013", "ldn_4.5mg", 84 -> "Younger 2013\nLDN 4.5 mg (day 84)"; an arm name without a dose gives "LDN".
arm_label <- function(trial, arm, day) {
  t <- trimws(gsub("\\s+", " ", gsub("(\\d{4})", " \\1", tools::toTitleCase(gsub("([a-z])(\\d)", "\\1 \\2", trial)))))
  dose <- vapply(seq_along(arm), function(i) {
    m <- regmatches(arm[i], regexpr("[0-9]+([.][0-9]+)?", arm[i]))
    if (length(m)) paste0(" ", m, " mg") else "" }, "")
  a <- ifelse(arm == "placebo", "Placebo", paste0("LDN", dose))
  paste0(t, "\n", a, " (day ", day, ")")
}

# D50 values at which the adopted dose-dependent PK's refit sits on its parameter boundary (fp0 or floss0 -> 1), read from
# pk_dose_model_comparison.csv (R4 verification: such a refit has a near-singular Hessian and must be flagged, not shown
# as an ordinary sensitivity).
d50_boundary <- function(res_dir, dose_model) {
  f <- file.path(res_dir, "pk_dose_model_comparison.csv")
  if (!file.exists(f)) stop("missing input: ", f, call. = FALSE)
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d$fp_D50_mg[d$dose_model == dose_model & as.logical(d$at_boundary) %in% TRUE & is.finite(d$fp_D50_mg)]
}
