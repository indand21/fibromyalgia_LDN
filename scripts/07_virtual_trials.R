options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
stage <- readLines("output/results/binding_stage.txt")[1]
cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
tt <- read_test_data("data/test/trials_test.csv", "B")
om_all <- read.csv("output/results/omegas.csv")
N_DRAWS <- 200; N_REP <- 1000  # reduce N_DRAWS to 100 if a hypothesis x trial takes > 15 min; recorded in the summary
tt$se <- ifelse(is.na(tt$se), (tt$hi95 - tt$lo95) / 3.92, tt$se)
ch <- tt[tt$endpoint == "change_nrs", ]
designs <- lapply(split(ch, ch$trial), function(d) {
  ldn <- d[d$arm != "placebo", ][1, ]; pbo <- d[d$arm == "placebo", ][1, ]
  list(trial = ldn$trial, schedule_ldn = ldn$schedule, day = min(d$day), n_ldn = ldn$n, n_pbo = pbo$n,
       base_mean = mean(d$baseline_nrs), base_sd = mean(d$baseline_sd))
})
for (dz in designs) stopifnot(is.finite(c(dz$base_mean, dz$base_sd, dz$n_ldn, dz$n_pbo, dz$day)))
vt <- list()
for (h in c("H0", "H1", "H2", "H3")) {
  fit <- readRDS(sprintf("output/results/posterior_%s.rds", h)); om <- as.list(om_all[om_all$hypothesis == h, ])
  for (dz in designs) for (mult in c(1, 2)) {
    om2 <- modifyList(om, list(omega_p = om$omega_p * mult, omega_e = om$omega_e * mult))
    v <- virtual_trials(fit, dz, cache, om2, n_rep = N_REP, n_draws = N_DRAWS); v$hypothesis <- h; v$trial <- dz$trial; v$omega_mult <- mult
    v$omega_p_used <- om2$omega_p; v$omega_e_used <- om2$omega_e
    vt[[length(vt) + 1]] <- v
  }
}
vt <- do.call(rbind, vt); write.csv(vt, "output/results/virtual_trials.csv", row.names = FALSE)

# ---- summary: per hypothesis x trial x omega multiplier ----
obs_rate <- function(d, cut, pbo, trial) {  # observed responder proportion for an arm at the design day, NA (with warning) if not reported
  ep <- paste0("resp", cut)
  r <- d[d$endpoint == ep & ((d$arm == "placebo") == pbo), ]
  if (nrow(r) == 0) { warning(trial, ": observed endpoint ", ep, " (", if (pbo) "placebo" else "LDN", ") is missing; recording NA"); return(NA_real_) }
  m <- r$mean[1]; if (m > 1) m / 100 else m
}
summ <- list()
keys <- unique(vt[c("hypothesis", "trial", "omega_mult")])
for (k in seq_len(nrow(keys))) {
  key <- keys[k, ]
  v <- vt[vt$hypothesis == key$hypothesis & vt$trial == key$trial & vt$omega_mult == key$omega_mult, ]
  dz <- designs[[key$trial]]; d <- tt[tt$trial == key$trial & tt$day == dz$day, ]; od <- d[d$endpoint == "diff_nrs", ][1, ]
  if (nrow(d[d$endpoint == "diff_nrs", ]) == 0) warning(key$trial, ": observed diff_nrs missing at day ", dz$day, "; observed columns are NA")
  # Same estimand as script 06: the model predicts a difference in CHANGE. Where the reported diff row is an LMM-adjusted
  # difference in LEVEL (INNOVA), compare with the unadjusted change difference from the two arm rows; the adjusted value is
  # carried only as a labelled sensitivity column. SE (used only for reporting) is the adjusted row's, an approximation.
  obs_adj <- NA_real_; obs_est <- od$estimand
  obs <- od$mean; se_obs <- od$se
  if (!is.na(od$estimand) && grepl("level", od$estimand, ignore.case = TRUE)) {
    ca <- d[d$endpoint == "change_nrs", ]
    stopifnot(sum(ca$arm != "placebo") == 1, sum(ca$arm == "placebo") == 1)
    obs_adj <- od$mean; obs <- ca$mean[ca$arm != "placebo"] - ca$mean[ca$arm == "placebo"]
    obs_est <- "unadjusted difference in change from baseline (arm change_nrs rows); SE borrowed from the adjusted 95% CI"
  }
  tail_p <- if (is.na(obs)) NA_real_ else if (obs >= 0) mean(v$diff >= obs) else mean(v$diff <= obs)
  # omega_p is of order 2 (0-100 scale) but omega_e is only 1e-5 to 3e-3, so multiplying both by 2 effectively doubles the
  # PLACEBO-RESPONSE variability term; the drug-effect variability term changes by a negligible absolute amount.
  om_lab <- if (key$omega_mult == 1) "as calibrated (script 05)" else
    "sensitivity: both terms x2, which effectively doubles the placebo-response variability term (omega_p ~2); omega_e ~1e-5 to 3e-3 so the drug-effect variability change is negligible"
  row <- data.frame(hypothesis = key$hypothesis, trial = key$trial, omega_mult = key$omega_mult,
    omega_p_used = v$omega_p_used[1], omega_e_used = v$omega_e_used[1], omega_label = om_lab,
    n_draws = N_DRAWS, n_rep = N_REP, obs_diff = obs, obs_estimand = obs_est,
    obs_diff_adjusted_level_sensitivity = obs_adj, se_obs = se_obs,
    diff_median = stats::median(v$diff), diff_lo90 = unname(stats::quantile(v$diff, 0.05)),
    diff_hi90 = unname(stats::quantile(v$diff, 0.95)),
    p_tail = tail_p)
  for (cut in c(15, 30, 50)) {
    row[[paste0("resp", cut, "_ldn_median")]] <- stats::median(v[[paste0("resp", cut, "_ldn")]])
    row[[paste0("resp", cut, "_pbo_median")]] <- stats::median(v[[paste0("resp", cut, "_pbo")]])
    row[[paste0("resp", cut, "_ldn_obs")]] <- obs_rate(d, cut, FALSE, key$trial)
    row[[paste0("resp", cut, "_pbo_obs")]] <- obs_rate(d, cut, TRUE, key$trial)
  }
  summ[[length(summ) + 1]] <- row
}
write.csv(do.call(rbind, summ), "output/results/virtual_summary.csv", row.names = FALSE)

# ---- rejection-rate curve for the best-supported hypothesis (highest PRIMARY elpd, one endpoint per trial) ----
# Rejection rates come from the EMPIRICAL distribution of the simulated difference (LDN minus placebo change) across
# replicates. A replicate rejects (one-sided) when its diff falls below the empirical 5% quantile of the no-drug
# reference distribution (posterior H0) simulated at the same n per arm, under the SAME variant, with an independent
# seed. Two variants: model-structure uncertainty included (sigma_rel as fitted) and sampling variability only
# (sigma_rel = 0). Under a no-drug hypothesis this is a type-I-error check; under a drug hypothesis it is power.
el <- read.csv("output/results/test_elpd.csv"); best <- el$hypothesis[which.max(el$elpd_primary)]
fit <- readRDS(sprintf("output/results/posterior_%s.rds", best)); om <- as.list(om_all[om_all$hypothesis == best, ])
fit0 <- readRDS("output/results/posterior_H0.rds"); om0 <- as.list(om_all[om_all$hypothesis == "H0", ])
no_shift <- function(f) { if ("log_sigma_rel" %in% colnames(f$draws)) f$draws[, "log_sigma_rel"] <- -Inf; f }  # exp(-Inf) = 0
DEF_UNC <- paste("rej_rate_with_model_unc: share of replicates (200 replicates x 50 posterior draws) whose simulated diff in NRS change",
                 "(LDN minus placebo) is below the empirical 5% quantile of the H0 (no-drug) diff distribution at the same n,",
                 "with model-structure uncertainty (arm-level sigma_rel as fitted) in both distributions; independent seeds")
DEF_SAMP <- paste("rej_rate_sampling_only: as rej_rate_with_model_unc but with sigma_rel = 0 in both distributions",
                  "(sampling variability and between-subject variability only)")
DEF_ALL <- paste(DEF_UNC, "|", DEF_SAMP, "| one-sided, nominal 0.05; under a no-drug hypothesis this is a type-I-error check, under a drug hypothesis it is power")
pw <- list(); chk <- list()
for (dz in designs) for (n in c(50, 100, 150, 200, 300, 500)) {
  d2 <- modifyList(dz, list(n_ldn = n, n_pbo = n))
  vr <- function(f, o, samp_only, seed) virtual_trials(if (samp_only) no_shift(f) else f, d2, cache, o, n_rep = 200, n_draws = 50, seed = seed)
  rr <- list()
  for (variant in c("unc", "samp")) {
    so <- variant == "samp"
    crit <- unname(stats::quantile(vr(fit0, om0, so, 1001)$diff, 0.05))     # empirical null critical value, seed 1001
    v <- vr(fit, om, so, 99)                                                # candidate hypothesis, seed 99
    rr[[variant]] <- list(rate = mean(v$diff < crit), crit = crit, n_reps = nrow(v), sd = stats::sd(v$diff), mean = mean(v$diff))
    nul <- vr(fit0, om0, so, 2002)                                          # sanity: H0 (no drug) vs the same critical value, fresh seed
    chk[[length(chk) + 1]] <- data.frame(hypothesis = "H0", trial = dz$trial, n_per_arm = n,
      variant = if (so) "sampling_only (sigma_rel = 0)" else "with_model_unc", rej_rate_null = mean(nul$diff < crit),
      crit_value = crit, n_reps = nrow(nul), definition = "no drug effect (H0); critical value = empirical 5% quantile from an independent seed; expected about 0.05")
  }
  pw[[length(pw) + 1]] <- data.frame(hypothesis = best, best_is_H0 = best == "H0", trial = dz$trial, n_per_arm = n,
    rej_rate_with_model_unc = rr$unc$rate, rej_rate_sampling_only = rr$samp$rate,
    crit_diff_with_model_unc = rr$unc$crit, crit_diff_sampling_only = rr$samp$crit,
    mean_diff_with_model_unc = rr$unc$mean, sd_diff_with_model_unc = rr$unc$sd, sd_diff_sampling_only = rr$samp$sd,
    n_reps = rr$unc$n_reps, definition = DEF_ALL)
}
write.csv(do.call(rbind, pw), "output/results/power_curve.csv", row.names = FALSE)
write.csv(do.call(rbind, chk), "output/results/power_null_check.csv", row.names = FALSE)
if (best == "H0") message("Best-supported hypothesis is H0: the curve reflects placebo-only effects (type I error, not LDN power).")

# ---- downsampled virtual trials (the full file is ~500 MB and is not committed): fixed-seed random 1% sample ----
set.seed(20260930); keep <- sample.int(nrow(vt), ceiling(0.01 * nrow(vt)))
write.csv(vt[sort(keep), ], "output/results/virtual_trials_sample.csv", row.names = FALSE)
