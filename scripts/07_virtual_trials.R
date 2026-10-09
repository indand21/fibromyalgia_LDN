options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
stage <- readLines("output/results/binding_stage.txt")[1]
cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
tt <- read_test_data("data/test/trials_test.csv", "B")
om_all <- read.csv("output/results/omegas.csv")
if (!file.exists("output/results/omegas_sd_sensitivity.csv"))
  stop("missing output/results/omegas_sd_sensitivity.csv; run scripts/05_calibrate_variability.R first (R3a item 4)")
om_sd <- read.csv("output/results/omegas_sd_sensitivity.csv")
om_sd <- om_sd[startsWith(om_sd$calibration, "SENSITIVITY"), ]
N_DRAWS <- 200; N_REP <- 1000          # main virtual trials (unchanged)
EFFECT_THRESHOLD_PP <- 0.5             # magnitude guard for power and Sobol (review I2): day-84 net effect, pp
RESP_COHORT <- 2000; RESP_REP <- 50; RESP_SEED <- 4242; RESP_BINOM <- 20   # responder scoring (review I6)
POWER_NS <- c(50, 100, 150, 200, 300, 500); POWER_REP <- 200; POWER_DRAWS <- 50
HYPS <- c("H0", "H1", "H2", "H3")
# Run-driver switch (R3b): V5_SCRIPT07_MODE=summary, used by the binding-stage sensitivity passes, writes the virtual-trial
# summary and the responder scores only; it SKIPS the power curves and the SD-sensitivity calibration runs.
MODE <- Sys.getenv("V5_SCRIPT07_MODE", "full")
if (!MODE %in% c("full", "summary")) stop("V5_SCRIPT07_MODE must be 'full' or 'summary' (got '", MODE, "')")
SUMMARY_MODE <- MODE == "summary"
tt$se <- ifelse(is.na(tt$se), (tt$hi95 - tt$lo95) / 3.92, tt$se)
ch <- tt[tt$endpoint == "change_nrs", ]
bounds <- read_baseline_bounds("data/params/baseline_bounds.csv")   # per-trial enrolment bounds (R2 fix round 1)
designs <- lapply(split(ch, ch$trial), function(d) {
  ldn <- d[d$arm != "placebo", ][1, ]; pbo <- d[d$arm == "placebo", ][1, ]
  list(trial = ldn$trial, schedule_ldn = ldn$schedule, day = min(d$day), n_ldn = ldn$n, n_pbo = pbo$n,
       base_mean = mean(d$baseline_nrs), base_sd = mean(d$baseline_sd),
       base_bounds = baseline_bounds_for(bounds, ldn$trial))
})
for (dz in designs) stopifnot(is.finite(c(dz$base_mean, dz$base_sd, dz$n_ldn, dz$n_pbo, dz$day)))
fits <- lapply(stats::setNames(nm = HYPS), function(h) readRDS(sprintf("output/results/posterior_%s.rds", h)))
om_of <- function(tab, h) { r <- tab[tab$hypothesis == h, ]
  if (nrow(r) != 1) stop("expected one omega row for ", h, " (found ", nrow(r), ")"); as.list(r) }
oms <- lapply(stats::setNames(nm = HYPS), function(h) om_of(om_all, h))
oms_sd <- lapply(stats::setNames(nm = HYPS), function(h) om_of(om_sd, h))

# ---- magnitude guard (review I2): which hypotheses have an effect worth a power curve or a Sobol analysis ----
el <- effect_eligibility(fits, cache, day = 84, threshold_pp = EFFECT_THRESHOLD_PP)
write.csv(el, "output/results/effect_eligibility.csv", row.names = FALSE)
print(el[, c("hypothesis", "net_effect_pp", "threshold_pp", "eligible", "power_role")])
eligible <- el$hypothesis[el$eligible]
if (!length(eligible)) message("No hypothesis exceeds the ", EFFECT_THRESHOLD_PP, " pp threshold: the power file holds only the H0 type-I-error check.")

# ---- every independent simulation goes to one PSOCK cluster (R3a item 2) ----
# Each task carries a posterior thinned to its own draws and its own seed; run_vt_task installs a Mersenne-Twister
# stream from that seed, so the numbers are identical to a sequential run and do not depend on the worker count.
vt_task <- function(h, dz, om, seed, n_rep, n_draws, meta)
  list(fit = thin_fit(fits[[h]], n_draws), design = dz, cache = cache, om = om, n_rep = n_rep, n_draws = n_draws,
       seed = seed, meta = meta)
tasks <- list()
for (h in HYPS) for (dz in designs) for (mult in c(1, 2)) {           # main virtual trials, as before (seed 99)
  om2 <- utils::modifyList(oms[[h]], list(omega_p = oms[[h]]$omega_p * mult, omega_e = oms[[h]]$omega_e * mult))
  tasks[[length(tasks) + 1]] <- vt_task(h, dz, om2, 99, N_REP, N_DRAWS, list(set = "main", h = h, trial = dz$trial, mult = mult))
}
if (!SUMMARY_MODE) for (h in HYPS) for (dz in designs)                 # SD sensitivity (review I10), same seed: paired
  tasks[[length(tasks) + 1]] <- vt_task(h, dz, oms_sd[[h]], 99, N_REP, N_DRAWS, list(set = "sd_sens", h = h, trial = dz$trial, mult = 1))
for (h in HYPS) for (dz in designs) {                                  # responder RATES from a large cohort (review I6)
  dz_big <- utils::modifyList(dz, list(n_ldn = RESP_COHORT, n_pbo = RESP_COHORT))
  tasks[[length(tasks) + 1]] <- vt_task(h, dz_big, oms[[h]], RESP_SEED, RESP_REP, N_DRAWS, list(set = "resp", h = h, trial = dz$trial, mult = 1))
}
n_workers <- default_n_workers()
timing <- list()
run_jobs <- function(map) {
  # largest jobs first for the load balancer; results go back to their own slots, so order cannot change them
  ord <- order(-vapply(tasks, function(t) t$n_rep * (t$design$n_ldn + t$design$n_pbo), 0))
  t0 <- Sys.time(); out <- vector("list", length(tasks)); out[ord] <- map(tasks[ord], run_vt_task)
  timing$virtual <<- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  roles <- stats::setNames(el$power_role, el$hypothesis)
  t0 <- Sys.time(); pw <- pw_sd <- NULL
  if (!SUMMARY_MODE) {
  pw <- power_curves(fits, oms, designs, cache, hyps = eligible, ns = POWER_NS, n_rep = POWER_REP, n_draws = POWER_DRAWS,
                     map = map, roles = roles)
  pw_sd <- power_curves(fits, oms_sd, designs, cache, hyps = eligible, ns = POWER_NS, n_rep = POWER_REP, n_draws = POWER_DRAWS,
                        map = map, roles = roles)
  }
  timing$power <<- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  list(vt = out, pw = pw, pw_sd = pw_sd)
}
t_all <- Sys.time()
R <- if (n_workers > 1) with_cluster(n_workers, function(cl) run_jobs(function(X, FUN) parallel::parLapplyLB(cl, X, FUN))) else
  run_jobs(lapply)
cat(sprintf("script 07 simulations on %d worker(s): virtual trials %.1f min, power %.1f min, total %.1f min\n", n_workers,
            timing$virtual, timing$power, as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
set_of <- vapply(tasks, function(t) t$meta$set, "")
tag <- function(v, meta, om) { v$hypothesis <- meta$h; v$trial <- meta$trial; v$omega_mult <- meta$mult
  v$omega_p_used <- om$omega_p; v$omega_e_used <- om$omega_e; v }

# ---- main virtual trials ----
vt <- list(); clip_rows <- list()
for (i in which(set_of == "main")) {
  m <- tasks[[i]]$meta; v <- tag(R$vt[[i]], m, tasks[[i]]$om); dz <- designs[[m$trial]]
  # how often the physical pain_end clip at 0 / 10 binds, per arm (review I3: report, do not hide)
  cl <- attr(v, "pain_end_clip"); bp <- attr(v, "baseline_parent")
  clip_rows[[length(clip_rows) + 1]] <- do.call(rbind, lapply(split(cl, cl$arm), function(z) data.frame(
    hypothesis = m$h, trial = m$trial, omega_mult = m$mult, arm = z$arm[1],
    clip_low_mean = mean(z$low), clip_low_max = max(z$low), clip_high_mean = mean(z$high), clip_high_max = max(z$high),
    base_bounds = sprintf("[%g, %g]", dz$base_bounds[1], dz$base_bounds[2]),
    base_reported = sprintf("%.3f/%.3f", dz$base_mean, dz$base_sd),
    base_parent = sprintf("%.3f/%.3f", bp$mu, bp$sigma), base_parent_matched = bp$matched)))
  vt[[length(vt) + 1]] <- v
}
vt <- do.call(rbind, vt); write.csv(vt, "output/results/virtual_trials.csv", row.names = FALSE)
write.csv(do.call(rbind, clip_rows), "output/results/virtual_pain_end_clip.csv", row.names = FALSE)

# ---- summary: per hypothesis x trial x omega multiplier ----
obs_rate <- function(d, cut, pbo, trial) {  # observed responder proportion for an arm at the design day, NA (with warning) if not reported
  ep <- paste0("resp", cut)
  r <- d[d$endpoint == ep & ((d$arm == "placebo") == pbo), ]
  if (nrow(r) == 0) { warning(trial, ": observed endpoint ", ep, " (", if (pbo) "placebo" else "LDN", ") is missing; recording NA"); return(NA_real_) }
  m <- r$mean[1]; if (m > 1) m / 100 else m
}
summarise_vt <- function(vt, om_label_fun) {
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
    # Review M3: the old p_tail was one-sided toward the sign of the observation, so FINAL (obs < 0) and INNOVA
    # (obs > 0) were not comparable. It is kept under an explicit name next to the two-sided value.
    pv <- ppp_values(v$diff, obs)
    row <- data.frame(hypothesis = key$hypothesis, trial = key$trial, omega_mult = key$omega_mult,
      omega_p_used = v$omega_p_used[1], omega_e_used = v$omega_e_used[1], omega_label = om_label_fun(key$omega_mult),
      n_draws = N_DRAWS, n_rep = N_REP, obs_diff = obs, obs_estimand = obs_est,
      obs_diff_adjusted_level_sensitivity = obs_adj, se_obs = se_obs,
      diff_median = stats::median(v$diff), diff_lo90 = unname(stats::quantile(v$diff, 0.05)),
      diff_hi90 = unname(stats::quantile(v$diff, 0.95)),
      ppp_one_sided_toward_sign_of_obs = pv$p_one_sided_toward_sign_of_obs, ppp_two_sided = pv$p_two_sided,
      ppp_note = paste("ppp_one_sided_toward_sign_of_obs = P(sim diff >= obs) if obs >= 0, else P(sim diff <= obs) (the former",
                       "p_tail; its direction differs between trials); ppp_two_sided = 2 x the smaller tail, capped at 1"))
    for (cut in c(15, 30, 50)) {
      row[[paste0("resp", cut, "_ldn_median")]] <- stats::median(v[[paste0("resp", cut, "_ldn")]])
      row[[paste0("resp", cut, "_pbo_median")]] <- stats::median(v[[paste0("resp", cut, "_pbo")]])
      row[[paste0("resp", cut, "_ldn_obs")]] <- obs_rate(d, cut, FALSE, key$trial)
      row[[paste0("resp", cut, "_pbo_obs")]] <- obs_rate(d, cut, TRUE, key$trial)
    }
    summ[[length(summ) + 1]] <- row
  }
  do.call(rbind, summ)
}
# omega_p is of order 2 (0-100 scale) but omega_e is only 1e-5 to 3e-3, so multiplying both by 2 effectively doubles the
# PLACEBO-RESPONSE variability term; the drug-effect variability term changes by a negligible absolute amount.
write.csv(summarise_vt(vt, function(mult) if (mult == 1) "as calibrated (script 05)" else
  "sensitivity: both terms x2, which effectively doubles the placebo-response variability term (omega_p ~2); omega_e ~1e-5 to 3e-3 so the drug-effect variability change is negligible"),
  "output/results/virtual_summary.csv", row.names = FALSE)

# ---- SD sensitivity (review I10): virtual trials and responder rates under the Bested et al. 2023 SDs ----
if (!SUMMARY_MODE) {
vt_sd <- do.call(rbind, lapply(which(set_of == "sd_sens"), function(i) tag(R$vt[[i]], tasks[[i]]$meta, tasks[[i]]$om)))
sum_sd <- summarise_vt(vt_sd, function(mult) "SENSITIVITY: omegas calibrated to the SDs quoted by Bested et al. 2023 (14.6 placebo, 12.5 LDN); see omegas_sd_sensitivity.csv")
prim <- read.csv("output/results/virtual_summary.csv"); prim <- prim[prim$omega_mult == 1, ]
cmpc <- c("diff_median", "diff_lo90", "diff_hi90", paste0("resp", rep(c(15, 30, 50), each = 2), c("_ldn_median", "_pbo_median")))
pm <- prim[, c("hypothesis", "trial", cmpc)]; names(pm)[-(1:2)] <- paste0(cmpc, "_primary")
sum_sd <- merge(sum_sd, pm, by = c("hypothesis", "trial"))
write.csv(sum_sd, "output/results/virtual_summary_sd_sensitivity.csv", row.names = FALSE)
}

# ---- responder scoring (review I6) ----
# Rates come from a LARGE simulated cohort per replicate trial (RESP_COHORT per arm, RESP_REP replicates per posterior
# draw), so each (draw, replicate) pair is the population responder probability of one replicate trial; counts are then
# sampled Binomial(actual arm n, rate). (The main virtual trials simulate the actual n, so their per-replicate rates
# already contain the binomial sampling once; applying a binomial to them would count it twice. They are used only for
# the *_direct_check columns, scored as counts round(rate x n).)
resp_rows <- list(); ri <- 0
for (i in which(set_of == "resp")) {
  m <- tasks[[i]]$meta; big <- R$vt[[i]]; dz <- designs[[m$trial]]
  oc <- observed_responder_counts(tt, m$trial, dz$day)
  main <- vt[vt$hypothesis == m$h & vt$trial == m$trial & vt$omega_mult == 1, ]
  for (j in seq_len(nrow(oc))) {
    ri <- ri + 1; cut <- oc$cut[j]
    s <- responder_scores(big[[paste0("resp", cut, "_ldn")]], big[[paste0("resp", cut, "_pbo")]],
                          oc$k_ldn[j], oc$n_ldn[j], oc$k_pbo[j], oc$n_pbo[j], seed = RESP_SEED + 100 * match(m$h, HYPS) + ri,
                          n_binom = RESP_BINOM)
    kd <- round(main[[paste0("resp", cut, "_ldn")]] * dz$n_ldn) / dz$n_ldn - round(main[[paste0("resp", cut, "_pbo")]] * dz$n_pbo) / dz$n_pbo
    resp_rows[[ri]] <- cbind(data.frame(hypothesis = m$h, trial = m$trial, day = dz$day, cut_pct = cut), s,
      data.frame(pred_diff_median_direct_check = stats::median(kd), ppp_two_sided_direct_check = ppp_values(kd, s$obs_diff)$p_two_sided,
                 rate_source = sprintf("%d posterior draws x %d replicate trials, %d simulated subjects per arm each", N_DRAWS, RESP_REP, RESP_COHORT),
                 definition = paste("EXPLORATORY (added after the held-out results were seen; review I6).",
                                    "pp_prob_*: posterior predictive probability of the observed count(s) under Binomial(actual n, rate);",
                                    "ppp_one_sided_upper: P(sim LDN-minus-placebo difference in proportions >= observed), the direction of an",
                                    "LDN benefit for every trial; ppp_two_sided: 2 x the smaller tail; fisher_*: the source trial's own",
                                    "Fisher exact test on the observed counts (source-trial test, descriptive; not a model score)")))
  }
}
resp_scores <- do.call(rbind, resp_rows)
write.csv(resp_scores, "output/results/responder_scores.csv", row.names = FALSE)
print(resp_scores[, c("hypothesis", "trial", "cut_pct", "obs_diff", "pred_diff_median", "pp_prob_counts_joint",
                      "ppp_one_sided_upper", "ppp_two_sided", "fisher_p_two_sided")], digits = 3)

# ---- rejection-rate curves (review I2): genuine power for each ELIGIBLE hypothesis, H0 as the type-I check ----
if (!SUMMARY_MODE) {
pw <- R$pw; pw$calibration <- "primary"
write.csv(pw, "output/results/power_curve.csv", row.names = FALSE)
h0 <- pw[pw$hypothesis == "H0", ]
write.csv(rbind(
  data.frame(hypothesis = "H0", trial = h0$trial, n_per_arm = h0$n_per_arm, variant = "with_model_unc",
             rej_rate_null = h0$rej_rate_with_model_unc, crit_value = h0$crit_diff_with_model_unc, n_reps = h0$n_reps),
  data.frame(hypothesis = "H0", trial = h0$trial, n_per_arm = h0$n_per_arm, variant = "sampling_only (sigma_rel = 0)",
             rej_rate_null = h0$rej_rate_sampling_only, crit_value = h0$crit_diff_sampling_only, n_reps = h0$n_reps)) |>
  transform(definition = "no drug effect (H0); critical value = empirical 5% quantile from an independent seed; expected about 0.05"),
  "output/results/power_null_check.csv", row.names = FALSE)
pw_sd <- R$pw_sd; pw_sd$calibration <- "SENSITIVITY: omegas calibrated to the Bested et al. 2023 SDs (14.6 placebo, 12.5 LDN)"
write.csv(pw_sd, "output/results/power_curve_sd_sensitivity.csv", row.names = FALSE)
}

# ---- downsampled virtual trials (the full file is ~500 MB and is not committed): fixed-seed random 1% sample ----
set.seed(20260930); keep <- sample.int(nrow(vt), ceiling(0.01 * nrow(vt)))
write.csv(vt[sort(keep), ], "output/results/virtual_trials_sample.csv", row.names = FALSE)
