options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
stage <- readLines("output/results/binding_stage.txt")[1]
cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
tt <- read_test_data("data/test/trials_test.csv", "B")
tt$se <- ifelse(is.na(tt$se), (tt$hi95 - tt$lo95) / 3.92, tt$se)
# INNOVA's change_nrs rows were derived by hand as Table 3 level minus baseline (review I11). Recompute them from the
# curated level file, read through read_test_data so it is hashed under stage T, and stop on any disagreement.
innova_chk <- check_innova_change(read_test_data("data/test/innova_nrs_levels.csv", "B"), tt)
write.csv(innova_chk, "output/results/innova_change_derivation.csv", row.names = FALSE)
cat(sprintf("INNOVA change rows re-derived from innova_nrs_levels.csv: %d rows agree within rounding\n", nrow(innova_chk)))
# Primary endpoint: exactly ONE scored point per trial (INNOVA's three days come from one correlated LMM and
# must not be summed as independent observations). FINAL: day 84. INNOVA: day 90 (its primary endpoint).
PRIMARY_DAY <- c(final = 84, innova = 90)
trial_key <- function(tr) { k <- tolower(tr); if (!k %in% names(PRIMARY_DAY)) stop("no primary endpoint day defined for trial '", tr, "'"); k }
EST_UNADJ <- "unadjusted difference in change from baseline (LDN minus placebo; from the two arm change_nrs rows)"
res <- list(); warns <- character(0); dstore <- list()
note <- function(...) warns <<- c(warns, sprintf(...))  # collected, emitted once after the hypothesis loop
# One predictive seed per TRIAL (R2 fix round 1, Minor 4): different trials get INDEPENDENT trial-level draws (the
# new-trial placebo ceiling and the shared model-structure residual of review I7), while the days of ONE trial
# (INNOVA d90/180/365) share that trial's draw, as the same patients do. The seeds are the same for every
# hypothesis, so the draws are paired across hypotheses (used for the contrast MCSE below).
pred_seed <- local({
  key <- unique(tolower(tt$trial))
  function(tr, dy) 7L + match(tolower(tr), key)
})
interval_seed <- local({
  key <- unique(tolower(tt$trial))
  function(tr, dy, endpoint) 900000L + 10000L * match(tolower(tr), key) + 10L * as.integer(dy) +
    match(sub(":.*$", "", endpoint), c("change_nrs", "diff_nrs", "diff_nrs_adjusted_level")) +
    if (grepl("placebo$", endpoint)) 5L else 0L
})
for (h in c("H0", "H1", "H2", "H3")) {
  fit <- readRDS(sprintf("output/results/posterior_%s.rds", h))
  for (tr in unique(tt$trial)) for (dy in unique(tt$day[tt$trial == tr])) {
    arms <- tt[tt$trial == tr & tt$day == dy & tt$endpoint == "change_nrs", ]
    if (nrow(arms) == 0) {
      note("trial %s day %s: skipped, no change_nrs arm rows", tr, dy); next }
    p <- predict_trial(fit, arms, cache, seed = pred_seed(tr, dy))  # one shared residual draw per trial-draw
    pm <- split(p$pred, p$row)
    mk <- function(endpoint, obs, draws, lpd, estimand, role, obs_se = NA_real_, se_note = "") {
      # keep the draws so the Monte Carlo SE of each elpd can be formed after the loop (review I8)
      dstore[[length(dstore) + 1L]] <<- list(hypothesis = h, trial = tr, day = dy, endpoint = endpoint,
                                             role = role, obs = obs, se = obs_se, draws = draws)
      # R3a item 1: the draws describe the TRUE value (the trial-shared residual cancels in a difference), so the
      # observed value is compared with the posterior predictive interval for the OBSERVED value (draws plus
      # N(0, obs_se)); the true-value interval is kept as true_lo90/true_hi90. Rows without an obs_se (INNOVA's
      # derived arm changes) get no predictive interval and inside90 = NA. The observation-error seed depends on
      # (trial, day, endpoint) only, so the hypotheses share it.
      iv <- heldout_intervals(obs, draws, obs_se, seed = interval_seed(tr, dy, endpoint))
      data.frame(hypothesis = h, trial = tr, endpoint = endpoint, day = dy, observed = obs,
                 pred_median = stats::median(draws), pred_lo90 = iv$pred_lo90, pred_hi90 = iv$pred_hi90,
                 true_lo90 = iv$true_lo90, true_hi90 = iv$true_hi90, inside90 = iv$inside90, lpd = lpd,
                 lpd_mcse = if (is.finite(obs_se)) elpd_mcse(obs, cbind(draws), obs_se)$mcse else NA_real_,
                 estimand = estimand, role = role, obs_se = obs_se, se_note = se_note,
                 interval_definition = iv$interval_definition)
    }
    # Arm rows with se = NA are kept (figures can show them) but get lpd = NA and stay out of the lpd sums.
    rows <- lapply(seq_len(nrow(arms)), function(i)
      mk(paste0("change_nrs:", arms$arm[i]), arms$mean[i], pm[[i]],
         if (is.finite(arms$se[i])) lpd_normal(arms$mean[i], pm[[i]], arms$se[i]) else NA_real_,
         arms$estimand[i], "arm", arms$se[i], ""))
    d <- tt[tt$trial == tr & tt$day == dy & tt$endpoint == "diff_nrs", ]
    ldn <- which(arms$arm != "placebo"); pbo <- which(arms$arm == "placebo")
    is_prim <- dy == PRIMARY_DAY[[trial_key(tr)]]
    if (nrow(d) != 1) {
      note("trial %s day %s: diff_nrs not scored, %d diff_nrs rows (need exactly 1)", tr, dy, nrow(d))
    } else if (length(ldn) != 1 || length(pbo) != 1) {
      note("trial %s day %s: diff_nrs not scored, arms not simulable (%d LDN, %d placebo arm rows; need 1 each)",
           tr, dy, length(ldn), length(pbo))
    } else {
      dd <- pm[[ldn]] - pm[[pbo]]  # same draw index in both arms (predicted difference in CHANGE)
      if (!is.finite(d$se)) note("trial %s day %s: diff_nrs predicted but lpd = NA, no finite se or CI on the diff row", tr, dy)
      role <- if (is_prim) "primary" else "secondary"
      if (grepl("level", d$estimand, ignore.case = TRUE)) {
        # Estimand mismatch: the reported row is an LMM-adjusted difference in LEVEL, the model predicts a difference in
        # CHANGE. Score against the like-for-like UNADJUSTED difference in change (LDN change minus placebo change,
        # from the two arm rows). Point estimate is unadjusted; the SE is borrowed from the adjusted row's 95% CI.
        obs_u <- arms$mean[ldn] - arms$mean[pbo]
        se_note <- "point estimate unadjusted (difference of arm changes); SE borrowed from the LMM-adjusted 95% CI of the level difference (approximation)"
        rows[[length(rows) + 1]] <- mk("diff_nrs", obs_u, dd, if (is.finite(d$se)) lpd_normal(obs_u, dd, d$se) else NA_real_,
                                       EST_UNADJ, role, d$se, se_note)
        rows[[length(rows) + 1]] <- mk("diff_nrs_adjusted_level", d$mean, dd,
                                       if (is.finite(d$se)) lpd_normal(d$mean, dd, d$se) else NA_real_,
                                       paste0("SENSITIVITY: ", d$estimand, " (model predicts a difference in change)"),
                                       "sensitivity", d$se, "as reported (adjusted, level); prediction is a difference in change")
      } else {
        rows[[length(rows) + 1]] <- mk("diff_nrs", d$mean, dd, if (is.finite(d$se)) lpd_normal(d$mean, dd, d$se) else NA_real_,
                                       d$estimand, role, d$se, "as reported")
      }
    }
    res[[length(res) + 1]] <- do.call(rbind, rows)
  }
}
for (w in unique(warns)) warning(w, call. = FALSE)
out <- do.call(rbind, res)   # inside90 is set per row by heldout_intervals() against the PREDICTIVE interval
write.csv(out, "output/results/test_predictions.csv", row.names = FALSE)

# elpd per hypothesis. PRIMARY = one endpoint per trial (FINAL d84 + INNOVA d90, unadjusted difference in change).
prim_txt <- paste(sprintf("%s d%d", names(PRIMARY_DAY), PRIMARY_DAY), collapse = " + ")
# R2 item 6 (review I8): every elpd also gets its MONTE CARLO standard error, formed from the variability across
# the posterior-predictive draws it was computed from. It is distinct from the observation-level SE (n = 2 held-out
# endpoints) that test_elpd_contrasts.csv reports as se_diff: the Monte Carlo SE says how far the number would move
# with different predictive draws from the same posterior; the observation-level SE says how far it would move with
# different trials. Neither covers a refit with different MCMC chains or different inputs.
set_elpd <- function(h, pick) {
  rs <- Filter(function(r) r$hypothesis == h && is.finite(r$se) && pick(r), dstore)
  if (!length(rs)) return(list(elpd = 0, mcse = NA_real_, n_draws = NA_integer_, n_obs = 0L))
  elpd_mcse(vapply(rs, function(r) r$obs, 0), do.call(cbind, lapply(rs, function(r) r$draws)),
            vapply(rs, function(r) r$se, 0))
}
is_primary_diff <- function(r) r$endpoint == "diff_nrs" && r$role == "primary"
elpd <- do.call(rbind, lapply(c("H0", "H1", "H2", "H3"), function(h) {
  o <- out[out$hypothesis == h, ]
  sm <- function(x) sum(x$lpd, na.rm = TRUE)
  pr <- set_elpd(h, is_primary_diff)
  al <- set_elpd(h, function(r) r$endpoint == "diff_nrs")
  sn <- set_elpd(h, function(r) (is_primary_diff(r) && !grepl("^innova$", r$trial, ignore.case = TRUE)) ||
                   (r$endpoint == "diff_nrs_adjusted_level" && r$role == "sensitivity" && r$day == PRIMARY_DAY[["innova"]]))
  data.frame(hypothesis = h,
    elpd_primary = sm(o[o$endpoint == "diff_nrs" & o$role == "primary", ]),
    mcse_elpd_primary = pr$mcse,
    primary_endpoints = prim_txt,
    n_primary_points = sum(o$endpoint == "diff_nrs" & o$role == "primary" & is.finite(o$lpd)),
    n_draws_mcse = pr$n_draws,
    elpd_all_days_secondary_nonindependent = sm(o[o$endpoint == "diff_nrs", ]),
    mcse_elpd_all_days_secondary_nonindependent = al$mcse,
    all_days_note = "SECONDARY, NON-INDEPENDENT: sums every scored day; INNOVA days 90/180/365 are correlated, so this overweights INNOVA. Not for model selection.",
    elpd_sensitivity_adjusted_innova = sm(o[(o$endpoint == "diff_nrs" & o$role == "primary" & !grepl("^innova$", o$trial, ignore.case = TRUE)) |
                                             (o$endpoint == "diff_nrs_adjusted_level" & o$role == "sensitivity" & o$day == PRIMARY_DAY[["innova"]]), ]),
    mcse_elpd_sensitivity_adjusted_innova = sn$mcse,
    sensitivity_note = "SENSITIVITY: primary endpoints with INNOVA d90 scored against the LMM-adjusted level difference (estimand mismatch)",
    mcse_note = "mcse_* columns are MONTE CARLO (posterior-draw) standard errors of the elpd itself, from the variability across the draws used; they are NOT the observation-level SE across held-out endpoints (see se_diff in test_elpd_contrasts.csv) and do not cover run-to-run MCMC variation.",
    check_elpd_primary = pr$elpd)
}))
stopifnot(all(elpd$n_primary_points == length(PRIMARY_DAY)))
stopifnot(all(abs(elpd$elpd_primary - elpd$check_elpd_primary) < 1e-8))   # the MCSE uses the SAME draws
elpd$check_elpd_primary <- NULL
write.csv(elpd, "output/results/test_elpd.csv", row.names = FALSE)

# Monte Carlo SE of each pairwise PRIMARY contrast from the per-draw PAIRED difference (R2 fix round 1, Minor 3).
# All hypotheses use the same per-trial predictive seeds, so their draw errors are positively correlated and the
# independence formula sqrt(mcse_a^2 + mcse_b^2) is conservative; elpd_mcse_diff() is the paired estimate.
prim_recs <- function(h) {
  rs <- Filter(function(r) r$hypothesis == h && is.finite(r$se) && is_primary_diff(r), dstore)
  rs[order(vapply(rs, function(r) paste(tolower(r$trial), r$day), ""))]
}
hh <- c("H0", "H1", "H2", "H3"); cmc <- list()
for (i in seq_along(hh)) for (j in seq_along(hh)) if (i < j) {
  ra <- prim_recs(hh[i]); rb <- prim_recs(hh[j])
  stopifnot(identical(vapply(ra, function(r) paste(r$trial, r$day), ""), vapply(rb, function(r) paste(r$trial, r$day), "")))
  m <- elpd_mcse_diff(vapply(ra, function(r) r$obs, 0), do.call(cbind, lapply(ra, function(r) r$draws)),
                      do.call(cbind, lapply(rb, function(r) r$draws)), vapply(ra, function(r) r$se, 0))
  cmc[[length(cmc) + 1]] <- data.frame(hypothesis_a = hh[i], hypothesis_b = hh[j], elpd_diff_a_minus_b = m$diff,
    mcse_diff = m$mcse, n_draws = m$n_draws,
    mcse_note = "Monte Carlo SE of the contrast from the per-draw paired difference (same predictive seeds across hypotheses)")
}
write.csv(do.call(rbind, cmc), "output/results/test_elpd_contrast_mcse.csv", row.names = FALSE)
print(elpd[, c("hypothesis", "elpd_primary", "mcse_elpd_primary",
               "elpd_all_days_secondary_nonindependent", "elpd_sensitivity_adjusted_innova")])
# INNOVA check: every hypothesis should have all three days scored (diff_nrs) and shown (arms)
inn <- out[grepl("innova", out$trial, ignore.case = TRUE) & out$endpoint == "diff_nrs", ]
cat("INNOVA diff_nrs days scored per hypothesis:\n"); print(tapply(inn$day, inn$hypothesis, function(x) paste(sort(x), collapse = ",")))
