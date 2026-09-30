options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
stage <- readLines("output/results/binding_stage.txt")[1]
cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
tt <- read_test_data("data/test/trials_test.csv", "B")
tt$se <- ifelse(is.na(tt$se), (tt$hi95 - tt$lo95) / 3.92, tt$se)
# Primary endpoint: exactly ONE scored point per trial (INNOVA's three days come from one correlated LMM and
# must not be summed as independent observations). FINAL: day 84. INNOVA: day 90 (its primary endpoint).
PRIMARY_DAY <- c(final = 84, innova = 90)
trial_key <- function(tr) { k <- tolower(tr); if (!k %in% names(PRIMARY_DAY)) stop("no primary endpoint day defined for trial '", tr, "'"); k }
EST_UNADJ <- "unadjusted difference in change from baseline (LDN minus placebo; from the two arm change_nrs rows)"
res <- list(); warns <- character(0)
note <- function(...) warns <<- c(warns, sprintf(...))  # collected, emitted once after the hypothesis loop
for (h in c("H0", "H1", "H2", "H3")) {
  fit <- readRDS(sprintf("output/results/posterior_%s.rds", h))
  for (tr in unique(tt$trial)) for (dy in unique(tt$day[tt$trial == tr])) {
    arms <- tt[tt$trial == tr & tt$day == dy & tt$endpoint == "change_nrs", ]
    if (nrow(arms) == 0) {
      note("trial %s day %s: skipped, no change_nrs arm rows", tr, dy); next }
    p <- predict_trial(fit, arms, cache)  # each arm prediction already carries its own residual draw
    pm <- split(p$pred, p$row)
    mk <- function(endpoint, obs, draws, lpd, estimand, role, obs_se = NA_real_, se_note = "")
      data.frame(hypothesis = h, trial = tr, endpoint = endpoint, day = dy, observed = obs,
                 pred_median = stats::median(draws), pred_lo90 = unname(stats::quantile(draws, 0.05)),
                 pred_hi90 = unname(stats::quantile(draws, 0.95)), lpd = lpd, estimand = estimand, role = role,
                 obs_se = obs_se, se_note = se_note)
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
out <- do.call(rbind, res); out$inside90 <- out$observed >= out$pred_lo90 & out$observed <= out$pred_hi90
write.csv(out, "output/results/test_predictions.csv", row.names = FALSE)

# elpd per hypothesis. PRIMARY = one endpoint per trial (FINAL d84 + INNOVA d90, unadjusted difference in change).
prim_txt <- paste(sprintf("%s d%d", names(PRIMARY_DAY), PRIMARY_DAY), collapse = " + ")
elpd <- do.call(rbind, lapply(c("H0", "H1", "H2", "H3"), function(h) {
  o <- out[out$hypothesis == h, ]
  sm <- function(x) sum(x$lpd, na.rm = TRUE)
  data.frame(hypothesis = h,
    elpd_primary = sm(o[o$endpoint == "diff_nrs" & o$role == "primary", ]),
    primary_endpoints = prim_txt,
    n_primary_points = sum(o$endpoint == "diff_nrs" & o$role == "primary" & is.finite(o$lpd)),
    elpd_all_days_secondary_nonindependent = sm(o[o$endpoint == "diff_nrs", ]),
    all_days_note = "SECONDARY, NON-INDEPENDENT: sums every scored day; INNOVA days 90/180/365 are correlated, so this overweights INNOVA. Not for model selection.",
    elpd_sensitivity_adjusted_innova = sm(o[(o$endpoint == "diff_nrs" & o$role == "primary" & !grepl("^innova$", o$trial, ignore.case = TRUE)) |
                                             (o$endpoint == "diff_nrs_adjusted_level" & o$role == "sensitivity" & o$day == PRIMARY_DAY[["innova"]]), ]),
    sensitivity_note = "SENSITIVITY: primary endpoints with INNOVA d90 scored against the LMM-adjusted level difference (estimand mismatch)")
}))
stopifnot(all(elpd$n_primary_points == length(PRIMARY_DAY)))
write.csv(elpd, "output/results/test_elpd.csv", row.names = FALSE)
print(elpd[, c("hypothesis", "elpd_primary", "elpd_all_days_secondary_nonindependent", "elpd_sensitivity_adjusted_innova")])
# INNOVA check: every hypothesis should have all three days scored (diff_nrs) and shown (arms)
inn <- out[grepl("innova", out$trial, ignore.case = TRUE) & out$endpoint == "diff_nrs", ]
cat("INNOVA diff_nrs days scored per hypothesis:\n"); print(tapply(inn$day, inn$hypothesis, function(x) paste(sort(x), collapse = ",")))
