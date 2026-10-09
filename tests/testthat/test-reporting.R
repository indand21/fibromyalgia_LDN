# R3b: Table 2 labels by computation, the held-out contrasts with both SEs, and the sensitivity summary.

priors <- read_dataset(file.path("..", "..", "data", "priors", "placebo_priors.csv"))

# SD of a 1-D density known up to a constant on a grid
grid_sd <- function(x, logd) { w <- exp(logd - max(logd)); w <- w / sum(w); m <- sum(w * x); sqrt(sum(w * (x - m)^2)) }

test_that("the closed-form prior SDs match the priors actually used by log_prior()", {
  for (h in c("H0", "H1", "H2", "H3")) {
    p0 <- init_par(h, "t1")
    # the marginal prior of log_tau integrates eta out (the N(0, tau) density of eta integrates to 1), so the eta term is
    # removed for it; for every other parameter tau is held fixed and the eta term is a constant
    lp_marg <- function(p) log_prior(p, h, priors, "t1") - stats::dnorm(p[["eta_t1"]], 0, exp(p[["log_tau"]]), log = TRUE)
    for (v in setdiff(names(p0), "eta_t1")) {
      x <- seq(-20, 12, length.out = 4001)   # long left tail for the log of a half-normal
      lp <- vapply(x, function(z) { p <- p0; p[[v]] <- z; lp_marg(p) }, 0)
      expect_equal(grid_sd(x[is.finite(lp)], lp[is.finite(lp)]), prior_sd_on_sampler_scale(v, priors), tolerance = 2e-3, info = paste(h, v))
    }
  }
  # eta_<trial>: marginal of N(0, tau) with tau ~ half-normal(0.5)
  set.seed(1); tau <- abs(stats::rnorm(4e5, 0, 0.5)); eta <- stats::rnorm(4e5, 0, tau)
  expect_equal(stats::sd(eta), prior_sd_on_sampler_scale("eta_final", priors), tolerance = 0.01)
  expect_equal(prior_sd_on_sampler_scale("log_kpbo", priors), 0.00833 / 0.106, tolerance = 1e-6)
  expect_error(prior_sd_on_sampler_scale("log_new_thing", priors), "no prior SD defined")
})

test_that("every sampled parameter of every hypothesis has a prior SD (labels are never typed by hand)", {
  for (h in c("H0", "H1", "H2", "H3"))
    expect_true(all(is.finite(prior_sd_on_sampler_scale(names(init_par(h, c("younger2013", "final"))), priors))), info = h)
})

test_that("posterior parameters are labelled prior-determined at a posterior/prior SD ratio of 0.9 or more", {
  lab <- posterior_source_type(c(0.123, 0.08, 0.43, 0.899, 0.9), c(0.1206, 0.0786, 0.5, 1, 1))
  expect_equal(lab[1], "prior-determined (posterior SD / prior SD = 1.02)")
  expect_equal(lab[2], "prior-determined (posterior SD / prior SD = 1.02)")
  expect_match(lab[3], "^fitted \\(training; posterior SD / prior SD = 0.86\\)$")
  expect_match(lab[4], "^fitted")                      # 0.899 < 0.9 even though it prints as 0.90
  expect_match(lab[5], "^prior-determined")
})

tp_fake <- function(lpd) {   # lpd: named list hypothesis -> c(final, innova)
  do.call(rbind, lapply(names(lpd), function(h) data.frame(hypothesis = h, trial = c("final", "innova"), endpoint = "diff_nrs",
    day = c(84, 90), role = "primary", lpd = lpd[[h]])))
}

test_that("contrasts carry the observation-level SE (n = 2) and the Monte Carlo SE, and flag a collapsed SE", {
  tp <- tp_fake(list(H1 = c(-1.0, -1.2), H3 = c(-1.532, -1.609), H0 = c(-0.3, -0.4)))
  tp <- rbind(tp, transform(tp[1, ], role = "secondary", day = 180, lpd = -99))     # secondary rows are ignored
  cm <- data.frame(hypothesis_a = c("H0", "H0", "H1"), hypothesis_b = c("H1", "H3", "H3"), mcse_diff = c(0.04, 0.07, 0.05))
  ct <- elpd_contrasts(tp, cm)
  r <- ct[ct$hypothesis_a == "H1" & ct$hypothesis_b == "H3", ]
  d <- c(0.532, 0.409)
  expect_equal(r$elpd_diff_a_minus_b, sum(d))
  expect_equal(r$se_diff_observation_level, abs(d[1] - d[2]))                     # n = 2: sd * sqrt(2) = |d1 - d2|
  expect_equal(r$mcse_diff_monte_carlo, 0.05)
  expect_equal(r$ratio_diff_over_mcse, sum(d) / 0.05)
  expect_true(r$se_collapsed_flag); expect_match(r$flag_note, "COLLAPSED SE")
  expect_match(r$per_trial_lpd_diff, "final d84 \\+0.532; innova d90 \\+0.409")
  # differences of opposite sign: large SE, no flag
  ct2 <- elpd_contrasts(tp_fake(list(H0 = c(0, 0), H1 = c(-0.5, 0.4))))
  expect_false(ct2$se_collapsed_flag); expect_true(is.na(ct2$mcse_diff_monte_carlo))
  expect_match(ct2$note, "unstable in both directions")
  expect_error(elpd_contrasts(tp[, -6]), "lack column")
})

test_that("sensitivity_summary reports occupancy, elpd with MCSE, contrasts and responder p-values for every stage", {
  root <- normalizePath(tempfile("ssum_"), winslash = "/", mustWork = FALSE); res <- file.path(root, "output", "results")
  dir.create(res, recursive = TRUE); writeLines("J", file.path(res, "binding_stage.txt"))
  w <- function(d, f) utils::write.csv(d, f, row.names = FALSE)
  w(data.frame(stage = c("A", "A2", "J"), dose_mg = 4.5, mean_RO = c(0.93, 0.43, 0.76), peak_RO = NA, trough_RO = NA, adopted = c(FALSE, FALSE, TRUE)),
    file.path(res, "steady_state_occupancy_by_stage.csv"))
  w(data.frame(stage = "J", dose_mg = c(4.5, 4.5), uncertainty_source = c("binding_only", "binding_and_pk"), mean_RO_lo95 = c(0.61, 0.60),
               mean_RO_hi95 = c(0.87, 0.88)), file.path(res, "occupancy_uncertainty_J.csv"))
  fill <- function(dir, shift) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    w(data.frame(hypothesis = c("H0", "H1"), elpd_primary = c(-0.7, -1.8) + shift, mcse_elpd_primary = c(0.001, 0.04),
                 primary_endpoints = "final d84 + innova d90"), file.path(dir, "test_elpd.csv"))
    w(tp_fake(list(H0 = c(-0.3, -0.4) + shift, H1 = c(-0.9, -0.9) + shift)), file.path(dir, "test_predictions.csv"))
    w(data.frame(hypothesis_a = "H0", hypothesis_b = "H1", mcse_diff = 0.03), file.path(dir, "test_elpd_contrast_mcse.csv"))
    w(data.frame(hypothesis = "H1", trial = "final", cut_pct = 30, ppp_two_sided = 0.45 + shift, obs_diff = 0.15, pred_diff_median = 0.07),
      file.path(dir, "responder_scores.csv"))
    w(data.frame(hypothesis = "H1", trial = "final", omega_mult = c(1, 2), ppp_two_sided = c(0.9, 0.8), obs_diff = -0.4,
                 diff_median = -0.5, diff_lo90 = -1, diff_hi90 = 0), file.path(dir, "virtual_summary.csv"))
  }
  fill(res, 0); fill(file.path(res, "sensitivity_stage_A2"), -1); fill(file.path(res, "sensitivity_stage_A"), -0.5)
  s <- sensitivity_summary(root)
  expect_setequal(unique(s$stage), c("J", "A2", "A"))
  occ <- s[s$metric == "mean_steady_state_occupancy_4.5mg", ]
  expect_equal(occ$value[occ$stage == "A2"], 0.43); expect_equal(occ$lo[occ$stage == "J"], 0.60)
  el <- s[s$metric == "elpd_primary" & s$stage == "A2", ]
  expect_equal(el$value, c(-1.7, -2.8)); expect_equal(el$mcse_monte_carlo, c(0.001, 0.04))
  ct <- s[s$metric == "elpd_contrast_a_minus_b" & s$stage == "J", ]
  expect_equal(ct$value, 1.1); expect_equal(ct$se_observation_level, 0.1); expect_equal(ct$mcse_monte_carlo, 0.03)
  expect_equal(s$value[s$metric == "responder_ppp_two_sided" & s$stage == "A"], -0.05)
  expect_equal(nrow(s[s$metric == "virtual_diff_ppp_two_sided", ]), 3)            # omega_mult 1 only
  expect_match(unique(s$stage_role[s$stage == "J"]), "^PRIMARY")
  unlink(file.path(res, "sensitivity_stage_A", "responder_scores.csv"))
  expect_error(sensitivity_summary(root), "missing input .*sensitivity_stage_A/responder_scores.csv")
})

# ---- scripts 10 and 11 (R3b): adopted stage only, clear stops, no Jansen S1, required outputs ----
test_that("scripts 10 and 11 stop with a message naming a missing input", {
  repo <- normalizePath(file.path("..", ".."), winslash = "/")
  root <- normalizePath(tempfile("figtab_"), winslash = "/", mustWork = FALSE); dir.create(root)
  file.copy(file.path(repo, c("R", "scripts", "data")), root, recursive = TRUE)
  dir.create(file.path(root, "output", "results"), recursive = TRUE)
  for (s in c("10_figures.R", "11_tables.R")) {
    r <- run_step(s, root, file.path(root, paste0(s, ".log")))
    expect_false(r$status == 0L, info = s)
    expect_true(any(grepl("missing input file[(]s[)]:.*binding_stage.txt", readLines(r$log, warn = FALSE))), info = s)
  }
  writeLines("J", file.path(root, "output", "results", "binding_stage.txt"))
  r <- run_step("11_tables.R", root, file.path(root, "t.log"))
  expect_true(any(grepl("missing input file[(]s[)]:.*stage_decision.csv", readLines(r$log, warn = FALSE))))
})

test_that("figures and tables use the adopted stage's files, never Jansen's S1, and write every required output", {
  s10 <- paste(readLines(file.path("..", "..", "scripts", "10_figures.R")), collapse = "\n")
  s11 <- paste(readLines(file.path("..", "..", "scripts", "11_tables.R")), collapse = "\n")
  for (s in list(s10, s11)) {
    expect_match(s, 'readLines(file.path(res, "binding_stage.txt"))[1]', fixed = TRUE)
    expect_false(grepl("$S1_jansen|[[\"S1_jansen", s, fixed = TRUE))                      # never read for display
    expect_false(grepl('"S1_jansen"', s, fixed = TRUE))
  }
  for (f in c("fig1_schematic", "fig2_pk_pet", "fig3_training_fits", "fig4_heldout", "fig5_cytokines", "fig6_virtual_trials",
              "fig7_occupancy", "fig8_power", "figS1_sobol", "figS2_power_sd_sensitivity", "figure_notes.csv",
              "sobol_seed_agreement.csv", "occupancy_uncertainty_%s.csv", "steady_state_occupancy_by_stage.csv", "true_lo90"))
    expect_true(grepl(f, s10, fixed = TRUE), info = f)
  expect_match(s10, 'compression = "lzw"', fixed = TRUE); expect_match(s10, "dpi = 300", fixed = TRUE)
  expect_match(s10, "Sobol 2007", fixed = TRUE); expect_match(s10, "Saltelli", fixed = TRUE)
  for (f in c("table1_module_screen.csv", "table2_parameters.csv", "table2_footnotes.csv", "table3_heldout.csv", "test_elpd_contrasts.csv",
              "table4_virtual.csv", "tableS_occupancy_sensitivity.csv", "tableS_identifiability.csv", "tableS_sd_sensitivity.csv",
              "tableS_responder_calibrations.csv", "tableS_pk_structure.csv", "tableS_bruun_check.csv", "tableS_sobol.csv",
              "posterior_source_type(", "elpd_contrasts(", "pbo_emax_logit", "consistency check, not a candidate calibration",
              "koff_provenance.csv", "pk_fit_notes.txt", "stage_decision.csv"))
    expect_true(grepl(f, s11, fixed = TRUE), info = f)
})

# ---- R4 fix round 1 (review I4, I9): Fig 3 labels and Fig 6 layers ----
test_that("arm_label gives observed and predicted rows of Fig 3 the same category", {
  arms <- data.frame(trial = "younger2013", arm = c("placebo", "ldn_4.5mg"), day = c(28, 84), stringsAsFactors = FALSE)
  lab <- arm_label(arms$trial, arms$arm, arms$day)
  expect_equal(lab, c("Younger 2013\nPlacebo (day 28)", "Younger 2013\nLDN 4.5 mg (day 84)"))
  pred <- arms[rep(1:2, 4), ]                                     # one row per hypothesis, as script 10 builds them
  expect_setequal(unique(arm_label(pred$trial, pred$arm, pred$day)), lab)
  expect_equal(arm_label("final2023", "ldn", 84), "Final 2023\nLDN (day 84)")   # arm name without a dose: no crash
  s10 <- paste(readLines(file.path("..", "..", "scripts", "10_figures.R")), collapse = "\n")
  expect_match(s10, "day = a$day", fixed = TRUE)
  expect_match(s10, "arms$arm_id <- arm_label(arms$trial, arms$arm, arms$day)", fixed = TRUE)
  expect_match(s10, "pred_rows$arm_id <- arm_label(pred_rows$trial, pred_rows$arm, pred_rows$day)", fixed = TRUE)
  expect_match(s10, "identical(sort(unique(pred_rows$arm_id)), sort(unique(arms$arm_id)))", fixed = TRUE)
})

test_that("Fig 6: the zero line has its own data, the simulated count is computed, and the unused legend frame is gone", {
  s10 <- paste(readLines(file.path("..", "..", "scripts", "10_figures.R")), collapse = "\n")
  expect_match(s10, "geom_hline(data = data.frame(y = 0)", fixed = TRUE)
  expect_false(grepl("5000 simulated", s10, fixed = TRUE))
  expect_false(grepl("ob6_leg", s10, fixed = TRUE))
})
