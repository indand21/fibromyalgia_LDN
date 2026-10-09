# R3b item 4 and section 6: the run driver's sensitivity-path logic (output directories, no overwrite of primary
# results, the freeze behaviour) and script 07's summary mode. Every pass runs real Rscript processes on small fake
# scripts in a temporary model root; nothing touches the repository's output/.

repo <- normalizePath(file.path("..", ".."), winslash = "/")

wcsv <- function(d, f) { dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); utils::write.csv(d, f, row.names = FALSE); f }
prov <- function(d) cbind(d, source = "fake", location = "fake", verified = "yes")

# A temporary model root: the real R/ (the fake scripts use source_all() and freeze_stage()), fake data, fake scripts,
# and a primary output/results as the primary pass would leave it (stage J adopted, a primary manifest with stage T pins).
fake_root <- function(scripts) {
  root <- normalizePath(tempfile("v5root_"), winslash = "/", mustWork = FALSE); dir.create(root)
  file.copy(file.path(repo, "R"), root, recursive = TRUE)
  wcsv(prov(data.frame(x = 1:3)), file.path(root, "data", "test", "fake_test.csv"))
  wcsv(prov(data.frame(x = 1)), file.path(root, "data", "training", "fake_train.csv"))
  dir.create(file.path(root, "scripts"))
  for (nm in names(scripts)) writeLines(scripts[[nm]], file.path(root, "scripts", nm))
  res <- file.path(root, "output", "results"); dir.create(res, recursive = TRUE)
  writeLines("J", file.path(res, "binding_stage.txt"))
  for (sg in c("A", "A2", "J")) { saveRDS(list(stage = sg), file.path(res, sprintf("exposure_cache_%s.rds", sg)))
                                  saveRDS(list(stage = sg), file.path(res, sprintf("binding_params_%s.rds", sg))) }
  wcsv(data.frame(stage = "A", RMSE = 30), file.path(res, "pet_metrics.csv"))
  wcsv(data.frame(hypothesis = "H0", elpd_primary = -0.5, who = "PRIMARY"), file.path(res, "test_elpd.csv"))
  man <- file.path(res, "freeze_manifest.csv")
  freeze_stage("B", file.path(root, "data", "training", "fake_train.csv"), manifest = man)
  read_test_data(file.path(root, "data", "test", "fake_test.csv"), "B", manifest = man)   # pins stage T
  root
}

HDR <- 'options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()'
fake_steps <- list(
  "04_fit_hypotheses.R" = c(HDR,
    'stage <- readLines("output/results/binding_stage.txt")[1]',
    'cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))',
    'saveRDS(list(stage = stage, cache_stage = cache$stage), "output/results/posterior_H0.rds")',
    'freeze_stage("B", c("output/results/posterior_H0.rds", "output/results/binding_stage.txt"))'),
  "06_test_trials.R" = c(HDR,
    'tt <- read_test_data("data/test/fake_test.csv", "B")',
    'write.csv(data.frame(hypothesis = "H0", elpd_primary = -nrow(tt), who = "SENSITIVITY"), "output/results/test_elpd.csv", row.names = FALSE)'),
  "07_virtual_trials.R" = c(
    'writeLines(Sys.getenv("V5_SCRIPT07_MODE", "full"), "output/results/script07_mode.txt")',
    'write.csv(data.frame(x = 1), "output/results/virtual_summary.csv", row.names = FALSE)',
    'write.csv(data.frame(x = 1), "output/results/virtual_trials.csv", row.names = FALSE)'))
STEPS3 <- names(fake_steps)

test_that("sensitivity outputs go to output/results/sensitivity_stage_<S>/ and only stages A2 and A are allowed", {
  expect_equal(sensitivity_dir("/x", "A2"), file.path("/x", "output", "results", "sensitivity_stage_A2"))
  expect_equal(sensitivity_dir("/x", "A"), file.path("/x", "output", "results", "sensitivity_stage_A"))
  expect_error(sensitivity_dir("/x", "J"), "must be one of")
  expect_error(sensitivity_dir("/x", "B"), "must be one of")
  expect_setequal(SENSITIVITY_STAGES, c("A2", "A"))
})

test_that("the sensitivity passes run 04, 05, 06, 08 and 07 in summary mode, never 09, 10 or 11", {
  expect_setequal(SENSITIVITY_STEPS, c("04_fit_hypotheses.R", "05_calibrate_variability.R", "06_test_trials.R",
                                       "07_virtual_trials.R", "08_test_cytokines.R"))
  expect_false(any(grepl("^(09|10|11)_", SENSITIVITY_STEPS)))
  expect_equal(unname(SENSITIVITY_STEP_ENV[["07_virtual_trials.R"]]), "summary")
  expect_equal(names(SENSITIVITY_STEP_ENV[["07_virtual_trials.R"]]), "V5_SCRIPT07_MODE")
  for (s in c(SENSITIVITY_STEPS, PRIMARY_STEPS)) expect_true(file.exists(file.path(repo, "scripts", s)), info = s)
  expect_equal(PRIMARY_STEPS[1:5], c("01_fit_pk.R", "02_binding_bottomup.R", "03_test_pet.R", "04_fit_hypotheses.R", "05_calibrate_variability.R"))
  expect_equal(which(PRIMARY_STEPS == "05b_younger_floor_study.R"), which(PRIMARY_STEPS == "05_calibrate_variability.R") + 1L)
  expect_setequal(DEFERRED_FAIL_STEPS, c("10_figures.R", "11_tables.R"))   # 05b is a required step (its output is cited)
})

test_that("the driver refuses to start over a freeze manifest or existing sensitivity output, and says how to archive", {
  root <- normalizePath(tempfile("clean_"), winslash = "/", mustWork = FALSE); dir.create(file.path(root, "output", "results"), recursive = TRUE)
  expect_silent(check_clean_start(root))
  writeLines("x", file.path(root, "output", "results", "freeze_manifest.csv"))
  expect_error(check_clean_start(root), "freeze manifest already exists.*archiv")
  unlink(file.path(root, "output", "results", "freeze_manifest.csv"))
  dir.create(file.path(root, "output", "results", "sensitivity_stage_A2"))
  expect_error(check_clean_start(root), "Sensitivity output already exists.*sensitivity_stage_A2.*archiv")
})

test_that("run_all.R guards first, then runs the primary pass, the sensitivity passes and the summary", {
  s <- readLines(file.path(repo, "scripts", "run_all.R")); txt <- paste(s, collapse = "\n")
  first <- function(p) which(grepl(p, s, fixed = TRUE))[1]
  expect_true(is.finite(first("check_clean_start(root)")))
  expect_lt(first("check_clean_start(root)"), first("run_steps("))
  expect_lt(first("run_steps("), first("run_sensitivity_pass("))
  expect_lt(first("run_sensitivity_pass("), first("sensitivity_summary("))
  for (p in c("sensitivity_summary.csv", "run_timing.txt", "sessionInfo.txt", "setdiff(SENSITIVITY_STAGES, stage)"))
    expect_true(grepl(p, txt, fixed = TRUE), info = p)
})

test_that("a work copy holds only the stage's inputs, names the stage, and has no freeze manifest", {
  root <- fake_root(fake_steps); work <- file.path(root, "output", "sensitivity_work", "stage_A2")
  prepare_sensitivity_workspace(root, "A2", work)
  wres <- file.path(work, "output", "results")
  expect_setequal(list.files(wres), c("binding_stage.txt", sensitivity_inputs("A2")))
  expect_equal(readLines(file.path(wres, "binding_stage.txt")), "A2")
  expect_false(file.exists(file.path(wres, "freeze_manifest.csv")))
  expect_true(all(dir.exists(file.path(work, c("R", "scripts", "data")))))
  expect_equal(readRDS(file.path(wres, "exposure_cache_A2.rds"))$stage, "A2")       # never another stage's cache
  expect_error(prepare_sensitivity_workspace(root, "J", work), "unknown sensitivity stage")
  unlink(file.path(root, "output", "results", "binding_params_A2.rds"))
  expect_error(prepare_sensitivity_workspace(root, "A2", work), "missing input file.*binding_params_A2.rds")
  expect_error(prepare_sensitivity_workspace(root, "A", file.path(root, "output", "results", "w")), "must not live inside output/results")
  writeLines("A", file.path(root, "output", "results", "binding_stage.txt"))
  expect_error(prepare_sensitivity_workspace(root, "A", work), "ADOPTED stage")
})

test_that("a sensitivity pass writes only to its own directory and its own manifest; the primary results and manifest are untouched", {
  root <- fake_root(fake_steps); res <- file.path(root, "output", "results")
  before <- hash_dir_files(res); man_rows <- utils::read.csv(file.path(res, "freeze_manifest.csv"))
  tf <- file.path(root, "output", "logs", "run_timing.txt")
  run_sensitivity_pass(root, "A2", steps = STEPS3, timing_file = tf, echo = FALSE)
  dest <- sensitivity_dir(root, "A2")
  # primary outputs and the primary manifest are byte-identical; no file was added to output/results itself
  expect_identical(hash_dir_files(res), before)
  expect_equal(utils::read.csv(file.path(res, "freeze_manifest.csv")), man_rows)
  expect_equal(utils::read.csv(file.path(res, "test_elpd.csv"))$who, "PRIMARY")
  # the pass's own results, from its own stage
  expect_equal(utils::read.csv(file.path(dest, "test_elpd.csv"))$who, "SENSITIVITY")
  expect_equal(readRDS(file.path(dest, "posterior_H0.rds"))$cache_stage, "A2")
  expect_equal(readLines(file.path(dest, "binding_stage.txt")), "A2")
  expect_equal(readLines(file.path(dest, "script07_mode.txt")), "summary")            # the driver's switch reached 07
  # the pass froze stage B in its OWN manifest, pinned the same test file, and kept that manifest under a distinct name
  expect_false(file.exists(file.path(dest, "freeze_manifest.csv")))
  sm <- utils::read.csv(file.path(dest, SENSITIVITY_MANIFEST_NAME))
  expect_true(all(c("B", "T") %in% sm$stage))
  expect_equal(sm$sha256[sm$stage == "T"], man_rows$sha256[man_rows$stage == "T"])
  expect_true(all(grepl("sensitivity_work/stage_A2/", sm$file, fixed = TRUE)))        # every frozen/pinned path is in the work copy
  # inputs and the full virtual-trial file are not copied back; a README is written; the work copy is removed
  expect_false(any(c("exposure_cache_A2.rds", "virtual_trials.csv", "pet_metrics.csv") %in% list.files(dest)))
  expect_true(file.exists(file.path(dest, "README.md")))
  expect_false(dir.exists(file.path(root, "output", "sensitivity_work", "stage_A2")))
  # per-step timings and logs
  tm <- utils::read.delim(tf)
  expect_equal(tm$step, STEPS3); expect_true(all(tm$status == "ok")); expect_equal(unique(tm$pass), "sensitivity_A2")
  expect_true(all(file.exists(file.path(root, "output", "logs", "sensitivity_stage_A2", sub("[.]R$", ".log", STEPS3)))))
  # never overwrites an existing sensitivity directory, and the clean-start guard now refuses
  expect_error(run_sensitivity_pass(root, "A2", steps = STEPS3, timing_file = tf, echo = FALSE), "refusing to overwrite")
  unlink(file.path(res, "freeze_manifest.csv"))
  expect_error(check_clean_start(root), "Sensitivity output already exists")
})

test_that("a pass that writes into the primary output/results is stopped and nothing is collected", {
  bad <- fake_steps
  bad[["06_test_trials.R"]] <- c(bad[["06_test_trials.R"]],
    'write.csv(data.frame(x = 2), "../../results/test_elpd.csv", row.names = FALSE)')   # work copy is output/sensitivity_work/stage_A2
  root <- fake_root(bad)
  expect_error(run_sensitivity_pass(root, "A2", steps = STEPS3, timing_file = tempfile(), echo = FALSE), "changed primary output.*test_elpd.csv")
  expect_false(dir.exists(sensitivity_dir(root, "A2")))
})

test_that("a failing step stops the pass with its log; a figure/table failure is logged, the run continues and is reported", {
  st <- fake_steps; st[["06_test_trials.R"]] <- 'stop("deliberate failure in 06")'
  root <- fake_root(st); tf <- tempfile()
  expect_error(run_sensitivity_pass(root, "A2", steps = STEPS3, timing_file = tf, echo = FALSE), "06_test_trials.R failed.*deliberate failure")
  expect_false(dir.exists(sensitivity_dir(root, "A2")))
  writeLines('stop("figure failure")', file.path(root, "scripts", "10_figures.R"))
  writeLines('writeLines("ran", "after.txt")', file.path(root, "scripts", "11_tables.R"))
  writeLines('stop("05b failure")', file.path(root, "scripts", "05b_younger_floor_study.R"))
  tf2 <- tempfile()
  expect_warning(failed <- run_steps(c("10_figures.R", "11_tables.R"), root, "primary", file.path(root, "logs"), tf2, echo = FALSE),
                 "step 10_figures.R failed")
  expect_equal(failed, "10_figures.R")
  expect_true(file.exists(file.path(root, "after.txt")))
  expect_equal(utils::read.delim(tf2)$status, c("FAILED (deferred; run continued)", "ok"))
  expect_error(run_steps(c("05b_younger_floor_study.R", "11_tables.R"), root, "primary", file.path(root, "logs"), tempfile(), echo = FALSE),
               "05b_younger_floor_study.R failed.*05b failure")                          # 05b is required
  rs <- paste(readLines(file.path(repo, "scripts", "run_all.R")), collapse = " ")
  expect_false(grepl("V5_SKIP_05B", rs, fixed = TRUE))
  expect_match(rs, "these steps FAILED", fixed = TRUE)
})

test_that("test-file pins of a pass must match the primary pass's pins", {
  d <- tempfile(); dir.create(d); f1 <- file.path(d, "a", "t.csv"); f2 <- file.path(d, "b", "t.csv")
  wcsv(prov(data.frame(x = 1)), f1); wcsv(prov(data.frame(x = 2)), f2); tr <- wcsv(prov(data.frame(y = 1)), file.path(d, "tr.csv"))
  m1 <- file.path(d, "m1.csv"); m2 <- file.path(d, "m2.csv")
  freeze_stage("B", tr, manifest = m1); read_test_data(f1, "B", manifest = m1)
  freeze_stage("B", tr, manifest = m2); read_test_data(f2, "B", manifest = m2)
  expect_error(check_test_pins(m2, m1), "differ from the primary pass: t.csv")
  expect_silent(check_test_pins(m1, m1))
  expect_error(check_test_pins(file.path(d, "none.csv"), m1), "wrote no freeze manifest")
})

# ---- script 07 summary mode (the one permitted edit to R3a's scripts) ----
# The REAL script 07 runs in a temporary root whose R/ is the real one plus a stub file sourced last (zzz_), which replaces
# the simulation functions with cheap deterministic ones and records which task sets and power runs were requested.
stub_07 <- c(
  'read_test_data <- function(file, stage, manifest = NULL) readRDS(file.path(v5_root(), "fake_tt.rds"))',
  'default_n_workers <- function() 1L',
  'effect_eligibility <- function(fits, cache, day = 84, threshold_pp = 0.5) data.frame(hypothesis = c("H0", "H1", "H2", "H3"),',
  '  day = day, net_effect_pp = c(0, 5, 0, 5), threshold_pp = threshold_pp, eligible = c(FALSE, TRUE, FALSE, TRUE),',
  '  sobol_role = "stub", power_role = c("type I error check (no drug term)", "power", "not run", "power"))',
  'run_vt_task <- function(task) {',
  '  cat(task$meta$set, "\\n", file = file.path(v5_root(), "vt_calls.txt"), append = TRUE)',
  '  d <- (seq_len(40) - 20) / 20',
  '  v <- data.frame(draw = rep(1:4, each = 10), rep = rep(1:10, 4), diff = d, sd_ldn = 1, sd_pbo = 1)',
  '  for (cut in c(15, 30, 50)) { v[[paste0("resp", cut, "_ldn")]] <- 0.45 - cut / 200 + d / 20; v[[paste0("resp", cut, "_pbo")]] <- 0.35 - cut / 200 + d / 20 }',
  '  attr(v, "pain_end_clip") <- data.frame(arm = c("ldn", "placebo"), low = 0, high = 0, draw = 1)',
  '  attr(v, "baseline_parent") <- list(mu = 6, sigma = 1.5, matched = TRUE)',
  '  v }',
  'power_curves <- function(fits, oms, designs, cache, hyps, ns, n_rep, n_draws, map, roles) {',
  '  cat("power\\n", file = file.path(v5_root(), "power_calls.txt"), append = TRUE)',
  '  do.call(rbind, lapply(c("H0", hyps), function(h) data.frame(hypothesis = h, role = "stub", trial = "final", n_per_arm = ns,',
  '    rej_rate_with_model_unc = 0.05, rej_rate_sampling_only = 0.05, crit_diff_with_model_unc = -1, crit_diff_sampling_only = -1,',
  '    n_reps = n_rep))) }')

fake_tt <- function() {
  base <- data.frame(trial = "final", schedule = c("0:4.5", "0:0"), day = 84, baseline_nrs = c(6.3, 6.2), baseline_sd = c(1.3, 1.6),
                     n = c(50, 50), lo95 = NA_real_, hi95 = NA_real_, location = "", stringsAsFactors = FALSE)
  ch <- cbind(base, arm = c("ldn", "placebo"), endpoint = "change_nrs", estimand = "change from baseline", mean = c(-1.3, -0.9), se = 0.2)
  df <- cbind(base[1, ], arm = "ldn", endpoint = "diff_nrs", estimand = "difference in change", mean = -0.4, se = 0.28)
  rs <- do.call(rbind, lapply(c(15, 30, 50), function(cut) { k <- c(ldn = 20, placebo = 13) - (cut - 15) / 5
    cbind(base, arm = c("ldn", "placebo"), endpoint = paste0("resp", cut), estimand = "responders", mean = k / 50, se = NA_real_) }))
  rs$location <- sprintf("(%d of 50)", as.integer(round(rs$mean * 50)))
  rbind(ch, df, rs)
}

run_07 <- function(mode) {
  root <- normalizePath(tempfile("s07_"), winslash = "/", mustWork = FALSE); dir.create(root)
  file.copy(file.path(repo, "R"), root, recursive = TRUE)
  writeLines(stub_07, file.path(root, "R", "zzz_test_stubs.R"))
  dir.create(file.path(root, "scripts")); file.copy(file.path(repo, "scripts", "07_virtual_trials.R"), file.path(root, "scripts"))
  dir.create(file.path(root, "data", "params"), recursive = TRUE)
  file.copy(file.path(repo, "data", "params", "baseline_bounds.csv"), file.path(root, "data", "params"))
  res <- file.path(root, "output", "results"); dir.create(res, recursive = TRUE)
  writeLines("J", file.path(res, "binding_stage.txt")); saveRDS(list(), file.path(res, "exposure_cache_J.rds"))
  for (h in c("H0", "H1", "H2", "H3")) saveRDS(list(draws = matrix(0, 10, 1, dimnames = list(NULL, "mu_p"))), file.path(res, sprintf("posterior_%s.rds", h)))
  om <- data.frame(hypothesis = c("H0", "H1", "H2", "H3"), omega_p = 2, omega_e = 0.01)
  wcsv(om, file.path(res, "omegas.csv"))
  wcsv(rbind(cbind(calibration = "primary: x", om), cbind(calibration = "SENSITIVITY: y", transform(om, omega_p = 0.7))),
       file.path(res, "omegas_sd_sensitivity.csv"))
  saveRDS(fake_tt(), file.path(root, "fake_tt.rds"))
  r <- run_step("07_virtual_trials.R", root, file.path(root, "s07.log"), env = if (is.null(mode)) character(0) else c(V5_SCRIPT07_MODE = mode))
  list(root = root, res = res, status = r$status, log = readLines(r$log, warn = FALSE),
       calls = if (file.exists(file.path(root, "vt_calls.txt"))) trimws(readLines(file.path(root, "vt_calls.txt"))) else character(0),
       power = file.exists(file.path(root, "power_calls.txt")))
}
POWER_FILES <- c("power_curve.csv", "power_null_check.csv", "power_curve_sd_sensitivity.csv", "virtual_summary_sd_sensitivity.csv")

test_that("script 07 in summary mode writes the summary and responder files and no power or SD-sensitivity file", {
  r <- run_07("summary")
  expect_equal(r$status, 0L, info = paste(utils::tail(r$log, 20), collapse = "\n"))
  for (f in c("virtual_summary.csv", "responder_scores.csv", "virtual_trials_sample.csv", "effect_eligibility.csv"))
    expect_true(file.exists(file.path(r$res, f)), info = f)
  for (f in POWER_FILES) expect_false(file.exists(file.path(r$res, f)), info = f)
  expect_false(r$power)                                         # power_curves() never called
  expect_false("sd_sens" %in% r$calls)                          # no SD-sensitivity virtual trials
  expect_setequal(unique(r$calls), c("main", "resp"))
  rs <- utils::read.csv(file.path(r$res, "responder_scores.csv"))
  expect_equal(nrow(rs), 4 * 3)                                  # 4 hypotheses x 3 cuts for the one fake trial
  expect_true(all(c("ppp_two_sided", "obs_diff") %in% names(rs)))
})

test_that("script 07 without the switch (full mode) still writes the power and SD-sensitivity files; a bad mode stops", {
  r <- run_07(NULL)
  expect_equal(r$status, 0L, info = paste(utils::tail(r$log, 20), collapse = "\n"))
  for (f in c("virtual_summary.csv", "responder_scores.csv", POWER_FILES)) expect_true(file.exists(file.path(r$res, f)), info = f)
  expect_true(r$power); expect_true("sd_sens" %in% r$calls)
  b <- run_07("fast")
  expect_false(b$status == 0L)
  expect_true(any(grepl("V5_SCRIPT07_MODE must be 'full' or 'summary'", b$log, fixed = TRUE)))
})
