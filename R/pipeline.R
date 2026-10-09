# Run-driver logic for scripts/run_all.R (R3b item 4, narrowed by section 6 of the R3b instructions).
#
# PRIMARY pass: scripts 01-12 in the repository with the adopted binding stage (binding_stage.txt, J), and 05b after 05
#   (required: data/params/baseline_bounds.csv cites its output, younger_floor_stability.csv).
# SENSITIVITY passes: for stages A2 and A, scripts 04, 05, 06, 07 (summary mode) and 08 in a WORK COPY of R/, scripts/
#   and data/ whose output/results holds only the stage's exposure cache, its binding parameters, pet_metrics.csv and a
#   binding_stage.txt naming the stage. Every script sets v5.root to its working directory, so the work copy's freeze
#   manifest (and its stage B, test-file pins) is its OWN file: the sensitivity passes can never write to or trip the
#   primary manifest. Results are copied to output/results/sensitivity_stage_<S>/ afterwards; nothing in the primary
#   output/results is written by a pass, and the driver checks that by hashing before and after.
# SUMMARY: output/results/sensitivity_summary.csv (occupancy at 4.5 mg, primary elpd with its Monte Carlo SE, the
#   pairwise contrasts and the responder p-values, per stage).
#
# Each script runs in its own Rscript process with the working directory set to the pass's root, so that no state
# (global functions re-sourced by source_all(), options, the RNG kind) carries from one script into the next.

PRIMARY_STEPS <- c("01_fit_pk.R", "02_binding_bottomup.R", "03_test_pet.R", "04_fit_hypotheses.R",
                   "05_calibrate_variability.R", "05b_younger_floor_study.R", "06_test_trials.R",
                   "07_virtual_trials.R", "08_test_cytokines.R", "09_sobol.R", "10_figures.R", "11_tables.R",
                   "12_derived_quantities.R")   # 12 reads primary outputs only (derived_quantities.csv)
# Figures and tables only present results: if one fails, the failure is logged, the sensitivity passes still run, and
# the driver stops at the end naming it (so a presentation error never costs the long sensitivity passes).
DEFERRED_FAIL_STEPS <- c("10_figures.R", "11_tables.R")
SENSITIVITY_STAGES <- c("A2", "A")
SENSITIVITY_STEPS <- c("04_fit_hypotheses.R", "05_calibrate_variability.R", "06_test_trials.R",
                       "07_virtual_trials.R", "08_test_cytokines.R")
SENSITIVITY_STEP_ENV <- list("07_virtual_trials.R" = c(V5_SCRIPT07_MODE = "summary"))
SENSITIVITY_COPY_DIRS <- c("R", "scripts", "data")
# Files a sensitivity pass needs from the primary pass (written by scripts 02 and 03).
sensitivity_inputs <- function(stage) c(sprintf("exposure_cache_%s.rds", stage), sprintf("binding_params_%s.rds", stage),
                                        "pet_metrics.csv")
# Not copied back: the inputs (identical to the primary's copies) and the full virtual-trial file (~500 MB; the 1%
# sample virtual_trials_sample.csv is kept, as for the primary pass). The work copy's manifest is kept under a name
# that cannot be mistaken for the primary one.
SENSITIVITY_NO_COPY <- c("virtual_trials.csv", "pet_metrics.csv")
SENSITIVITY_MANIFEST_NAME <- "freeze_manifest_sensitivity_pass.csv"

results_dir <- function(root) file.path(root, "output", "results")
primary_manifest <- function(root) file.path(results_dir(root), "freeze_manifest.csv")

sensitivity_dir <- function(root, stage) {
  if (!stage %in% SENSITIVITY_STAGES)
    stop("sensitivity stage must be one of ", paste(SENSITIVITY_STAGES, collapse = ", "), " (got '", stage, "')", call. = FALSE)
  file.path(results_dir(root), paste0("sensitivity_stage_", stage))
}

ARCHIVE_HINT <- paste("A clean re-run requires archiving the previous outputs first: move output/results, output/figures",
                      "and output/logs to ../_archive/<name>_<YYYYMMDD_HHMMSS>/ (nothing is deleted), then run",
                      "Rscript scripts/run_all.R again from model_v5/.")

# Refuse to start over existing results: the primary freeze manifest (as before) and any sensitivity output.
check_clean_start <- function(root, stages = SENSITIVITY_STAGES) {
  if (file.exists(primary_manifest(root)))
    stop("A freeze manifest already exists (output/results/freeze_manifest.csv). ", ARCHIVE_HINT, call. = FALSE)
  old <- vapply(stages, function(s) dir.exists(sensitivity_dir(root, s)), TRUE)
  if (any(old))
    stop("Sensitivity output already exists (", paste(basename(vapply(stages[old], sensitivity_dir, "", root = root)), collapse = ", "),
         "). ", ARCHIVE_HINT, call. = FALSE)
  invisible(TRUE)
}

rscript_path <- function() file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

set_env_vars <- function(vars) {
  if (!length(vars)) return(list())
  old <- Sys.getenv(names(vars), unset = NA, names = TRUE)
  do.call(Sys.setenv, as.list(vars))
  as.list(old)
}
restore_env_vars <- function(old) for (nm in names(old)) if (is.na(old[[nm]])) Sys.unsetenv(nm) else
  do.call(Sys.setenv, stats::setNames(list(old[[nm]]), nm))

# Run one script in its own Rscript process, working directory `wd`, stdout and stderr to `log_file`.
run_step <- function(script, wd, log_file, env = character(0), rscript = rscript_path()) {
  if (!file.exists(file.path(wd, "scripts", script))) stop("script not found: ", file.path(wd, "scripts", script), call. = FALSE)
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log_file <- normalizePath(log_file, winslash = "/", mustWork = FALSE)
  old_wd <- setwd(wd); on.exit(setwd(old_wd), add = TRUE)
  old_env <- set_env_vars(env); on.exit(restore_env_vars(old_env), add = TRUE)
  t0 <- Sys.time()
  status <- system2(rscript, shQuote(file.path("scripts", script)), stdout = log_file, stderr = log_file)
  list(status = as.integer(status), seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), log = log_file)
}

append_timing <- function(file, pass, step, seconds, status, note = "") {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  row <- data.frame(pass = pass, step = step, seconds = round(seconds, 1), status = status, note = note,
                    finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
  utils::write.table(row, file, sep = "\t", row.names = FALSE, quote = FALSE,
                     col.names = !file.exists(file), append = file.exists(file))
}

# Run a list of steps in `wd`, logging each; stop at the first failure, except for a deferred (figure/table) step, whose
# failure is logged and returned so the caller can stop at the end.
run_steps <- function(steps, wd, pass, log_dir, timing_file, step_env = list(), deferred = DEFERRED_FAIL_STEPS,
                      rscript = rscript_path(), echo = TRUE) {
  failed <- character(0)
  for (s in steps) {
    if (echo) cat(sprintf("\n==> [%s] %s\n", pass, s))
    r <- run_step(s, wd, file.path(log_dir, sub("[.]R$", ".log", s)), env = if (is.null(step_env[[s]])) character(0) else step_env[[s]],
                  rscript = rscript)
    ok <- r$status == 0L
    note <- if (length(step_env[[s]])) paste(names(step_env[[s]]), step_env[[s]], sep = "=", collapse = " ") else ""
    append_timing(timing_file, pass, s, r$seconds, if (ok) "ok" else if (s %in% deferred) "FAILED (deferred; run continued)" else "FAILED", note)
    if (echo) cat(sprintf("    %s in %.1f min (log: %s)\n", if (ok) "done" else "FAILED", r$seconds / 60, r$log))
    if (!ok) {
      tail_txt <- utils::tail(readLines(r$log, warn = FALSE), 15)
      if (s %in% deferred) { warning("step ", s, " failed (", pass, "); continuing, the driver stops at the end. Log: ", r$log, call. = FALSE)
                             failed <- c(failed, s); next }
      stop("step ", s, " failed in the ", pass, " pass (exit status ", r$status, "). Log: ", r$log,
           "\nLast lines:\n", paste(tail_txt, collapse = "\n"), call. = FALSE)
    }
  }
  invisible(failed)
}

hash_dir_files <- function(dir) {
  f <- list.files(dir, full.names = TRUE, recursive = FALSE)
  f <- f[!dir.exists(f)]
  stats::setNames(vapply(f, sha_file, ""), basename(f))
}

# Build the work copy for one sensitivity stage. Stops naming any missing input; never uses another stage's files.
prepare_sensitivity_workspace <- function(root, stage, work) {
  if (!stage %in% SENSITIVITY_STAGES) stop("unknown sensitivity stage '", stage, "'", call. = FALSE)
  res <- results_dir(root)
  sf <- file.path(res, "binding_stage.txt")
  if (!file.exists(sf)) stop("sensitivity pass ", stage, ": missing input ", sf, " (run the primary pass first)", call. = FALSE)
  adopted <- readLines(sf, warn = FALSE)[1]
  if (identical(adopted, stage)) stop("sensitivity pass ", stage, ": this is the ADOPTED stage, which the primary pass already ran", call. = FALSE)
  inp <- file.path(res, sensitivity_inputs(stage)); miss <- inp[!file.exists(inp)]
  if (length(miss)) stop("sensitivity pass ", stage, ": missing input file(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (normalizePath(dirname(work), winslash = "/", mustWork = FALSE) == normalizePath(res, winslash = "/", mustWork = FALSE) ||
      startsWith(normalizePath(work, winslash = "/", mustWork = FALSE), paste0(normalizePath(res, winslash = "/", mustWork = FALSE), "/")))
    stop("the sensitivity work copy must not live inside output/results", call. = FALSE)
  if (dir.exists(work)) unlink(work, recursive = TRUE, force = TRUE)
  dir.create(file.path(work, "output", "results"), recursive = TRUE)
  for (d in SENSITIVITY_COPY_DIRS) {
    if (!dir.exists(file.path(root, d))) stop("sensitivity pass ", stage, ": missing directory ", file.path(root, d), call. = FALSE)
    if (!isTRUE(file.copy(file.path(root, d), work, recursive = TRUE)))
      stop("sensitivity pass ", stage, ": failed to copy directory ", d, " into the work copy", call. = FALSE)
  }
  cp <- file.copy(inp, file.path(work, "output", "results"))
  if (!all(cp)) stop("sensitivity pass ", stage, ": failed to copy input(s) ", paste(basename(inp[!cp]), collapse = ", "), call. = FALSE)
  writeLines(stage, file.path(work, "output", "results", "binding_stage.txt"))
  if (file.exists(file.path(work, "output", "results", "freeze_manifest.csv")))
    stop("sensitivity work copy unexpectedly holds a freeze manifest", call. = FALSE)
  writeLines(c(sprintf("Sensitivity work copy for binding stage %s, created %s from %s.", stage, format(Sys.time()), root),
               "Its freeze manifest is this copy's own; the primary manifest is never read or written by this pass."),
             file.path(work, "SENSITIVITY_WORK_COPY.txt"))
  invisible(work)
}

# The test files a pass pinned (stage T) must be byte-identical to the ones the primary pass pinned.
check_test_pins <- function(sens_manifest, primary_manifest_file) {
  if (!file.exists(sens_manifest)) stop("sensitivity pass wrote no freeze manifest (script 04 freezes stage B)", call. = FALSE)
  s <- utils::read.csv(sens_manifest, stringsAsFactors = FALSE); s <- s[s$stage == "T", ]
  if (!file.exists(primary_manifest_file)) return(invisible(s))
  p <- utils::read.csv(primary_manifest_file, stringsAsFactors = FALSE); p <- p[p$stage == "T", ]
  m <- merge(data.frame(f = basename(s$file), s = s$sha256), data.frame(f = basename(p$file), p = p$sha256), by = "f")
  bad <- m$f[m$s != m$p]
  if (length(bad)) stop("sensitivity pass read test file(s) that differ from the primary pass: ", paste(bad, collapse = ", "), call. = FALSE)
  invisible(m)
}

# Copy a finished pass's results to output/results/sensitivity_stage_<S>/ (never into output/results itself).
collect_sensitivity_outputs <- function(work, dest, stage, steps, root) {
  if (normalizePath(dirname(dest), winslash = "/", mustWork = FALSE) != normalizePath(results_dir(root), winslash = "/", mustWork = FALSE) ||
      !grepl("^sensitivity_stage_", basename(dest)))
    stop("sensitivity outputs must go to output/results/sensitivity_stage_<S>/ (got ", dest, ")", call. = FALSE)
  if (dir.exists(dest) && length(list.files(dest, all.files = TRUE, no.. = TRUE)))
    stop("refusing to overwrite existing sensitivity output in ", dest, ". ", ARCHIVE_HINT, call. = FALSE)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  wres <- file.path(work, "output", "results")
  f <- list.files(wres, full.names = TRUE); f <- f[!dir.exists(f)]
  f <- f[!basename(f) %in% c(SENSITIVITY_NO_COPY, sensitivity_inputs(stage)[1])]   # the cache is an input
  man <- f[basename(f) == "freeze_manifest.csv"]; f <- setdiff(f, man)
  ok <- file.copy(f, dest)
  if (!all(ok)) stop("could not copy: ", paste(basename(f[!ok]), collapse = ", "), call. = FALSE)
  if (length(man)) file.copy(man, file.path(dest, SENSITIVITY_MANIFEST_NAME))
  writeLines(c(
    sprintf("# Sensitivity pass: binding stage %s (not the adopted stage)", stage),
    "",
    sprintf("Written by scripts/run_all.R on %s.", format(Sys.time(), "%Y-%m-%d %H:%M")),
    sprintf("Scripts run, in a work copy of R/, scripts/ and data/: %s.", paste(steps, collapse = ", ")),
    "Script 07 ran in SUMMARY mode (V5_SCRIPT07_MODE=summary): virtual-trial summary and responder scores only;",
    "no power curves and no SD-sensitivity calibration runs. Script 09 (Sobol) is not run in sensitivity passes.",
    sprintf("Inputs taken from the primary pass: %s; binding_stage.txt set to %s.", paste(sensitivity_inputs(stage), collapse = ", "), stage),
    sprintf("The work copy's own freeze manifest is kept here as %s; the primary manifest was never touched.", SENSITIVITY_MANIFEST_NAME),
    "virtual_trials.csv (full) is not kept; virtual_trials_sample.csv is."),
    file.path(dest, "README.md"))
  invisible(dest)
}

run_sensitivity_pass <- function(root, stage, steps = SENSITIVITY_STEPS, step_env = SENSITIVITY_STEP_ENV,
                                 work_base = file.path(root, "output", "sensitivity_work"),
                                 log_dir = file.path(root, "output", "logs", paste0("sensitivity_stage_", stage)),
                                 timing_file = file.path(root, "output", "logs", "run_timing.txt"),
                                 keep_work = FALSE, rscript = rscript_path(), echo = TRUE) {
  dest <- sensitivity_dir(root, stage)
  if (dir.exists(dest)) stop("refusing to overwrite existing sensitivity output in ", dest, ". ", ARCHIVE_HINT, call. = FALSE)
  work <- file.path(work_base, paste0("stage_", stage))
  prepare_sensitivity_workspace(root, stage, work)
  before <- hash_dir_files(results_dir(root))
  man_before <- if (file.exists(primary_manifest(root))) sha_file(primary_manifest(root)) else NA_character_
  run_steps(steps, work, paste0("sensitivity_", stage), log_dir, timing_file, step_env = step_env, deferred = character(0),
            rscript = rscript, echo = echo)
  after <- hash_dir_files(results_dir(root))
  changed <- union(setdiff(names(after), names(before)), names(before)[before != after[names(before)] | is.na(after[names(before)])])
  if (length(changed)) stop("sensitivity pass ", stage, " changed primary output(s): ", paste(changed, collapse = ", "), call. = FALSE)
  man_after <- if (file.exists(primary_manifest(root))) sha_file(primary_manifest(root)) else NA_character_
  if (!identical(man_before, man_after)) stop("sensitivity pass ", stage, " changed the primary freeze manifest", call. = FALSE)
  check_test_pins(file.path(work, "output", "results", "freeze_manifest.csv"), primary_manifest(root))
  collect_sensitivity_outputs(work, dest, stage, steps, root)
  if (!keep_work) unlink(work, recursive = TRUE, force = TRUE)
  invisible(dest)
}

# ---- sensitivity_summary.csv ----
STAGE_ROLE <- c(J = "PRIMARY: adopted; Kd and koff scales calibrated jointly on all 11 PET rows (post hoc)",
                A2 = "SENSITIVITY: Kd scale calibrated on the Rabiner rows only (koff in vitro)",
                A = "SENSITIVITY: bottom-up in vitro binding (the held-out PET test of this stage failed)")

sensitivity_summary <- function(root, stages = SENSITIVITY_STAGES) {
  res <- results_dir(root)
  rd <- function(f) { if (!file.exists(f)) stop("sensitivity summary: missing input ", f, call. = FALSE)
    utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }
  adopted <- { f <- file.path(res, "binding_stage.txt"); if (!file.exists(f)) stop("sensitivity summary: missing input ", f, call. = FALSE)
               readLines(f, warn = FALSE)[1] }
  occ <- rd(file.path(res, "steady_state_occupancy_by_stage.csv"))
  unc <- if (file.exists(file.path(res, sprintf("occupancy_uncertainty_%s.csv", adopted))))
    rd(file.path(res, sprintf("occupancy_uncertainty_%s.csv", adopted))) else NULL
  row <- function(stage, metric, item, value, se_obs = NA_real_, mcse = NA_real_, lo = NA_real_, hi = NA_real_, note = "")
    data.frame(stage = stage, stage_role = if (stage %in% names(STAGE_ROLE)) STAGE_ROLE[[stage]] else "", metric = metric,
               item = item, value = value, se_observation_level = se_obs, mcse_monte_carlo = mcse, lo = lo, hi = hi,
               note = note, stringsAsFactors = FALSE)
  out <- list()
  for (sg in c(adopted, stages)) {
    dir <- if (sg == adopted) res else sensitivity_dir(root, sg)
    if (!dir.exists(dir)) stop("sensitivity summary: missing directory ", dir, call. = FALSE)
    o <- occ[occ$stage == sg & abs(occ$dose_mg - 4.5) < 1e-9, ]
    if (nrow(o) != 1) stop("sensitivity summary: steady_state_occupancy_by_stage.csv has no single 4.5 mg row for stage ", sg, call. = FALSE)
    b <- if (sg == adopted && !is.null(unc)) unc[unc$uncertainty_source == "binding_and_pk" & abs(unc$dose_mg - 4.5) < 1e-9, ] else NULL
    out[[length(out) + 1]] <- row(sg, "mean_steady_state_occupancy_4.5mg", "4.5 mg once daily", o$mean_RO,
      lo = if (length(b) && nrow(b) == 1) b$mean_RO_lo95 else NA_real_, hi = if (length(b) && nrow(b) == 1) b$mean_RO_hi95 else NA_real_,
      note = if (length(b) && nrow(b) == 1) "lo/hi: binding + PK Laplace 95% band" else "point value (no band computed for this stage)")
    el <- rd(file.path(dir, "test_elpd.csv"))
    out[[length(out) + 1]] <- row(sg, "elpd_primary", el$hypothesis, el$elpd_primary, mcse = el$mcse_elpd_primary,
                                  note = el$primary_endpoints)
    tp <- rd(file.path(dir, "test_predictions.csv"))
    cmf <- file.path(dir, "test_elpd_contrast_mcse.csv")
    ct <- elpd_contrasts(tp, if (file.exists(cmf)) rd(cmf) else NULL)
    out[[length(out) + 1]] <- row(sg, "elpd_contrast_a_minus_b", paste(ct$hypothesis_a, "-", ct$hypothesis_b), ct$elpd_diff_a_minus_b,
                                  se_obs = ct$se_diff_observation_level, mcse = ct$mcse_diff_monte_carlo,
                                  note = ifelse(ct$se_collapsed_flag, ct$flag_note, sprintf("diff/SE %.2f (n = %d)", ct$ratio_diff_over_se, ct$n_obs)))
    rs <- rd(file.path(dir, "responder_scores.csv"))
    out[[length(out) + 1]] <- row(sg, "responder_ppp_two_sided", sprintf("%s %s %d%% cut", rs$hypothesis, rs$trial, rs$cut_pct),
                                  rs$ppp_two_sided, note = sprintf("observed diff %+.3f; predicted median %+.3f", rs$obs_diff, rs$pred_diff_median))
    vs <- rd(file.path(dir, "virtual_summary.csv")); vs <- vs[vs$omega_mult == 1, ]
    out[[length(out) + 1]] <- row(sg, "virtual_diff_ppp_two_sided", sprintf("%s %s", vs$hypothesis, vs$trial), vs$ppp_two_sided,
                                  lo = vs$diff_lo90, hi = vs$diff_hi90,
                                  note = sprintf("observed %+.3f; simulated median %+.3f (lo/hi: 90%% interval)", vs$obs_diff, vs$diff_median))
  }
  do.call(rbind, out)
}
