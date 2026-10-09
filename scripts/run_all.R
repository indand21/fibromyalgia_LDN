# One-command pipeline: Rscript scripts/run_all.R (run from model_v5/).
#   1. PRIMARY pass, scripts 01-12 with the adopted binding stage (binding_stage.txt, stage J), with 05b after 05 (the
#      curated data/params/baseline_bounds.csv cites its output). If figures (10) or tables (11) fail, the failure is
#      logged, the sensitivity passes still run, and the driver stops at the end naming the failed step.
#   2. SENSITIVITY passes for binding stages A2 and A: scripts 04, 05, 06, 07 (summary mode) and 08 in a work copy, results
#      in output/results/sensitivity_stage_<S>/, with the work copy's own freeze manifest (never the primary one).
#   3. SUMMARY: output/results/sensitivity_summary.csv.
# Per-script wall-clock times go to output/logs/run_timing.txt, per-script logs to output/logs/, and the session
# information to output/results/sessionInfo.txt. The logic lives in R/pipeline.R (tested in test-pipeline.R).
stopifnot(file.exists("R/utils.R"))
options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
root <- v5_root()
Sys.setenv(V5_SCRIPT07_MODE = "full")   # primary pass always runs script 07 in full; a leftover shell value must not skip power
check_clean_start(root)   # refuses to run over an existing freeze manifest or sensitivity output (with archiving advice)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)
timing <- file.path(root, "output", "logs", "run_timing.txt")
if (file.exists(timing)) stop("output/logs/run_timing.txt already exists. ", ARCHIVE_HINT)
t_start <- Sys.time()

# 1. primary pass
failed <- run_steps(PRIMARY_STEPS, root, "primary", file.path(root, "output", "logs"), timing)
stage <- readLines("output/results/binding_stage.txt")[1]
if (stage %in% SENSITIVITY_STAGES) warning("the adopted stage is ", stage, ", which is also a sensitivity stage; that pass is skipped")

# 2. sensitivity passes (the adopted stage is never re-run)
for (sg in setdiff(SENSITIVITY_STAGES, stage)) {
  cat(sprintf("\n==== sensitivity pass: binding stage %s ====\n", sg))
  run_sensitivity_pass(root, sg, timing_file = timing)
}

# 3. summary across binding stages
t0 <- Sys.time()
ss <- sensitivity_summary(root, setdiff(SENSITIVITY_STAGES, stage))
write.csv(ss, "output/results/sensitivity_summary.csv", row.names = FALSE)
append_timing(timing, "summary", "sensitivity_summary.csv", as.numeric(difftime(Sys.time(), t0, units = "secs")), "ok")
print(ss[ss$metric %in% c("mean_steady_state_occupancy_4.5mg", "elpd_primary"), c("stage", "metric", "item", "value", "mcse_monte_carlo")],
      digits = 4, row.names = FALSE)

append_timing(timing, "total", "run_all.R", as.numeric(difftime(Sys.time(), t_start, units = "secs")), "ok")
writeLines(capture.output(sessionInfo()), "output/results/sessionInfo.txt")
if (length(failed)) stop("the pipeline ran to the end, but these steps FAILED and must be rerun after a fix: ",
                         paste(failed, collapse = ", "), " (logs in output/logs/)")
cat(sprintf("\nPipeline finished in %.1f h. Timings: %s\n", as.numeric(difftime(Sys.time(), t_start, units = "hours")), timing))
