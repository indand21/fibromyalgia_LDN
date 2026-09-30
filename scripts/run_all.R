# One-command pipeline: Rscript scripts/run_all.R (run from model_v5/).
stopifnot(file.exists("R/utils.R"))
steps <- c("01_fit_pk.R", "02_binding_bottomup.R", "03_test_pet.R", "04_fit_hypotheses.R",
           "05_calibrate_variability.R", "06_test_trials.R", "07_virtual_trials.R",
           "08_test_cytokines.R", "09_sobol.R", "10_figures.R", "11_tables.R")
if (file.exists("output/results/freeze_manifest.csv"))
  stop("A freeze manifest already exists. A clean re-run requires archiving output/results first (see README).")
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
for (s in steps) { cat("\n==>", s, "\n"); sys.source(file.path("scripts", s), envir = new.env()) }
writeLines(capture.output(sessionInfo()), "output/results/sessionInfo.txt")
