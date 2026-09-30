options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
# Training data only: Younger 2013 arms with a reported SD of % reduction (sd_pct).
stage <- readLines("output/results/binding_stage.txt")[1]
cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
arms <- read_dataset("data/training/efficacy_training.csv")
priors <- read_dataset("data/priors/placebo_priors.csv")
sd_res <- priors$value[priors$parameter == "sd_res_nrs"]
if (length(sd_res) != 1 || !is.finite(sd_res)) stop("priors need exactly one finite sd_res_nrs row (cited NRS test-retest SD)")
rows <- lapply(c("H0", "H1", "H2", "H3"), function(h) {
  fit <- readRDS(sprintf("output/results/posterior_%s.rds", h))
  om <- calibrate_omegas(fit, arms, cache, sd_res = sd_res)
  data.frame(hypothesis = h, omega_p = om$omega_p, omega_e = om$omega_e, sd_res = om$sd_res,
             achieved_sd_pbo = om$achieved_sd_pbo, target_sd_pbo = om$target_sd_pbo,
             achieved_sd_ldn = om$achieved_sd_ldn, target_sd_ldn = om$target_sd_ldn,
             at_bound = om$at_bound, scale_used = om$scale_used, base_sd_source = om$base_sd_source)
})
omegas <- do.call(rbind, rows)
write.csv(omegas, "output/results/omegas.csv", row.names = FALSE)
print(omegas)
