options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
# Training data only: Younger 2013 arms with a reported SD of % reduction (sd_pct).
stage <- readLines("output/results/binding_stage.txt")[1]
cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
arms <- read_dataset("data/training/efficacy_training.csv")
priors <- read_dataset("data/priors/placebo_priors.csv")
sd_res <- priors$value[priors$parameter == "sd_res_nrs"]
if (length(sd_res) != 1 || !is.finite(sd_res)) stop("priors need exactly one finite sd_res_nrs row (cited NRS test-retest SD)")
bounds <- read_baseline_bounds("data/params/baseline_bounds.csv")   # per-trial enrolment bounds (R2 fix round 1)
om_row <- function(h, om) {
  data.frame(hypothesis = h, omega_p = om$omega_p, omega_e = om$omega_e, sd_res = om$sd_res,
             achieved_sd_pbo = om$achieved_sd_pbo, target_sd_pbo = om$target_sd_pbo,
             achieved_sd_ldn = om$achieved_sd_ldn, target_sd_ldn = om$target_sd_ldn,
             at_bound = om$at_bound, scale_used = om$scale_used, base_sd_source = om$base_sd_source,
             baseline_bounds = om$baseline_bounds, baseline_parent_mu = om$baseline_parent_mu,
             baseline_parent_sigma = om$baseline_parent_sigma, baseline_parent_matched = om$baseline_parent_matched,
             baseline_sampling = om$baseline_sampling,
             pain_end_clip_low_pbo = om$pain_end_clip_low_pbo, pain_end_clip_high_pbo = om$pain_end_clip_high_pbo,
             pain_end_clip_low_ldn = om$pain_end_clip_low_ldn, pain_end_clip_high_ldn = om$pain_end_clip_high_ldn,
             noise_sd_pbo = om$noise_sd_pbo)
}
fits <- lapply(stats::setNames(nm = c("H0", "H1", "H2", "H3")), function(h) readRDS(sprintf("output/results/posterior_%s.rds", h)))
omegas <- do.call(rbind, lapply(names(fits), function(h) om_row(h, calibrate_omegas(fits[[h]], arms, cache, sd_res = sd_res, bounds = bounds))))
write.csv(omegas, "output/results/omegas.csv", row.names = FALSE)
print(omegas)

# ---- SD sensitivity (R3a item 4, review I10) ----
# calibrate_omegas() moment-matches the SD of % reduction derived from Younger 2013's reported 95% CI (27.9 placebo,
# 24.0 LDN). Bested et al. 2023 quotes the same statistics with SDs of 14.6 and 12.5. Re-run the calibration with those
# values, in memory only (the curated file and the PRIMARY calibration above are unchanged), same seed and cohort.
SD_SENS <- c(placebo = 14.6, ldn = 12.5)
arms_sens <- override_sd_pct(arms, pbo = SD_SENS[["placebo"]], ldn = SD_SENS[["ldn"]])
sens <- do.call(rbind, lapply(names(fits), function(h)
  om_row(h, suppressWarnings(calibrate_omegas(fits[[h]], arms_sens, cache, sd_res = sd_res, bounds = bounds)))))
lab_p <- "primary: SD of % reduction derived from the reported 95% CI (27.9 placebo, 24.0 LDN)"
lab_s <- sprintf("SENSITIVITY: SDs quoted by Bested et al. 2023 (%.1f placebo, %.1f LDN)", SD_SENS[["placebo"]], SD_SENS[["ldn"]])
both <- rbind(cbind(calibration = lab_p, omegas), cbind(calibration = lab_s, sens))
key <- c("hypothesis", "omega_p", "omega_e", "achieved_sd_pbo", "target_sd_pbo", "achieved_sd_ldn", "target_sd_ldn", "at_bound")
cmp <- merge(omegas[, key], sens[, key], by = "hypothesis", suffixes = c("_primary", "_sens"))
cmp$omega_p_ratio_sens_over_primary <- cmp$omega_p_sens / cmp$omega_p_primary
cmp$omega_e_ratio_sens_over_primary <- cmp$omega_e_sens / cmp$omega_e_primary
both <- merge(both, cmp[, c("hypothesis", "omega_p_ratio_sens_over_primary", "omega_e_ratio_sens_over_primary")], by = "hypothesis")
both$note <- paste("SENSITIVITY rows: calibration repeated with the Bested et al. 2023 SDs in memory; omegas.csv (primary)",
                   "is unchanged. A target below the noise-only SD of % reduction (noise_sd_pbo) cannot be reached and",
                   "pins omega at its lower bound (at_bound = TRUE).")
both <- both[order(both$hypothesis, both$calibration != lab_p), ]
write.csv(both, "output/results/omegas_sd_sensitivity.csv", row.names = FALSE)
print(cmp, digits = 4)
