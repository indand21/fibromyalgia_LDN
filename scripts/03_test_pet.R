options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
pk <- readRDS("output/results/pk_fit.rds"); bpA <- readRDS("output/results/binding_params_A.rds")
pet <- read_test_data("data/test/pet_occupancy.csv", "A")

NORM_REF <- PET_NORM_REF   # predict_pet() and the reference conventions live in R/pet.R (shared with script 10)
if (!"reference" %in% names(pet)) stop("pet_occupancy.csv lacks the `reference` column")
if (!all(pet$reference %in% c("absolute", NORM_REF))) stop("unknown `reference` value in pet_occupancy.csv")
predict_pet_rows <- function(bp, rows) predict_pet(bp, rows, pk$par)
metrics <- function(obs, pred, stage, rows = "all") data.frame(stage = stage, n = length(obs), RMSE = sqrt(mean((pred - obs)^2)),
                                                  MAE = mean(abs(pred - obs)), n_within_15 = sum(abs(pred - obs) <= 15), rows = rows)

pet$predicted <- predict_pet_rows(bpA, pet); pet$stage <- "A (bottom-up)"
mA <- metrics(pet$occupancy_pct, pet$predicted, "A")
out <- pet; m <- mA; stage <- "A"
dec <- data.frame(stage_A_rmse_heldout = NA_real_, stage_A2_rmse_heldout = NA_real_, stage_A_rmse_all = mA$RMSE, stage_A2_rmse_all = NA_real_,
                  n_heldout = NA_integer_, n_all = nrow(pet), stage_adopted = "A",
                  reason = sprintf("stage A RMSE %.2f <= 15: fallback not triggered, A2 not fitted", mA$RMSE), stringsAsFactors = FALSE)
if (mA$RMSE > 15) {
  msg <- sprintf("bottom-up binding failed PET test (RMSE %.1f). Calibrating Kd scale on Rabiner 2011 only (stage A2).", mA$RMSE)
  cl <- if (file.exists("CHANGELOG.md")) readLines("CHANGELOG.md", warn = FALSE) else character(0)
  if (!any(grepl(msg, cl, fixed = TRUE)))   # append only if an entry with this text (any date) is not already present
    cat(paste0("
## ", Sys.Date(), " - ", msg, "
"), file = "CHANGELOG.md", append = TRUE)
  rab <- pet[grepl("Rabiner", pet$study), ]; rest <- pet[!grepl("Rabiner", pet$study), ]
  sse <- function(ls) { b <- binding_params(bpA$Kd_NTX * exp(ls), bpA$Kd_BN * exp(ls), bpA$koff_NTX, bpA$koff_BN, bpA$keq)
                        sum((predict_pet_rows(b, rab) - rab$occupancy_pct)^2) }
  ls <- stats::optimize(sse, c(-6, 6))$minimum
  bpA2 <- binding_params(bpA$Kd_NTX * exp(ls), bpA$Kd_BN * exp(ls), bpA$koff_NTX, bpA$koff_BN, bpA$keq)
  bpA2$kd_scale <- exp(ls)
  saveRDS(bpA2, "output/results/binding_params_A2.rds")
  saveRDS(exposure_cache(as.numeric(names(readRDS("output/results/exposure_cache_A.rds"))), pk$par, bpA2),
          "output/results/exposure_cache_A2.rds")
  freeze_stage("A2", c("output/results/binding_params_A2.rds", "data/training/pk_summary.csv"))
  rest$predicted <- predict_pet_rows(bpA2, rest); rest$stage <- "A2 (Kd calibrated on Rabiner; held-out)"
  rab$predicted <- predict_pet_rows(bpA2, rab); rab$stage <- "A2 (calibration data)"
  restA <- pet[!grepl("Rabiner", pet$study), ]; rabA <- pet[grepl("Rabiner", pet$study), ]
  out <- rbind(out, rab, rest)
  # Corrected selection rule (run 5): adopt A2 only if it improves RMSE on the held-out (non-calibration) rows relative to
  # stage A on the SAME rows. Matched-row metrics are reported alongside the original A (n = all) and A2 held-out rows.
  dec <- select_binding_stage(rest$occupancy_pct, restA$predicted, rest$predicted,
                              c(rab$occupancy_pct, rest$occupancy_pct), c(rabA$predicted, restA$predicted), c(rab$predicted, rest$predicted))
  m <- rbind(mA,
             metrics(rest$occupancy_pct, rest$predicted, "A2 held-out", "held-out rows only (n = 5)"),
             metrics(restA$occupancy_pct, restA$predicted, "A held-out (matched)", "same held-out rows as A2 held-out"),
             metrics(c(rab$occupancy_pct, rest$occupancy_pct), c(rab$predicted, rest$predicted), "A2 all", "all rows (6 calibration + 5 held-out)"),
             metrics(rab$occupancy_pct, rab$predicted, "A2 calibration", "calibration rows (Rabiner; fitted)"),
             metrics(rabA$occupancy_pct, rabA$predicted, "A calibration", "calibration rows (Rabiner; not fitted under A)"))
  m$rows[1] <- "all rows"
  stage <- dec$stage_adopted
}
write.csv(dec, "output/results/stage_decision.csv", row.names = FALSE)
write.csv(out, "output/results/pet_test.csv", row.names = FALSE)
write.csv(m, "output/results/pet_metrics.csv", row.names = FALSE)
writeLines(stage, "output/results/binding_stage.txt")
# Occupancy of the ACTIVE stage (from its own exposure cache); every figure/table reads this stage, never stage A by default.
occ <- active_occupancy("output/results")
if (file.exists("output/results/exposure_cache_A2.rds")) {   # both stages' occupancy, for comparison and the sensitivity analysis
  both <- do.call(rbind, lapply(c("A", "A2"), function(sg) { o <- steady_state_occupancy(readRDS(sprintf("output/results/exposure_cache_%s.rds", sg)))
    o$stage <- sg; o[o$dose_mg %in% c(1.5, 4.5, 6, 50), c("stage", "dose_mg", "mean_RO", "peak_RO", "trough_RO")] }))
  write.csv(both, "output/results/steady_state_occupancy_both_stages.csv", row.names = FALSE)
  if (stage != "A2") write.csv(steady_state_occupancy(readRDS("output/results/exposure_cache_A2.rds")), "output/results/steady_state_occupancy_A2.csv", row.names = FALSE)
}
write.csv(occ, sprintf("output/results/steady_state_occupancy_%s.csv", stage), row.names = FALSE)
write.csv(occ[occ$dose_mg %in% c(1.5, 4.5, 6, 50), c("stage", "dose_mg", "mean_RO", "peak_RO", "trough_RO")],
          "output/results/steady_state_occupancy_active.csv", row.names = FALSE)
# DOR/KOR descriptive equilibrium occupancy for the active stage (Kd values scaled by the stage's Kd scale; 1 for stage A)
bc <- read_dataset("data/params/binding_constants.csv"); acache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
kd_sc <- if (stage == "A2") readRDS("output/results/binding_params_A2.rds")$kd_scale else 1
ki <- function(an, tg) kd_sc * pool_binding(bc, an, tg, "Ki")$value
dk <- do.call(rbind, lapply(c(4.5, 50), function(d) { pr <- acache[[as.character(d)]]
  do.call(rbind, lapply(c("DOR", "KOR"), function(tg) data.frame(stage = stage, dose_mg = d, target = tg,
    mean_occupancy = equilibrium_occupancy(mean(pr$Cu_NTX), mean(pr$Cu_BN), ki("NTX", tg), ki("BN", tg))))) }))
write.csv(dk, "output/results/dor_kor_occupancy_active.csv", row.names = FALSE)
writeLines(c(sprintf("Stage-A (bottom-up) PET RMSE: %.2f percentage points (n = %d; MAE %.2f; %d within 15 points)",
                     mA$RMSE, mA$n, mA$MAE, mA$n_within_15),
             sprintf("Stage A2 fitted (fallback triggered at RMSE > 15): %s", if (mA$RMSE > 15) "yes" else "no"),
             sprintf("Held-out rows (matched): stage A RMSE %.2f vs stage A2 RMSE %.2f (n = %s); all rows: A %.2f vs A2 %.2f",
                     dec$stage_A_rmse_heldout, dec$stage_A2_rmse_heldout, dec$n_heldout, dec$stage_A_rmse_all, dec$stage_A2_rmse_all),
             sprintf("Decision rule: adopt A2 only if it improves held-out RMSE relative to A on the same rows. %s", dec$reason),
             sprintf("Binding stage carried forward: %s", stage),
             sprintf("Reference convention: %d row(s) 'absolute' (prediction = 100 * RO(t)); %d row(s) '%s' (prediction = 100 * RO(t) / RO(1 h), same dose and regimen); the convention used per row is in the `reference` column of pet_test.csv",
                     sum(pet$reference == "absolute"), sum(pet$reference == NORM_REF), NORM_REF)),
           "output/results/pet_test_note.txt")
print(m)
