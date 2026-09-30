options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
pk <- readRDS("output/results/pk_fit.rds")
bc <- read_dataset("data/params/binding_constants.csv")
# Select on analyte (NTX/BN); pool human rows by geometric mean (see pool_binding in R/data_io.R).
pick <- function(an, q, tgt = "MOR") pool_binding(bc, an, tgt, q)$value
kr <- bc[bc$quantity == "koff", ]
if (nrow(kr) && any(is.na(kr$unit) | kr$unit != "1/h"))
  stop("binding_constants.csv koff row(s) must have unit '1/h' (the model integrates in hours); found: ",
       paste(unique(kr$unit), collapse = ", "))
koff_ntx <- pick("NTX", "koff"); koff_bn <- pick("BN", "koff")
bp <- binding_params(Kd_NTX = pick("NTX", "Ki"), Kd_BN = pick("BN", "Ki"),
                     koff_NTX = ifelse(is.na(koff_ntx), 2, koff_ntx),
                     koff_BN = ifelse(is.na(koff_bn), 2, koff_bn), keq = 2)
stopifnot(is.finite(bp$Kd_NTX), is.finite(bp$Kd_BN))
bp$koff_source <- ifelse(is.na(koff_ntx), "assumed: rapid equilibrium (no human koff found)", "binding_constants.csv")
ASSUMED_KOFF_TXT <- "assumed: koff = 2 per h (no curated koff row for this analyte; rapid-equilibrium assumption)"
bp$koff_NTX_source <- if (is.na(koff_ntx)) ASSUMED_KOFF_TXT else pool_binding(bc, "NTX", "MOR", "koff")$source
bp$koff_BN_source <- if (is.na(koff_bn)) ASSUMED_KOFF_TXT else pool_binding(bc, "BN", "MOR", "koff")$source
write.csv(data.frame(analyte = c("NTX", "BN"), koff_per_h = c(bp$koff_NTX, bp$koff_BN), koff_assumed = c(is.na(koff_ntx), is.na(koff_bn)),
                     source = c(bp$koff_NTX_source, bp$koff_BN_source)), "output/results/koff_provenance.csv", row.names = FALSE)
bp$Kd_NTX_source <- pool_binding(bc, "NTX", "MOR", "Ki")$source
bp$Kd_BN_source <- pool_binding(bc, "BN", "MOR", "Ki")$source
doses <- c(0, 0.75, 1.5, 2, 3, 3.88, 4.5, 5, 5.4, 6, 10, 15, 20, 50)
cache <- exposure_cache(doses, pk$par, bp)
saveRDS(bp, "output/results/binding_params_A.rds"); saveRDS(cache, "output/results/exposure_cache_A.rds")
ss <- steady_state_occupancy(cache)
write.csv(ss, "output/results/steady_state_occupancy_A.csv", row.names = FALSE)
# DOR/KOR: descriptive equilibrium occupancy at the steady-state mean unbound concentrations (spec §4.2)
ki <- function(an, tgt) pick(an, "Ki", tgt)
dk <- do.call(rbind, lapply(c(4.5, 50), function(d) { pr <- cache[[as.character(d)]]
  do.call(rbind, lapply(c("DOR", "KOR"), function(tg) data.frame(dose_mg = d, target = tg,
    mean_occupancy = equilibrium_occupancy(mean(pr$Cu_NTX), mean(pr$Cu_BN), ki("NTX", tg), ki("BN", tg))))) }))
write.csv(dk, "output/results/dor_kor_occupancy_A.csv", row.names = FALSE)
freeze_stage("A", c("data/training/pk_summary.csv", "data/params/binding_constants.csv",
                    "output/results/pk_fit.rds", "output/results/binding_params_A.rds"))
print(ss)
