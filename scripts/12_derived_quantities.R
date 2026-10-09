# Derived quantities the manuscript cites that no earlier script writes (first computed by hand in the R4 run
# verification, .superpowers/sdd/2026-09-29-ldn-qsp-v5/R4-run-verification-report.md). Reads ONLY primary outputs in
# output/results (never a held-out test-data file, never the freeze manifest) and writes one file:
#   output/results/derived_quantities.csv   quantity, value, lo, hi, unit, definition, source_files
# Rows:
#   (a) PK Laplace draws of the binding + PK occupancy band (script 03: laplace_draws(pk$theta, PK covariance, n,
#       seed = UNC_SEED + 1), with n and UNC_SEED read back from occupancy_uncertainty_J.csv): fraction with fp0 < 0.95
#       (reversed dose dependence of first pass) and the oral bioavailability at 4.5 mg over those draws;
#   (b) minimum tail ESS (posterior::ess_tail) per hypothesis H0-H3 over the sampled parameters, and which parameter;
#   (c) posterior distribution of the day-84 net effect (net_effect_par, the effect_eligibility.csv definition) over
#       400 evenly thinned draws, for H1 and H3 (deterministic subset; no seed involved);
#   (d) mean steady-state unbound plasma naltrexone (nM) over one dosing interval at 1.5, 4.5 and 6 mg once daily, from
#       the adopted stage's exposure cache (unbound NTX depends on the PK only, so it is the same in every stage cache).
#   (e) absolute MOR occupancy (%) of the adopted binding stage 1, 2, 3 and 4 h after a single 50 mg dose under the
#       adopted PK (predict_pet), i.e. the absolute occupancy at the Lee 1988 1 h normalisation reference.
options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
res <- "output/results"
t0 <- Sys.time()

need <- function(...) {
  f <- file.path(res, c(...)); miss <- f[!file.exists(f)]
  if (length(miss)) stop("Derived-quantities script cannot run; missing input file(s): ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(f)
}
FP0_THRESHOLD <- 0.95; F_DOSE <- 4.5            # (a)
NET_DAY <- 84; NET_HYPS <- c("H1", "H3"); NET_N_DRAWS <- 400   # (c)
CU_DOSES <- c(1.5, 4.5, 6)                      # (d)
LEE_DOSE <- 50; LEE_TIMES <- 1:4                # (e)
HYPS <- paste0("H", 0:3)

need("pk_fit.rds", "occupancy_uncertainty_J.csv", "binding_stage.txt", sprintf("posterior_%s.rds", HYPS))
stage <- readLines(file.path(res, "binding_stage.txt"), warn = FALSE)[1]
cache_file <- sprintf("exposure_cache_%s.rds", stage); need(cache_file)
pk <- readRDS(file.path(res, "pk_fit.rds"))
cache <- readRDS(file.path(res, cache_file))
fits <- lapply(stats::setNames(nm = HYPS), function(h) readRDS(file.path(res, sprintf("posterior_%s.rds", h))))
rows <- list()

# ---- (a) PK Laplace draws of the binding + PK band ----
unc <- utils::read.csv(file.path(res, "occupancy_uncertainty_J.csv"), stringsAsFactors = FALSE)
bpk <- unc[!is.na(unc$uncertainty_source) & unc$uncertainty_source == "binding_and_pk", , drop = FALSE]
src_a <- c("pk_fit.rds", "occupancy_uncertainty_J.csv")
if (nrow(bpk) && isTRUE(pk$hessian_summary$pd)) {
  unc_seed <- unique(bpk$seed); n_draws <- unique(bpk$n_draws_requested)
  if (length(unc_seed) != 1 || length(n_draws) != 1) stop("occupancy_uncertainty_J.csv: binding_and_pk rows disagree on seed or n_draws_requested")
  pkd <- pk_laplace_draw_summary(pk, n_draws, seed = unc_seed + 1, dose_mg = F_DOSE)
  how <- sprintf("the %d PK Laplace draws of script 03's binding + PK band (laplace_draws(pk$theta, pk$hessian_summary$cov, %d, seed = %d + 1 = %d), back-transformed with pk_from_theta)",
                 n_draws, n_draws, unc_seed, unc_seed + 1)
  qF <- stats::quantile(pkd$F, c(0.5, 0.025, 0.975), names = FALSE)
  rows <- c(rows, list(
    dq_row("pk_draws_fraction_fp0_below_0.95", mean(pkd$fp0 < FP0_THRESHOLD), unit = "fraction",
           definition = paste0("Fraction of ", how, " with fp0 < ", FP0_THRESHOLD, " (first pass at low dose below the 100 mg value: reversed dose dependence)"),
           source_files = src_a),
    dq_row("pk_draws_bioavailability_4.5mg", qF[1], qF[2], qF[3], unit = "fraction",
           definition = paste0("Oral bioavailability F at 4.5 mg (pk_bioavailability) over ", how, ": median, lo/hi = 2.5th/97.5th percentiles"),
           source_files = src_a)))
} else {
  rows <- c(rows, list(
    dq_row("pk_draws_fraction_fp0_below_0.95", NA, unit = "fraction", definition = "NOT COMPUTED: script 03 drew no PK Laplace draws (PK Hessian not positive definite)", source_files = src_a),
    dq_row("pk_draws_bioavailability_4.5mg", NA, unit = "fraction", definition = "NOT COMPUTED: script 03 drew no PK Laplace draws (PK Hessian not positive definite)", source_files = src_a)))
}

# ---- (b) minimum tail ESS per hypothesis ----
for (h in HYPS) {
  e <- min_ess_tail(fits[[h]])
  rows <- c(rows, list(dq_row(sprintf("min_ess_tail_%s", h), e$value, unit = "draws",
    definition = sprintf("Minimum tail ESS (posterior::ess_tail, %d chains x %d draws) over the %d sampled parameters of %s; attained by parameter %s. Not part of the convergence gate (rhat <= 1.05, bulk ESS >= 400)",
                         length(unique(fits[[h]]$chain)), nrow(fits[[h]]$draws) / length(unique(fits[[h]]$chain)), ncol(fits[[h]]$draws), h, e$parameter),
    source_files = sprintf("posterior_%s.rds", h))))
}

# ---- (c) posterior distribution of the day-84 net effect ----
for (h in NET_HYPS) {
  ne <- net_effect_draws(fits[[h]], cache, NET_DAY, NET_N_DRAWS)
  q <- stats::quantile(ne, c(0.5, 0.05, 0.95), names = FALSE)
  rows <- c(rows, list(dq_row(sprintf("net_effect_d%d_posterior_%s", NET_DAY, h), q[1], q[2], q[3], unit = "pp",
    definition = sprintf("Net effect at day %d (placebo minus LDN 4.5 mg relative pain x 100; net_effect_par, the effect_eligibility.csv definition, which evaluates it at the posterior-median parameters) over %d evenly thinned posterior draws of %s (rows round(seq(1, %d, length.out = %d)), chains in stored order): median, lo/hi = 5th/95th percentiles",
                         NET_DAY, length(ne), h, nrow(fits[[h]]$draws), NET_N_DRAWS),
    source_files = c(sprintf("posterior_%s.rds", h), cache_file))))
}

# ---- (d) mean steady-state unbound plasma naltrexone ----
for (d in CU_DOSES) {
  key <- as.character(d)
  if (is.null(cache[[key]])) stop(cache_file, " has no profile for ", d, " mg")
  rows <- c(rows, list(dq_row(sprintf("mean_unbound_ntx_%smg", key), mean_unbound_ntx(cache[[key]]), unit = "nM",
    definition = sprintf("Mean unbound plasma naltrexone (Cu_NTX) over one steady-state dosing interval [0, 24) h at %s mg once daily under the adopted PK (exposure cache of the adopted binding stage %s; Cu_NTX does not depend on the binding stage)", key, stage),
    source_files = c(cache_file, "binding_stage.txt"))))
}

# ---- (e) absolute MOR occupancy after a single 50 mg dose (Lee 1988 1 h normalisation reference) ----
bp_file <- sprintf("binding_params_%s.rds", stage); need(bp_file)
bp <- readRDS(file.path(res, bp_file))
lee_rows <- data.frame(study = "Lee 1988 1 h reference", dose_mg = LEE_DOSE, regimen = "single", time_h = LEE_TIMES,
                       reference = "absolute", n_doses = 1, stringsAsFactors = FALSE)
lee_occ <- predict_pet(bp, lee_rows, pk$par)
for (j in seq_along(LEE_TIMES)) {
  rows <- c(rows, list(dq_row(sprintf("mor_occupancy_abs_%gmg_single_%dh", LEE_DOSE, LEE_TIMES[j]), lee_occ[j], unit = "%",
    definition = sprintf("Absolute MOR occupancy (stage %s binding, adopted PK; predict_pet with reference = absolute) %d h after a single %g mg dose: the absolute occupancy at Lee 1988's 1 h normalisation reference (Lee rows are reported relative to the 1 h scan)",
                         stage, LEE_TIMES[j], LEE_DOSE),
    source_files = c("pk_fit.rds", bp_file, "binding_stage.txt"))))
}

out <- do.call(rbind, rows)
stopifnot(identical(names(out), DERIVED_COLS))
utils::write.csv(out, file.path(res, "derived_quantities.csv"), row.names = FALSE)
print(out[, c("quantity", "value", "lo", "hi", "unit")], digits = 4, row.names = FALSE)
cat(sprintf("12_derived_quantities: %d rows written to %s in %.1f s\n", nrow(out), file.path(res, "derived_quantities.csv"),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
