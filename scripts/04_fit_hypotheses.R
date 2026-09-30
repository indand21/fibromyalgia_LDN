options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
# Task 7 (PET test) must have run first: it decides the binding stage.
for (req in c("output/results/binding_stage.txt", "output/results/pet_metrics.csv"))
  if (!file.exists(req)) stop("missing ", req, "; run Task 7 (PET test) before fitting")
stage <- readLines("output/results/binding_stage.txt")[1]
cache_file <- sprintf("output/results/exposure_cache_%s.rds", stage)
cache <- readRDS(cache_file)
arms_all <- read_dataset("data/training/efficacy_training.csv")
# Arms without a reported dispersion cannot enter a likelihood that needs an SE: exclude (never impute).
use <- is.finite(arms_all$se)
used_tab <- data.frame(trial = arms_all$trial, arm = arms_all$arm, day = arms_all$day, mean = arms_all$mean, se = arms_all$se,
                       status = ifelse(use, "included", "excluded"),
                       reason = ifelse(use, "", "no dispersion reported in source"))
print(used_tab)
write.csv(used_tab, "output/results/training_arms_used.csv", row.names = FALSE)
arms <- arms_all[use, , drop = FALSE]
if (nrow(arms) == 0) stop("no training arms with finite se remain")
priors <- read_dataset("data/priors/placebo_priors.csv")
bc <- read_dataset("data/params/binding_constants.csv")
# Select on analyte (NTX); pool by geometric mean, preferring human rows (fallback to non-human is
# stated in the recorded source text). Every contributing source string is recorded.
pool_tlr4 <- pool_binding(bc, "NTX", "TLR4", c("IC50", "Kd"))
pool_ogfr <- pool_binding(bc, "NTX", "OGFr", c("Kd", "Ki"))
if (pool_tlr4$n == 0) stop("no TLR4 IC50/Kd rows curated for NTX")
fixed <- list(IC50_TLR4 = pool_tlr4$value,
              Kd_OGFr = if (pool_ogfr$n > 0) pool_ogfr$value else NA)
fixed_sources <- list(IC50_TLR4 = pool_tlr4$source, Kd_OGFr = pool_ogfr$source)
kd_fixed <- is.finite(fixed$Kd_OGFr)
trials <- unique(arms$trial)
fits <- list(); loo_rows <- list()
for (h in c("H0", "H1", "H2", "H3")) {
  lp <- function(p) { pr <- log_prior(p, h, priors, trials); if (!is.finite(pr)) return(-Inf)
                      pr + log_lik_arms(p, h, arms, cache, fixed) }
  init <- init_par(h, trials)
  if (h == "H3" && kd_fixed) init <- init[setdiff(names(init), "log10_Kd_OGFr")]
  f <- run_chains(lp, init, n_chains = 4, n_iter = 60000, n_burn = 30000)
  f$hyp <- h; f$fixed <- fixed; f$fixed_sources <- fixed_sources; f$trials <- trials; f$stage <- stage
  keep <- round(seq(1, nrow(f$draws), length.out = 2000))
  ll <- t(vapply(keep, function(i) log_lik_arms(f$draws[i, ], h, arms, cache, fixed, pointwise = TRUE), numeric(nrow(arms))))
  reff <- loo::relative_eff(exp(ll), chain_id = f$chain[keep])
  lo <- loo::loo(ll, r_eff = reff)
  k <- lo$diagnostics$pareto_k
  loo_rows[[h]] <- data.frame(hypothesis = h, elpd_loo = lo$estimates["elpd_loo", 1], se = lo$estimates["elpd_loo", 2],
                              p_loo = lo$estimates["p_loo", 1], max_pareto_k = max(k),
                              n_k_gt_0.7 = sum(k > 0.7), ranking_valid = sum(k > 0.7) == 0)
  saveRDS(f, sprintf("output/results/posterior_%s.rds", h)); fits[[h]] <- f
}
write.csv(do.call(rbind, loo_rows), "output/results/loo_training.csv", row.names = FALSE)
fit_tab <- do.call(rbind, lapply(fits, function(f) cbind(hypothesis = f$hyp, f$summary)))
write.csv(fit_tab, "output/results/training_fit.csv", row.names = FALSE)

# Per-hypothesis, per-chain diagnostics (written before the gate so they exist even if the gate fails)
diag_tab <- do.call(rbind, lapply(fits, function(f) data.frame(hypothesis = f$hyp, chain = seq_along(f$accept),
  acceptance_rate = f$accept, min_ess_bulk = min(f$summary$ess_bulk, na.rm = TRUE), max_rhat = max(f$summary$rhat, na.rm = TRUE))))
write.csv(diag_tab, "output/results/mcmc_diagnostics.csv", row.names = FALSE)

# Convergence gate: no freeze unless every hypothesis converged
bad <- fit_tab[is.na(fit_tab$rhat) | is.na(fit_tab$ess_bulk) | fit_tab$rhat > 1.05 | fit_tab$ess_bulk < 400, ]
if (nrow(bad)) stop("convergence failure (rhat > 1.05, NA, or ess_bulk < 400); stage B NOT frozen. Failing: ",
                    paste(sprintf("%s:%s", bad$hypothesis, bad$variable), collapse = ", "),
                    ". Increase n_iter in run_chains and rerun.")

# Observed ED50/ED95 come from the estimation-side parameter file (not the held-out PET ED50 file)
dr <- read_dataset("data/params/dose_response.csv")
obs_dr <- function(q) { v <- dr$value[dr$quantity == q]
  if (length(v) != 1 || !is.finite(v)) stop("data/params/dose_response.csv needs exactly one finite ", q, " row (found ", length(v), ")")
  v }
obs_ED50 <- obs_dr("ED50"); obs_ED95 <- obs_dr("ED95")
# Bruun-Plesner consistency check: effect vs dose at day 365 (all grid doses are in the cache)
doses <- c(0, 0.75, 1.5, 2, 3, 3.88, 4.5, 5, 5.4, 6)
grid_rows <- list(); bp_rows <- list()
for (h in c("H1", "H2", "H3")) {
  f <- fits[[h]]; med <- apply(f$draws, 2, stats::median); th <- theta_from(med, h, fixed)
  eff <- vapply(doses, function(d) 1 - tail(simulate_arm(h, th, parse_schedule(paste0("0:", d)), 365, cache, 0)$rel, 1), 0)
  grid_rows[[h]] <- data.frame(hypothesis = h, dose_mg = doses, effect = eff)
  mx <- max(eff); ed50 <- NA_real_
  if (mx > 0) {
    thr <- 0.5 * mx
    for (i in seq_along(doses)[-1]) if (eff[i] >= thr && eff[i - 1] < thr) {
      ed50 <- doses[i - 1] + (thr - eff[i - 1]) / (eff[i] - eff[i - 1]) * (doses[i] - doses[i - 1]); break }
  }
  bp_rows[[h]] <- data.frame(hypothesis = h, peak_dose_mg = doses[which.max(eff)], model_ED50_mg = ed50,
                             observed_ED50_mg = obs_ED50, observed_ED95_mg = obs_ED95)
}
write.csv(do.call(rbind, grid_rows), "output/results/bruun_check_grid.csv", row.names = FALSE)
write.csv(do.call(rbind, bp_rows), "output/results/bruun_check.csv", row.names = FALSE)
freeze_stage("B", c("data/training/efficacy_training.csv", "data/priors/placebo_priors.csv",
                    sprintf("output/results/posterior_%s.rds", c("H0", "H1", "H2", "H3")),
                    cache_file, "output/results/binding_stage.txt"))
