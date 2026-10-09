options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
# Task 7 (PET test) must have run first: it decides the binding stage.
for (req in c("output/results/binding_stage.txt", "output/results/pet_metrics.csv"))
  if (!file.exists(req)) stop("missing ", req, "; run Task 7 (PET test) before fitting")
stage <- readLines("output/results/binding_stage.txt")[1]
cache_file <- sprintf("output/results/exposure_cache_%s.rds", stage)
cache <- readRDS(cache_file)
arms_all <- read_dataset("data/training/efficacy_training.csv")
# Arms without a reported dispersion cannot enter a likelihood that needs an SE: exclude (never impute).
# Every exclusion reason is recorded, including an endpoint mismatch stated in the data's own location column (review M2).
use <- is.finite(arms_all$se)
used_tab <- training_arms_table(arms_all)
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
HYPS <- c("H0", "H1", "H2", "H3")
N_CHAINS <- 4; N_ITER <- 60000; N_BURN <- 30000; MCMC_SEED <- 2026
init_for <- function(h) {
  init <- init_par(h, trials)
  if (h == "H3" && kd_fixed) init <- init[setdiff(names(init), "log10_Kd_OGFr")]
  init
}
# The 4 hypotheses x 4 chains are 16 INDEPENDENT runs, so they go out to a PSOCK cluster together (R2 item 5).
# Seeding is per CHAIN (run_chain installs its own Mersenne-Twister stream from MCMC_SEED + k), never per worker,
# so the fit does not depend on the number of workers or on how tasks were scheduled, and reproduces the
# sequential run bit for bit. chain_job is built inside local() so that every object it needs travels with the
# closure to the workers (a closure whose environment is the global one would arrive with nothing).
chain_job <- local({
  .arms <- arms; .cache <- cache; .priors <- priors; .fixed <- fixed; .trials <- trials
  .kd_fixed <- kd_fixed; .n_iter <- N_ITER; .n_burn <- N_BURN; .seed <- MCMC_SEED
  function(task) {
    h <- task[["h"]]; k <- as.integer(task[["k"]])
    lp <- function(p) { pr <- log_prior(p, h, .priors, .trials); if (!is.finite(pr)) return(-Inf)
                        pr + log_lik_arms(p, h, .arms, .cache, .fixed) }
    init <- init_par(h, .trials)
    if (h == "H3" && .kd_fixed) init <- init[setdiff(names(init), "log10_Kd_OGFr")]
    run_chain(lp, init, .n_iter, .n_burn, .seed, k)
  }
})
tasks <- do.call(rbind, lapply(HYPS, function(h) data.frame(h = h, k = seq_len(N_CHAINS), stringsAsFactors = FALSE)))
task_list <- lapply(seq_len(nrow(tasks)), function(i) as.list(tasks[i, ]))
n_workers <- min(nrow(tasks), default_n_workers())
t_mcmc <- Sys.time()
chains <- if (n_workers > 1)
  with_cluster(n_workers, function(cl) parallel::parLapply(cl, task_list, chain_job)) else lapply(task_list, chain_job)
cat(sprintf("%d chains (%d hypotheses x %d chains) on %d worker(s) in %.1f min\n",
            nrow(tasks), length(HYPS), N_CHAINS, n_workers,
            as.numeric(difftime(Sys.time(), t_mcmc, units = "mins"))))
fits <- list(); loo_rows <- list()
for (h in HYPS) {
  f <- combine_chains(chains[tasks$h == h])
  f$hyp <- h; f$fixed <- fixed; f$fixed_sources <- fixed_sources; f$trials <- trials; f$stage <- stage
  f$n_params <- length(init_for(h))
  stopifnot(f$n_params == ncol(f$draws))
  # pointwise terms are the leave-one-out conditionals p(y_i | y_-i, theta) (R2 fix round 1)
  loo_rows[[h]] <- training_loo(f, h, arms, cache, fixed, n_keep = 2000)
  saveRDS(f, sprintf("output/results/posterior_%s.rds", h)); fits[[h]] <- f
}
loo_tab <- do.call(rbind, loo_rows)
write.csv(loo_tab, "output/results/loo_training.csv", row.names = FALSE)

# R2 item 7 (review I5): p_loo against the free-parameter count. p_loo is the effective number of parameters the
# TWO training observations actually constrain, so 1 - p_loo/n_params is the fraction of the parameter space that
# the posterior takes from the priors rather than from the data. A citable file of its own; pure arithmetic.
p_loo_tab <- p_loo_table(loo_tab, n_obs = nrow(arms))   # per-row reliability note (R2 fix round 2)
write.csv(p_loo_tab, "output/results/p_loo_identifiability.csv", row.names = FALSE)
print(p_loo_tab[, c("hypothesis", "n_free_parameters", "p_loo", "prior_determined_fraction", "max_pareto_k", "p_loo_reliable")])
fit_tab <- do.call(rbind, lapply(fits, function(f) cbind(hypothesis = f$hyp, f$summary)))
write.csv(fit_tab, "output/results/training_fit.csv", row.names = FALSE)

# R2 item 2 (review I9): prior vs posterior SD of the placebo parameters, so the manuscript can state plainly
# which of them the data move and which the priors determine. mu_p carries the meta-analytic pooled mean with its
# SE (a prior LOCATION); between-trial heterogeneity sits in tau ~ half-normal(0.5), and a trial's ceiling is
# inv_logit(mu_p + eta_trial). A half-normal(s) has mean s*sqrt(2/pi) and SD s*sqrt(1 - 2/pi).
hn <- function(s) c(mean = s * sqrt(2 / pi), sd = s * sqrt(1 - 2 / pi))
kp_prior <- log_scale_prior(priors, "kpbo_per_day")
mu_prior_sd <- priors$value[priors$parameter == "pbo_logit_sd"]
pp <- do.call(rbind, lapply(fits, function(f) {
  d <- f$draws; tau <- exp(d[, "log_tau"]); sig <- exp(d[, "log_sigma_rel"])
  data.frame(hypothesis = f$hyp,
    mu_p_prior_sd = mu_prior_sd, mu_p_post_mean = mean(d[, "mu_p"]), mu_p_post_sd = stats::sd(d[, "mu_p"]),
    tau_prior_mean = hn(0.5)[["mean"]], tau_prior_sd = hn(0.5)[["sd"]],
    tau_post_mean = mean(tau), tau_post_sd = stats::sd(tau),
    log_kpbo_prior_mean = kp_prior$mean, log_kpbo_prior_sd = kp_prior$sd,
    log_kpbo_post_mean = mean(d[, "log_kpbo"]), log_kpbo_post_sd = stats::sd(d[, "log_kpbo"]),
    sigma_rel_prior_sd = hn(0.05)[["sd"]], sigma_rel_post_mean = mean(sig), sigma_rel_post_sd = stats::sd(sig),
    kpbo_prior_scale_note = sprintf("log-scale SD %.4f converted from the per-day SE %.5f / value %.3f (review C2)",
                                    kp_prior$sd, priors$sd[priors$parameter == "kpbo_per_day"],
                                    priors$value[priors$parameter == "kpbo_per_day"]))
}))
pp$mu_p_post_over_prior_sd <- pp$mu_p_post_sd / pp$mu_p_prior_sd
pp$tau_post_over_prior_sd <- pp$tau_post_sd / pp$tau_prior_sd
pp$log_kpbo_post_over_prior_sd <- pp$log_kpbo_post_sd / pp$log_kpbo_prior_sd
write.csv(pp, "output/results/placebo_prior_vs_posterior.csv", row.names = FALSE)
print(pp[, c("hypothesis", "mu_p_post_sd", "mu_p_prior_sd", "tau_post_mean", "tau_post_sd",
             "log_kpbo_post_sd", "log_kpbo_prior_sd")], digits = 4)

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
# Bruun-Plesner consistency check: effect vs dose at day 365, with the model ED50 AND ED95 and the position of the
# peak on the grid (R3a, review M7). The grid reaches down to 0.1 mg because H1's effect peaked at 0.75 mg, the
# lowest non-zero dose of the old grid; every grid dose has an exposure-cache profile (scripts/02 builds them).
BRUUN_DOSES <- c(0, 0.1, 0.25, 0.5, 0.75, 1.5, 2, 3, 3.88, 4.5, 5, 5.4, 6)
bp <- bruun_check(fits[c("H1", "H2", "H3")], cache, BRUUN_DOSES, obs_ED50, obs_ED95, day = 365)
write.csv(bp$grid, "output/results/bruun_check_grid.csv", row.names = FALSE)
write.csv(bp$summary, "output/results/bruun_check.csv", row.names = FALSE)
print(bp$summary[, c("hypothesis", "peak_dose_mg", "peak_position", "model_ED50_mg", "model_ED95_mg", "observed_ED95_mg")])
freeze_stage("B", c("data/training/efficacy_training.csv", "data/priors/placebo_priors.csv",
                    sprintf("output/results/posterior_%s.rds", c("H0", "H1", "H2", "H3")),
                    cache_file, "output/results/binding_stage.txt"))
