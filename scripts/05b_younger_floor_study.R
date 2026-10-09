# OPTIONAL step after script 05 (R3a item 7): the seed-stability study that chose the Younger 2013 baseline floor
# recorded in data/params/baseline_bounds.csv (controller ruling, R2 fix round 2). It writes the file that
# data/README.md and baseline_bounds.csv cite, output/results/younger_floor_stability.csv (+ _raw.csv).
# Training data only: it reads no test data.
#
# Rule: for each candidate floor and each of N_SEEDS seeds, run the pipeline's own calibrate_omegas() on the H0
# posterior with that floor and seed, and record (a) the SD of the measurement-noise term 100 e / pain0 of the
# simulated placebo cohort (noise_sd_pbo, taken from calibrate_omegas' own cohort, not re-created by hand) and
# (b) the calibrated omega_p. A floor passes when the 5-95% range of both lies within +/-10% of their medians; the
# chosen floor is the smallest passing candidate in {0.5, 1, 1.5, 2} (0 is the reference that motivated the rule).
# The truncated-normal parent of each floor comes from match_trunc_parent() (calibrate_omegas draws the cohort
# with rtnorm_trunc() from that parent). Seeds 1001.. and n = 4000 per calibration, as in R2.
options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
N_SEEDS <- as.integer(Sys.getenv("N_SEEDS", "200")); FLOORS <- c(0, 0.5, 1, 1.5, 2); N_COHORT <- 4000
YOUNGER <- "younger2013"
stage <- readLines("output/results/binding_stage.txt")[1]
cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
arms <- read_dataset("data/training/efficacy_training.csv")
priors <- read_dataset("data/priors/placebo_priors.csv"); sd_res <- priors$value[priors$parameter == "sd_res_nrs"]
if (length(sd_res) != 1 || !is.finite(sd_res)) stop("priors need exactly one finite sd_res_nrs row")
bt0 <- read_baseline_bounds("data/params/baseline_bounds.csv")
fit <- readRDS("output/results/posterior_H0.rds")
# the reported Younger 2013 baseline moments, on the 0-10 scale, exactly as calibrate_omegas uses them
y_pbo <- arms[grepl("Younger.*2013", arms$trial, ignore.case = TRUE) & arms$arm == "placebo" & !is.na(arms$sd_pct), ][1, ]
base_mean <- if (identical(y_pbo$scale, "0-100")) y_pbo$baseline_nrs / 10 else y_pbo$baseline_nrs
base_sd <- if ("baseline_sd" %in% names(y_pbo) && is.finite(y_pbo$baseline_sd)) y_pbo$baseline_sd else 1.5
tasks <- expand.grid(seed = seq_len(N_SEEDS) + 1000L, floor = FLOORS)
# Built inside local() so the objects travel with the closure to the workers. Each task depends only on its
# (floor, seed): calibrate_omegas seeds its own cohort, and nothing else draws a random number.
job <- local({
  .fit <- fit; .arms <- arms; .cache <- cache; .bt0 <- bt0; .sd_res <- sd_res; .tasks <- tasks
  .bm <- base_mean; .bs <- base_sd; .n <- N_COHORT; .y <- YOUNGER
  function(i) {
    fl <- .tasks$floor[i]; s <- .tasks$seed[i]
    bt <- .bt0; bt$lower_nrs[bt$trial == .y] <- fl
    p <- match_trunc_parent(.bm, .bs, fl, 10)
    om <- suppressWarnings(calibrate_omegas(.fit, .arms, .cache, sd_res = .sd_res, n = .n, bounds = bt, seed = s))
    c(floor = fl, seed = s, noise_sd = om$noise_sd_pbo, omega_p = om$omega_p, matched = p$matched,
      parent_mu = p$mu, parent_sigma = p$sigma, min_base = om$min_baseline_pbo)
  }
})
t0 <- Sys.time()
n_workers <- min(10L, default_n_workers())
rows <- if (n_workers > 1) with_cluster(n_workers, function(cl) parallel::parLapply(cl, seq_len(nrow(tasks)), job)) else
  lapply(seq_len(nrow(tasks)), job)
res <- as.data.frame(do.call(rbind, rows))
cat(sprintf("%d calibrations on %d worker(s) in %.1f min\n", nrow(res), n_workers,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write.csv(res, "output/results/younger_floor_stability_raw.csv", row.names = FALSE)
summ <- do.call(rbind, lapply(FLOORS, function(fl) { r <- res[res$floor == fl, ]
  q <- function(v) c(stats::median(v), stats::quantile(v, c(0.05, 0.95)))
  a <- q(r$noise_sd); b <- q(r$omega_p)
  data.frame(floor = fl, parent_mu = r$parent_mu[1], parent_sigma = r$parent_sigma[1], matched = r$matched[1] == 1,
             noise_sd_median = a[1], noise_sd_q05 = a[2], noise_sd_q95 = a[3],
             noise_rel_lo = a[2] / a[1] - 1, noise_rel_hi = a[3] / a[1] - 1,
             omega_p_median = b[1], omega_p_q05 = b[2], omega_p_q95 = b[3],
             omega_rel_lo = b[2] / b[1] - 1, omega_rel_hi = b[3] / b[1] - 1,
             passes = abs(a[2] / a[1] - 1) <= 0.10 && abs(a[3] / a[1] - 1) <= 0.10 &&
                      abs(b[2] / b[1] - 1) <= 0.10 && abs(b[3] / b[1] - 1) <= 0.10, row.names = NULL) }))
write.csv(summ, "output/results/younger_floor_stability.csv", row.names = FALSE)
print(summ, digits = 4, row.names = FALSE)
cand <- summ[summ$floor > 0 & summ$passes, ]
chosen <- if (nrow(cand)) min(cand$floor) else NA_real_
cat("chosen floor (smallest candidate in {0.5, 1, 1.5, 2} meeting the rule):",
    if (is.na(chosen)) "NONE MEETS THE RULE" else chosen, "\n")
curated <- baseline_bounds_for(bt0, YOUNGER)[1]
if (is.na(chosen) || chosen != curated)
  warning("the stability rule chooses ", chosen, " but baseline_bounds.csv records ", curated, " for ", YOUNGER,
          "; the curated floor and this study disagree")
