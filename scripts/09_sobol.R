options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
stage <- readLines("output/results/binding_stage.txt")[1]; cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
# Review I2: the hypothesis is no longer picked by elpd with a fallback only when it is NAMED H0 (which decomposed
# H2's 0.001% effect). Every hypothesis whose posterior-median net effect at day 84 exceeds the threshold is analysed,
# at a non-saturated (day 14) and a saturated (day 84) timepoint. A hypothesis below the threshold gets a row that says
# so instead of being skipped silently.
# R3a fix round 1: the outputs are heavy-tailed (kurtosis about 9-10 at H1 day 14), so N = 20000 and every run is
# repeated with a second seed; sobol_seed_agreement.csv flags any index whose two estimates differ by more than the
# larger bootstrap CI half-width. sobol.csv carries the first seed's estimates.
EFFECT_THRESHOLD_PP <- 0.5; SOBOL_DAYS <- c(14, 84); N <- 20000; NBOOT <- 200; SOBOL_SEEDS <- c(3, 17)
HYPS <- c("H0", "H1", "H2", "H3")
fits <- lapply(stats::setNames(nm = HYPS), function(h) readRDS(sprintf("output/results/posterior_%s.rds", h)))
el <- effect_eligibility(fits, cache, day = 84, threshold_pp = EFFECT_THRESHOLD_PP)
print(el[, c("hypothesis", "net_effect_pp", "eligible", "sobol_role")])
# Model evaluations are deterministic, so they go to a PSOCK cluster in chunks; the sample matrices and the bootstrap
# use their own seeded streams in the master (sobol_hypothesis), so the result does not depend on the worker count.
n_workers <- default_n_workers()
t0 <- Sys.time()
run <- function(map, seed) lapply(HYPS, function(h) {
  e <- el[el$hypothesis == h, ]
  if (!e$eligible) return(sobol_skipped(h, SOBOL_DAYS, net_effect_at(fits[[h]], cache, SOBOL_DAYS)))
  sobol_hypothesis(fits[[h]], cache, days = SOBOL_DAYS, N = N, nboot = NBOOT, seed = seed, map = map)
})
runs <- if (n_workers > 1) with_cluster(n_workers, function(cl) {
  map <- function(X, FUN) parallel::parLapplyLB(cl, X, FUN)
  lapply(SOBOL_SEEDS, function(s) do.call(rbind, run(map, s)))
}) else lapply(SOBOL_SEEDS, function(s) do.call(rbind, run(lapply, s)))
cat(sprintf("Sobol (N = %d, seeds %s) on %d worker(s) in %.1f min\n", N, paste(SOBOL_SEEDS, collapse = " and "), n_workers,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
agree <- sobol_seed_agreement(runs[[1]], runs[[2]], seeds = SOBOL_SEEDS)
agree$note <- sprintf("disagree = |estimate_seed%d - estimate_seed%d| > the larger of the two bootstrap 95%% CI half-widths (N = %d)",
                      SOBOL_SEEDS[1], SOBOL_SEEDS[2], N)
write.csv(agree, "output/results/sobol_seed_agreement.csv", row.names = FALSE)
cat(sprintf("seed agreement: %d of %d indices disagree beyond the larger CI half-width\n", sum(agree$disagree), nrow(agree)))
out <- runs[[1]]
flag <- function(ix) { a <- agree[agree$index == ix, ]
  a$disagree[match(paste(out$hypothesis, out$day, out$parameter), paste(a$hypothesis, a$day, a$parameter))] }
out$seed_disagree_S1 <- flag("S1"); out$seed_disagree_ST <- flag("ST")
out <- merge(out, el[, c("hypothesis", "net_effect_pp", "threshold_pp", "eligible", "sobol_role")], by = "hypothesis", sort = FALSE)
names(out)[names(out) == "net_effect_pp"] <- "net_effect_day84_pp"
out$note <- ifelse(out$eligible,
  paste("Sobol indices of the net effect (pp of pain reduction, posterior medians elsewhere) at this day: S1 by the Sobol 2007",
        "estimator, ST by Jansen, on the same evaluations (Jansen's S1, kept as S1_jansen, is high-variance here: V from one",
        "sample, heavy-tailed output); parameters uniform over their marginal 95% posterior intervals (posterior correlation",
        "ignored, review M4); seed_disagree_* compare with a second seed (sobol_seed_agreement.csv)"),
  out$sobol_role)
write.csv(out, "output/results/sobol.csv", row.names = FALSE)
print(out[, c("hypothesis", "day", "parameter", "S1", "ST", "seed_disagree_S1", "seed_disagree_ST", "eligible")], digits = 3)
