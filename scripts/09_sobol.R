options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
stage <- readLines("output/results/binding_stage.txt")[1]; cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
te <- read.csv("output/results/test_elpd.csv")   # selection uses the PRIMARY elpd column (one endpoint per trial)
best <- te$hypothesis[which.max(te$elpd_primary)]
note <- ""
if (best == "H0") {   # no drug effect to decompose under H0: use the best drug hypothesis
  d <- te[te$hypothesis != "H0", ]; best <- d$hypothesis[which.max(d$elpd_primary)]
  note <- "H0 best; Sobol on best drug hypothesis"
}
f <- readRDS(sprintf("output/results/posterior_%s.rds", best))
# Varied set: mechanism parameters only. Sigma is an observation-level term; mu_p, tau and eta are population/trial terms.
free <- setdiff(colnames(f$draws), c("mu_p", "log_tau", "log_sigma_rel", grep("^eta_", colnames(f$draws), value = TRUE)))
q <- apply(f$draws[, free, drop = FALSE], 2, stats::quantile, c(0.025, 0.975))
med <- apply(f$draws, 2, stats::median)                # full par from posterior medians; only `free` are overwritten
pm <- inv_logit(med[["mu_p"]])
net <- function(x) {
  par <- med; par[free] <- as.numeric(x)
  th <- theta_from(par, best, f$fixed)
  a <- simulate_arm(best, th, parse_schedule("0:4.5"), 84, cache, pm)$rel[85]
  b <- simulate_arm(best, th, parse_schedule("0:0"), 84, cache, pm)$rel[85]
  100 * (b - a)
}
# sensitivity::soboljansen needs d >= 2 (it fails on a single column). Every hypothesis always varies kout and kpbo plus Emax, so d >= 3.
if (length(free) < 2) stop("Sobol needs at least 2 varied parameters, got: ", paste(free, collapse = ", "))
set.seed(3); N <- 2000
box <- function() { x <- as.data.frame(lapply(free, function(p) stats::runif(N, q[1, p], q[2, p]))); names(x) <- free; x }
t1 <- system.time(for (i in 1:20) net(q[1, ]))[["elapsed"]] / 20
message(sprintf("hypothesis %s, d = %d, %.1f ms/evaluation, ~%.1f min for %d evaluations", best, length(free), 1000 * t1,
                t1 * N * (length(free) + 2) / 60, N * (length(free) + 2)))
sa <- sensitivity::soboljansen(model = NULL, X1 = box(), X2 = box(), nboot = 200)
sensitivity::tell(sa, apply(sa$X, 1, net))
out <- data.frame(parameter = free, S1 = sa$S$original, S1_lo = sa$S$`min. c.i.`, S1_hi = sa$S$`max. c.i.`,
                  ST = sa$T$original, ST_lo = sa$T$`min. c.i.`, ST_hi = sa$T$`max. c.i.`, hypothesis = best, note = note)
write.csv(out, "output/results/sobol.csv", row.names = FALSE)
