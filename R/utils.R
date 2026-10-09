MW <- c(NTX = 341.4, BN = 343.4)  # g/mol; naltrexone, 6-beta-naltrexol

ngml_to_nM <- function(x, species) x * 1000 / MW[[species]]
nM_to_ngml <- function(x, species) x * MW[[species]] / 1000
mg_to_nmol <- function(mg, species) mg * 1e6 / MW[[species]]

logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

v5_root <- function() getOption("v5.root", normalizePath(".", winslash = "/"))

# Draw from a properly normalised normal TRUNCATED to [lower, upper] by inverse CDF (review I3).
# Censoring with pmin/pmax instead piles probability mass onto the bounds: with N(5.08, 1.5) on [4, 10] that put
# 23.7% of simulated baselines at exactly 4.0, moved the mean to 5.29 and shrank the SD to 1.21. No new
# dependency: stats::qnorm on a uniform draw restricted to (F(lower), F(upper)) is exact.
# A zero (or negative) SD is degenerate and returns the mean clamped into the interval.
rtnorm_trunc <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  if (any(lower >= upper)) stop("rtnorm_trunc: lower must be below upper (got ", lower, " and ", upper, ")")
  if (any(sd <= 0)) return(rep_len(pmin(pmax(mean, lower), upper), n))
  lo <- stats::pnorm(lower, mean, sd); hi <- stats::pnorm(upper, mean, sd)
  if (any(hi - lo <= 0))
    stop("rtnorm_trunc: the truncation interval has no probability mass under N(", mean, ", ", sd, ")")
  x <- stats::qnorm(stats::runif(n, lo, hi), mean, sd)
  pmin(pmax(x, lower), upper)            # guards the extreme tails of qnorm only; it is never the distribution
}

# Mean and SD of N(mu, sigma) truncated to [lower, upper] (either bound may be infinite).
trunc_moments <- function(mu, sigma, lower = -Inf, upper = Inf) {
  a <- (lower - mu) / sigma; b <- (upper - mu) / sigma
  # probability mass, computed in whichever tail keeps it accurate
  Z <- if (is.finite(a) && a > 0) stats::pnorm(-a) - stats::pnorm(-b) else stats::pnorm(b) - stats::pnorm(a)
  da <- stats::dnorm(a); db <- stats::dnorm(b)
  ada <- if (is.finite(a)) a * da else 0; bdb <- if (is.finite(b)) b * db else 0
  m <- mu + sigma * (da - db) / Z
  v <- sigma^2 * (1 + (ada - bdb) / Z - ((da - db) / Z)^2)
  c(mean = m, sd = sqrt(max(v, 0)))
}

# The PARENT normal whose truncation to [lower, upper] has a given mean and SD (R2 fix round 1, IMPORTANT 1).
# A trial's reported baseline mean and SD describe the ENROLLED population, which has already passed the
# eligibility screen; truncating N(reported mean, reported SD) would apply the screen a second time and shift
# both moments (Younger 2013 5.08/1.50 -> 5.68/1.10 on [4, 10]). So the parent is solved for numerically.
# If no truncated normal on the interval has those moments (for example a mean too close to a bound for its
# SD), the reported moments are used as the parent, flagged `matched = FALSE` with the reason in `note`.
match_trunc_parent <- function(mean, sd, lower, upper, tol = 1e-6) {
  if (!(mean > lower && mean < upper))
    stop("match_trunc_parent: reported mean ", mean, " is outside the bounds [", lower, ", ", upper, "]")
  if (sd <= 0) return(list(mu = mean, sigma = sd, matched = TRUE, lower = lower, upper = upper, residual = 0,
                           note = "degenerate SD: no truncation to match"))
  obj <- function(p) { m <- trunc_moments(p[1], exp(p[2]), lower, upper)
    if (!all(is.finite(m))) return(1e10)
    ((m[["mean"]] - mean) / sd)^2 + ((m[["sd"]] - sd) / sd)^2 }
  best <- NULL
  for (st in list(c(mean, log(sd)), c(mean, log(1.5 * sd)), c(mean - sd, log(2 * sd)), c(mean + sd, log(2 * sd)))) {
    o <- stats::optim(st, obj, control = list(reltol = 1e-14, maxit = 5000))
    o <- stats::optim(o$par, obj, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
    if (is.null(best) || o$value < best$value) best <- o
  }
  m <- trunc_moments(best$par[1], exp(best$par[2]), lower, upper)
  res <- max(abs(m[["mean"]] - mean), abs(m[["sd"]] - sd))
  if (is.finite(res) && res <= tol * max(1, sd))
    return(list(mu = best$par[1], sigma = exp(best$par[2]), matched = TRUE, lower = lower, upper = upper,
                residual = res, note = "parent moment-matched to the reported (enrolled) mean and SD"))
  list(mu = mean, sigma = sd, matched = FALSE, lower = lower, upper = upper, residual = res,
       note = sprintf(paste("reported moments %.3f/%.3f could not be matched by any truncated normal on [%g, %g]",
                            "(closest achievable %.3f/%.3f); reported moments used as the parent"),
                      mean, sd, lower, upper, m[["mean"]], m[["sd"]]))
}

source_all <- function(dir = file.path(v5_root(), "R")) {
  for (f in list.files(dir, pattern = "\\.R$", full.names = TRUE)) source(f)
  invisible(TRUE)
}
