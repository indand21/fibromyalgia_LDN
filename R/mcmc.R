# Adaptive Metropolis (Haario, Saksman & Tamminen 2001, Bernoulli 7:223).
am_sample <- function(logpost, init, n_iter = 20000, n_burn = 10000, seed = 1, adapt_start = 1000) {
  set.seed(seed)
  d <- length(init); x <- init; lp <- logpost(x)
  if (!is.finite(lp)) stop("log posterior not finite at init")
  draws <- matrix(NA_real_, n_iter, d, dimnames = list(NULL, names(init)))
  mu <- x; C <- matrix(0, d, d); prop_cov <- diag(0.01, d); sc <- 2.38^2 / d; acc <- 0
  for (i in seq_len(n_iter)) {
    L <- t(chol(sc * prop_cov + diag(1e-6, d)))
    y <- x + drop(L %*% stats::rnorm(d)); names(y) <- names(init)
    lpy <- logpost(y)
    if (is.finite(lpy) && log(stats::runif(1)) < lpy - lp) { x <- y; lp <- lpy; acc <- acc + 1 }
    draws[i, ] <- x
    delta <- x - mu; mu <- mu + delta / i; C <- C + tcrossprod(delta, x - mu)
    if (i >= adapt_start) prop_cov <- C / (i - 1)
  }
  list(draws = draws[(n_burn + 1):n_iter, , drop = FALSE], accept = acc / n_iter)
}

# ONE chain, seeded per CHAIN (not per worker): its jitter comes from seed + k and its sampler from seed + 100 * k,
# both under an explicitly installed Mersenne-Twister / Inversion / Rejection stream whose state is restored afterwards. A chain therefore
# depends only on (seed, k) -- not on which worker ran it, how many workers there were, or what the ambient RNG
# state and kind happened to be. That is what makes the parallel fit of script 04 bit-identical to the sequential one
# (a per-worker L'Ecuyer stream could not be: results would then depend on the worker count and task scheduling).
run_chain <- function(logpost, init, n_iter, n_burn, seed, k, jitter = 0.5) {
  with_mt_seed(seed + k, {
    x0 <- init + stats::rnorm(length(init), 0, jitter)
    while (!is.finite(logpost(x0))) x0 <- init + stats::rnorm(length(init), 0, jitter / 2)
    am_sample(logpost, x0, n_iter, n_burn, seed = seed + 100 * k)
  })
}

combine_chains <- function(chains) {
  n_chains <- length(chains)
  arr <- simplify2array(lapply(chains, `[[`, "draws"))          # iter x par x chain
  arr <- aperm(arr, c(1, 3, 2))                                 # iter x chain x par
  dr <- posterior::as_draws_array(arr)
  # plain numeric matrix (not draws_matrix) so that draws[i, ] is a named numeric vector and par[["mu_p"]] works
  list(draws = do.call(rbind, lapply(chains, `[[`, "draws")),
       summary = as.data.frame(posterior::summarise_draws(dr, "mean", "sd", ~quantile(.x, c(0.025, 0.5, 0.975)), "rhat", "ess_bulk")),
       accept = vapply(chains, `[[`, 0, "accept"),
       chain = rep(seq_len(n_chains), times = vapply(chains, function(ch) nrow(ch$draws), 0L)))
}

# map_chains lets a caller run the independent chains in parallel (for example
# `function(X, FUN) parallel::parLapply(cl, X, FUN)`); the result is identical to the sequential default.
run_chains <- function(logpost, init, n_chains = 4, n_iter = 20000, n_burn = 10000, seed = 2026, jitter = 0.5,
                       map_chains = lapply) {
  # Force every argument BEFORE the closure goes to map_chains: an unforced promise points at the caller's
  # environment (often the global one), which a PSOCK worker does not have (R2 fix round 1, Minor 1).
  force(logpost); force(init); force(n_iter); force(n_burn); force(seed); force(jitter)
  combine_chains(map_chains(seq_len(n_chains),
                            function(k) run_chain(logpost, init, n_iter, n_burn, seed, k, jitter)))
}
