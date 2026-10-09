# Helpers for scripts/12_derived_quantities.R: quantities the manuscript cites that no earlier script writes (they were
# first computed by hand in the R4 run verification). Every function works on objects already saved in output/results;
# none reads test data, refits anything or changes a method.

DERIVED_COLS <- c("quantity", "value", "lo", "hi", "unit", "definition", "source_files")

# One row of derived_quantities.csv. `source_files` is a character vector, stored "; "-separated.
dq_row <- function(quantity, value, lo = NA_real_, hi = NA_real_, unit, definition, source_files)
  data.frame(quantity = quantity, value = as.numeric(value), lo = as.numeric(lo), hi = as.numeric(hi), unit = unit,
             definition = definition, source_files = paste(source_files, collapse = "; "), stringsAsFactors = FALSE)

# A saved posterior (draws stacked by chain, `chain` labels each row) as an iterations x chains x parameters array.
posterior_draws_array <- function(fit) {
  ch <- sort(unique(fit$chain))
  n <- unique(vapply(ch, function(k) sum(fit$chain == k), 0L))
  if (length(n) != 1) stop("posterior_draws_array: chains have unequal lengths", call. = FALSE)
  arr <- array(NA_real_, c(n, length(ch), ncol(fit$draws)), dimnames = list(NULL, NULL, colnames(fit$draws)))
  for (j in seq_along(ch)) arr[, j, ] <- fit$draws[fit$chain == ch[j], , drop = FALSE]
  posterior::as_draws_array(arr)
}

# Smallest tail ESS (posterior::ess_tail) over the sampled parameters, and the parameter that attains it.
min_ess_tail <- function(fit) {
  sm <- posterior::summarise_draws(posterior_draws_array(fit), "ess_tail")
  i <- which.min(sm$ess_tail)
  list(value = sm$ess_tail[[i]], parameter = sm$variable[[i]])
}

# The PK Laplace draws of script 03's binding + PK band (laplace_draws(pk$theta, cov, n, seed)), back-transformed:
# fp0 and the oral bioavailability at `dose_mg` for every draw.
pk_laplace_draw_summary <- function(pk, n, seed, dose_mg = 4.5) {
  d <- laplace_draws(pk$theta, pk$hessian_summary$cov, n, seed = seed)
  ps <- lapply(seq_len(n), function(i) pk_from_theta(d[i, ], pk$par, pk$pars))
  list(fp0 = vapply(ps, function(p) if (is.null(p$fp0)) NA_real_ else p$fp0, 0),
       F = vapply(ps, function(p) pk_bioavailability(dose_mg, p), 0))
}

# Evenly spaced (thinned) row indices over all saved draws, chains in their stored order.
thinned_index <- function(n_total, n_keep) round(seq(1, n_total, length.out = min(n_keep, n_total)))

# Net effect (net_effect_par, the definition behind effect_eligibility.csv) over thinned posterior draws.
net_effect_draws <- function(fit, cache, day, n_keep = 400)
  vapply(thinned_index(nrow(fit$draws), n_keep),
         function(i) net_effect_par(fit$draws[i, ], fit$hyp, fit$fixed, cache, day), 0)

# Mean unbound plasma naltrexone (nM) over one steady-state dosing interval [0, tau) of a cached exposure profile
# (the endpoint duplicates hour 0, as in profile_means()).
mean_unbound_ntx <- function(prof) mean(prof$Cu_NTX[prof$hour < max(prof$hour)])
