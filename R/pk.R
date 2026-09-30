# Starting values only; final values come from fit_pk() on data/training/pk_summary.csv.
# fu values are overwritten from data/params/binding_constants.csv in 01_fit_pk.R.
pk_default <- function() list(ka = 1.2, fp = 0.6, CL = 200, V1 = 400, Q = 100, V2 = 800,
                              fm = 1, CLm = 60, Vm = 1100, fu_NTX = 0.75, fu_BN = 0.8)

# Doses are the labelled naltrexone hydrochloride mass; all PK parameters are apparent
# (per mg labelled dose, oral), so no salt correction is applied.
# fm is fixed at 1 because naltrexone is predominantly converted to 6-beta-naltrexol,
# so CLm is an apparent CLm/fm.
# Apparent (oral) parameters. fp = fraction of absorbed dose converted to 6-beta-naltrexol
# presystemically; fm = fraction of systemic parent clearance forming 6-beta-naltrexol.
pk_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    absorb <- ka * Ad
    el <- CL / V1 * A1
    list(c(-absorb,
           absorb * (1 - fp) - el - Q / V1 * A1 + Q / V2 * A2,
           Q / V1 * A1 - Q / V2 * A2,
           absorb * fp + fm * el - CLm / Vm * Am))
  })
}

simulate_pk <- function(dose_mg, times, p = pk_default(), tau = 24, n_doses = 1) {
  dose_times <- (seq_len(n_doses) - 1) * tau
  times <- sort(unique(c(0, times, dose_times)))
  if (dose_mg == 0) {
    z <- rep(0, length(times))
    return(data.frame(time = times, Cp_NTX = z, Cp_BN = z, Cu_NTX = z, Cu_BN = z))
  }
  ev <- data.frame(var = "Ad", time = dose_times, value = mg_to_nmol(dose_mg, "NTX"), method = "add")
  out <- as.data.frame(deSolve::lsoda(c(Ad = 0, A1 = 0, A2 = 0, Am = 0), times, pk_ode, p,
                                      events = list(data = ev), rtol = 1e-8, atol = 1e-10))
  cp <- out$A1 / p$V1; cm <- out$Am / p$Vm
  data.frame(time = out$time, Cp_NTX = cp, Cp_BN = cm, Cu_NTX = p$fu_NTX * cp, Cu_BN = p$fu_BN * cm)
}

nca_summary <- function(time, conc, n_terminal = 3) {
  i <- which.max(conc)
  auc_last <- sum(diff(time) * (head(conc, -1) + tail(conc, -1)) / 2)
  idx <- tail(which(conc > 0), n_terminal)
  lz <- -unname(stats::coef(stats::lm(log(conc[idx]) ~ time[idx]))[2])
  c(Cmax = conc[i], Tmax = time[i], AUCinf = auc_last + tail(conc, 1) / lz, thalf = log(2) / lz)
}

pk_predict_stats <- function(p, obs) {
  bad_stat <- setdiff(unique(obs$stat), c("Cmax", "Tmax", "AUCinf", "thalf"))
  bad_an <- setdiff(unique(obs$analyte), c("NTX", "BN"))
  if (length(bad_stat)) stop("Unknown stat value(s): ", paste(bad_stat, collapse = ", "))
  if (length(bad_an)) stop("Unknown analyte value(s): ", paste(bad_an, collapse = ", "))
  pred <- numeric(nrow(obs))
  for (d in unique(obs$dose_mg)) {
    s <- simulate_pk(d, seq(0, 96, 0.25), p)
    nca <- list(NTX = nca_summary(s$time, nM_to_ngml(s$Cp_NTX, "NTX")),
                BN = nca_summary(s$time, nM_to_ngml(s$Cp_BN, "BN")))
    rows <- which(obs$dose_mg == d)
    pred[rows] <- vapply(rows, function(k) unname(nca[[obs$analyte[k]]][obs$stat[k]]), 0)
  }
  pred
}

PK_FIT_PARS <- c("ka", "fp", "CL", "V1", "Q", "V2", "CLm", "Vm")

fit_pk <- function(obs, start = pk_default(), pars = PK_FIT_PARS) {
  obs <- obs[obs$stat != "Tmax", , drop = FALSE]  # Tmax reflects sampling grids; used as a check only
  cv <- pmax(obs$sd / obs$mean, 0.1, na.rm = TRUE)
  to_p <- function(th) {
    p <- start
    for (nm in pars) p[[nm]] <- if (nm == "fp") inv_logit(th[[nm]]) else exp(th[[nm]])
    p
  }
  th0 <- vapply(pars, function(nm) if (nm == "fp") logit(start[[nm]]) else log(start[[nm]]), 0)
  obj <- function(th) {
    names(th) <- pars
    pred <- pk_predict_stats(to_p(th), obs)
    if (any(!is.finite(pred) | pred <= 0)) return(1e10)
    sum(((log(pred) - log(obs$mean)) / cv)^2)
  }
  o1 <- stats::optim(th0, obj, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-10))
  o2 <- stats::optim(o1$par, obj, method = "BFGS", hessian = TRUE)
  list(par = to_p(stats::setNames(o2$par, pars)), value = o2$value, convergence = o2$convergence,
       hessian = o2$hessian, pars = pars)
}
