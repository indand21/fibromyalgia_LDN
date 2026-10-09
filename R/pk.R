# Starting values only; final values come from fit_pk() on data/training/pk_summary.csv.
# fu values are overwritten from data/params/binding_constants.csv in 01_fit_pk.R.
# fp      = fraction of the absorbed dose converted to 6-beta-naltrexol presystemically AT THE REFERENCE DOSE fp_Dref
#           (100 mg, the Verebey 1976 dose); in the linear model it is the same at every dose.
# fp0     = presystemic conversion fraction as dose -> 0 (dose-dependent first pass, structure "fp_dose"); NULL = off.
# floss0  = fraction of the presystemically ESCAPING parent that is lost without forming 6-beta-naltrexol, as dose -> 0
#           (structure "parent_loss"); NULL = off.
# fp_D50  = dose (mg) of half-saturation in the dose shape h(D) = D / (D + D50); Inf = dose dependence switched off.
# fp_Dref = reference (anchor) dose, mg: fp(fp_Dref) = fp and no parent is lost at fp_Dref, exactly, for every D50.
pk_default <- function() list(ka = 1.2, fp = 0.6, fp0 = NULL, floss0 = NULL, fp_D50 = Inf, fp_Dref = FP_DREF_MG,
                              CL = 200, V1 = 400, Q = 100, V2 = 800,
                              fm = 1, CLm = 60, Vm = 1100, fu_NTX = 0.75, fu_BN = 0.8)

# ---- Dose-dependent first pass (R4; reparameterised in R4 fix round 1, review I1/I2/I6) ----
# Both alternatives to the linear model use the same saturable dose shape h(D) = D / (D + D50), NORMALISED AT THE
# REFERENCE DOSE: r(D) = h(D) / h(Dref), so r(Dref) = 1 for every D50 and r(D) -> 0 as D -> 0. The PK at 100 mg is
# therefore identical to the linear model at 100 mg for every D50 (fp at the Verebey dose stays the documented 0.95), and
# a D50 sensitivity changes only the shape between 0 and 100 mg.
#   fp_dose     (one extra parameter fp0): the presystemically escaping fraction F_esc(D) = 1 - fp(D) runs from
#               F_esc0 = 1 - fp0 at D -> 0 to F_esc_ref = 1 - fp at Dref:  F_esc(D) = F_esc0 + (F_esc_ref - F_esc0) r(D).
#               fp0 > fp means more conversion at low dose: LESS parent and slightly MORE 6-beta-naltrexol early on.
#   parent_loss (one extra parameter floss0): fp is unchanged; a fraction L(D) = floss0 (1 - r(D)) of the escaping parent
#               is lost presystemically WITHOUT forming 6-beta-naltrexol: less parent, and no extra metabolite.
# Above Dref the shape is extrapolated (r > 1): F_esc is clamped to (0, 1) and L to >= 0 (no loss at or above 100 mg).
# fp(D) and L(D) are computed ONCE per simulated dose, so the ODE stays linear in its states at a fixed dose and daily
# dosing by superposition remains valid. Every PK path (fits, exposure caches, PET predictions) goes through simulate_pk().
FP_DREF_MG <- 100
PK_DOSE_MODELS <- list(
  linear      = list(par = NULL,     label = "linear (fp dose-independent)"),
  fp_dose     = list(par = "fp0",    label = "dose-dependent presystemic conversion fp(D), anchored at 100 mg"),
  parent_loss = list(par = "floss0", label = "dose-dependent parent-only presystemic loss, anchored at 100 mg"))
pk_dose_shape <- function(dose_mg, D50, Dref = FP_DREF_MG) (dose_mg / (dose_mg + D50)) / (Dref / (Dref + D50))
pk_dref <- function(p) if (is.null(p$fp_Dref)) FP_DREF_MG else p$fp_Dref
pk_dd_on <- function(p, par) !is.null(p[[par]]) && !is.null(p$fp_D50) && is.finite(p$fp_D50) && p$fp_D50 > 0
fp_at_dose <- function(dose_mg, p) {
  if (!pk_dd_on(p, "fp0") || isTRUE(p$fp0 == p$fp)) return(p$fp)
  F_ref <- 1 - p$fp; F_0 <- 1 - p$fp0
  F_esc <- F_0 + (F_ref - F_0) * pk_dose_shape(dose_mg, p$fp_D50, pk_dref(p))
  1 - pmin(pmax(F_esc, 1e-12), 1 - 1e-12)
}
ploss_at_dose <- function(dose_mg, p) {
  if (!pk_dd_on(p, "floss0")) return(0)
  pmax(p$floss0 * (1 - pk_dose_shape(dose_mg, p$fp_D50, pk_dref(p))), 0)
}
# Systemic (oral) bioavailability of the parent: the fraction of the dose reaching the circulation as naltrexone.
pk_bioavailability <- function(dose_mg, p) (1 - fp_at_dose(dose_mg, p)) * (1 - ploss_at_dose(dose_mg, p))
pk_dose_model_of <- function(p) if (pk_dd_on(p, "fp0")) "fp_dose" else if (pk_dd_on(p, "floss0")) "parent_loss" else "linear"

# Doses are the labelled naltrexone hydrochloride mass; all PK parameters are apparent
# (per mg labelled dose, oral), so no salt correction is applied.
# fm is fixed at 1 because naltrexone is predominantly converted to 6-beta-naltrexol,
# so CLm is an apparent CLm/fm.
# Apparent (oral) parameters. fp = fraction of absorbed dose converted to 6-beta-naltrexol
# presystemically; fm = fraction of systemic parent clearance forming 6-beta-naltrexol.
# fpar (set by simulate_pk) = 1 - L(D), the fraction of the escaping parent NOT lost presystemically (1 unless the
# parent_loss structure is on).
pk_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    absorb <- ka * Ad
    el <- CL / V1 * A1
    list(c(-absorb,
           absorb * (1 - fp) * fpar - el - Q / V1 * A1 + Q / V2 * A2,
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
  # Dose-dependent fp: compute fp(D) once for this dose; the ODE is linear in state variables so
  # superposition for daily dosing remains valid (fp is constant within a given dose level).
  p_eff <- p; p_eff$fp <- fp_at_dose(dose_mg, p); p_eff$fpar <- 1 - ploss_at_dose(dose_mg, p)
  p_eff$fp0 <- NULL; p_eff$floss0 <- NULL
  ev <- data.frame(var = "Ad", time = dose_times, value = mg_to_nmol(dose_mg, "NTX"), method = "add")
  out <- as.data.frame(deSolve::lsoda(c(Ad = 0, A1 = 0, A2 = 0, Am = 0), times, pk_ode, p_eff,
                                      events = list(data = ev), rtol = 1e-8, atol = 1e-10))
  cp <- out$A1 / p_eff$V1; cm <- out$Am / p_eff$Vm
  data.frame(time = out$time, Cp_NTX = cp, Cp_BN = cm, Cu_NTX = p_eff$fu_NTX * cp, Cu_BN = p_eff$fu_BN * cm)
}

# Cmax/Tmax are refined by a parabola through the discrete maximum and its two neighbours, so that the predicted Tmax
# varies continuously with ka (a Tmax read off the simulation grid is a step function of ka and makes the objective flat).
refine_peak <- function(time, conc) {
  if (any(!is.finite(conc))) return(c(Cmax = NA_real_, Tmax = NA_real_))   # failed or truncated simulation
  i <- which.max(conc)
  if (i == 1 || i == length(conc)) return(c(Cmax = conc[i], Tmax = time[i]))
  x <- time[(i - 1):(i + 1)]; y <- conc[(i - 1):(i + 1)]
  b <- solve(cbind(1, x, x^2), y)
  if (!is.finite(b[3]) || b[3] >= 0) return(c(Cmax = conc[i], Tmax = time[i]))
  tm <- -b[2] / (2 * b[3])
  c(Cmax = unname(b[1] + b[2] * tm + b[3] * tm^2), Tmax = unname(tm))
}

# `thalf` is the slope of log concentration over the LAST `n_terminal` grid points of the simulated window, i.e. the final
# 2 * step hours (0.5 h at the 0.25 h spacing used here), not a wide late window. For these parameters it is the asymptotic
# terminal slope to 3 decimals, but it is NOT the same estimand as a half-life a paper derives from early sampling: Verebey
# 1976's 9.7-10.3 h naltrexone half-lives come from the "biexponential kinetics in the first 24 hr" and so mix the
# distribution and terminal phases. That mismatch, not a numerical problem, is what makes the terminal-phase sensitivity
# (R1 fix round 3) unable to fit both those half-lives and the 24-72 h tail.
nca_summary <- function(time, conc, n_terminal = 3) {
  if (any(!is.finite(conc))) return(c(Cmax = NA_real_, Tmax = NA_real_, AUCinf = NA_real_, thalf = NA_real_))
  pk <- refine_peak(time, conc)
  auc_last <- sum(diff(time) * (head(conc, -1) + tail(conc, -1)) / 2)
  idx <- tail(which(conc > 0), n_terminal)
  lz <- -unname(stats::coef(stats::lm(log(conc[idx]) ~ time[idx]))[2])
  c(Cmax = pk[["Cmax"]], Tmax = pk[["Tmax"]], AUCinf = auc_last + tail(conc, 1) / lz, thalf = log(2) / lz)
}

# Fine grid over the absorption phase (Tmax resolution), coarser afterwards; for `daily` rows the same grid is laid over
# the interval after the LAST dose (times relative to that dose).
PK_STAT_GRID <- sort(unique(c(seq(0, 12, 0.02), seq(12, 96, 0.25))))
PK_REGIMENS <- c("single", "daily")
PK_STATS <- c("Cmax", "Tmax", "AUCinf", "thalf", "C24", "C72", "C0.25", "C0.5", "C0.75", "C1", "Css13h")
# C0.25, C0.5, C0.75, C1: total plasma concentration exactly 15, 30, 45 and 60 min after the (last) dose, analogous to
#         C24 / C72 (Bested 2023 day-1 profile; R4 fix round 1 added 15, 45 and 60 min). All are points of PK_STAT_GRID.
# Css13h: mean total plasma concentration over the 1-3 h window after the (last) dose (trapezoid rule; defined for both
#         single and daily regimens but physically meaningful mainly for the daily steady-state case).
PK_MIN_DAILY_DOSES <- 14
PK_MAX_DAILY_DOSES <- 60
# Linear superposition for `daily` rows: the PK ODE is linear and dosing additive, so the profile after the last of n daily
# doses is the sum of the single-dose profile shifted by 0, tau, ..., (n-1) tau. One event-free single-dose solve therefore
# serves both the single and the daily rows (about 3x faster than a separate event-driven multi-dose solve). Gated by
# test-pk.R to 1e-5 relative against the event-driven path; set to FALSE to fall back to that path.
PK_SUPERPOSITION <- TRUE

# Number of once-daily doses simulated for a `daily` row: at least PK_MIN_DAILY_DOSES, and enough for 5 of the slowest
# disposition half-lives to elapse before the last dose, so the last interval is at steady state. The parent half-lives come
# from pk_parent_half_lives() (the true two-compartment beta phase, not Q/V2, which is not a half-life at all).
pk_n_daily_doses <- function(p, tau = 24) {
  hl <- c(pk_parent_half_lives(p), metabolite = log(2) / (p$CLm / p$Vm), absorption = log(2) / p$ka)
  t_half <- max(hl[is.finite(hl)], 0)
  as.integer(min(PK_MAX_DAILY_DOSES, max(PK_MIN_DAILY_DOSES, ceiling(5 * t_half / tau) + 1)))
}

# Summary statistics per row. `regimen` (optional column; absent means "single"):
#   single: one dose at time 0; NCA over 0-96 h.
#   daily : once-daily dosing (tau 24 h) to steady state; Cmax/Tmax over the interval after the last dose, terminal
#           half-life from the log-linear decay over the last grid points of the 0-96 h window after the last dose (the same
#           window as a single dose; see nca_summary above). AUCinf is not defined for a daily row and stops.
#   C24   : concentration exactly 24 h after the (last) dose, from the same simulation (24 h is a grid point); for a
#           daily row this is the steady-state trough. C72: likewise at 72 h (no further dose; used by the terminal-phase
#           sensitivity analysis only).
pk_predict_stats <- function(p, obs, superposition = PK_SUPERPOSITION, tau = 24) {
  reg <- if (is.null(obs$regimen)) rep("single", nrow(obs)) else obs$regimen
  bad_stat <- setdiff(unique(obs$stat), PK_STATS)
  bad_an <- setdiff(unique(obs$analyte), c("NTX", "BN"))
  bad_reg <- setdiff(unique(reg), PK_REGIMENS)
  if (length(bad_stat)) stop("Unknown stat value(s): ", paste(bad_stat, collapse = ", "))
  if (length(bad_an)) stop("Unknown analyte value(s): ", paste(bad_an, collapse = ", "))
  if (length(bad_reg) || anyNA(reg)) stop("Unknown regimen value(s): ", paste(c(bad_reg, if (anyNA(reg)) "NA"), collapse = ", "),
                                          " (expected ", paste(PK_REGIMENS, collapse = "/"), ")")
  if (any(reg == "daily" & obs$stat == "AUCinf")) stop("AUCinf is not defined for a daily-dosing row")
  pred <- numeric(nrow(obs))
  col <- c(NTX = "Cp_NTX", BN = "Cp_BN")
  # Css13h: mean over 1-3 h window (trapezoid rule); defined for any regimen but used mainly for daily (Bested 2023 day 21).
  css13h_fn <- function(time, conc) {
    w <- time >= 1 & time <= 3
    if (sum(w) < 2) return(NA_real_)
    tw <- time[w]; cw <- conc[w]
    sum(diff(tw) * (head(cw, -1) + tail(cw, -1)) / 2) / (tw[length(tw)] - tw[1])
  }
  stat_of <- function(time, conc) c(nca_summary(time, conc), C24 = stats::approx(time, conc, 24)$y,
                                    C72 = stats::approx(time, conc, 72)$y,
                                    "C0.25" = stats::approx(time, conc, 0.25)$y, "C0.5" = stats::approx(time, conc, 0.5)$y,
                                    "C0.75" = stats::approx(time, conc, 0.75)$y, C1 = stats::approx(time, conc, 1)$y,
                                    Css13h = css13h_fn(time, conc))
  for (d in unique(obs$dose_mg)) {
    at_dose <- obs$dose_mg == d
    want <- intersect(PK_REGIMENS, unique(reg[at_dose]))
    n <- if ("daily" %in% want) pk_n_daily_doses(p, tau) else 1L
    if (superposition) {
      # ONE event-free single-dose solve on the union of the shifted grids; it serves the single rows directly and the daily
      # rows by summation over the shifts (the PK ODE is linear).
      tg <- outer(PK_STAT_GRID, tau * (0:(n - 1)), "+")
      s1 <- simulate_pk(d, sort(unique(as.vector(tg))), p)
      ix <- matrix(match(as.vector(tg), s1$time), nrow(tg))
      if (anyNA(ix)) stop("pk_predict_stats: superposition grid lookup failed")
      profiles <- stats::setNames(lapply(want, function(rg) {
        lapply(stats::setNames(nm = intersect(c("NTX", "BN"), obs$analyte[at_dose & reg == rg])), function(an) {
          c1 <- nM_to_ngml(s1[[col[[an]]]], an)
          if (rg == "single") c1[ix[, 1]] else rowSums(matrix(c1[ix], nrow(ix)))
        })
      }), want)
      times <- PK_STAT_GRID
    } else {
      profiles <- stats::setNames(lapply(want, function(rg) {
        s <- if (rg == "single") simulate_pk(d, PK_STAT_GRID, p) else {
          t_last <- tau * (n - 1)
          z <- simulate_pk(d, t_last + PK_STAT_GRID, p, tau = tau, n_doses = n); z <- z[z$time >= t_last, ]
          z$time <- z$time - t_last; z
        }
        lapply(stats::setNames(nm = intersect(c("NTX", "BN"), obs$analyte[at_dose & reg == rg])),
               function(an) nM_to_ngml(s[[col[[an]]]], an))
      }), want)
      times <- PK_STAT_GRID
    }
    nca <- lapply(profiles, function(pr) lapply(pr, function(cc) stat_of(times, cc)))
    rows <- which(at_dose)
    pred[rows] <- vapply(rows, function(k) unname(nca[[reg[k]]][[obs$analyte[k]]][obs$stat[k]]), 0)
  }
  pred
}

PK_FIT_PARS <- c("ka", "fp", "CL", "V1", "Q", "V2", "CLm", "Vm")
# fp0 / floss0 (the one extra parameter of a dose-dependent structure; PK_DOSE_MODELS) are added to the free set in
# 01_fit_pk.R when that structure is fitted; fp itself (the 100 mg value) stays fixed.
PK_CV_FLOOR <- 0.1         # floor on the CV of Cmax / C24 / thalf / AUC summaries (C24 is a concentration at a fixed nominal
                           # time, so it takes the Cmax floor; the curated C24 rows report CVs of 0.22 and 0.28, above it)
PK_TMAX_CV_FLOOR <- 0.25   # observed Tmax sits on the sampling grid (e.g. 1 h, 2 h): coarse, so a wider floor

# Transformed (optimisation) scale: logit for the fractions fp, fp0 and floss0, log for everything else.
PK_LOGIT_PARS <- c("fp", "fp0", "floss0")
pk_to_theta <- function(p, pars) vapply(pars, function(nm) if (nm %in% PK_LOGIT_PARS) logit(p[[nm]]) else log(p[[nm]]), 0)
pk_from_theta <- function(th, start, pars) {
  p <- start
  for (nm in pars) p[[nm]] <- if (nm %in% PK_LOGIT_PARS) inv_logit(th[[nm]]) else exp(th[[nm]])
  p
}

# Weights (R4 fix round 1, review C2). Each row's CV is, in order of precedence:
#   sd / mean     where the source reports a dispersion. Convention (unchanged since R1): `sd` is the BETWEEN-SUBJECT SD
#                 (never an SE), so the weight does not grow with n (n is ignored); the Verebey 1976 rows are mean +/- SD
#                 of n = 4, and the Bested 2023 day-1 rows carry the SD digitised from the figure;
#   cv_assumed    (optional column) a documented assumed CV where none is reported (Bested 2023 day-21 medians: 0.8);
#   the floor     PK_CV_FLOOR, or PK_TMAX_CV_FLOOR for Tmax (Verebey Tmax rows), applied to every row as a lower bound
#                 (it binds only on rows without a reported or assumed CV; every Bested CV is far above it).
# obj_weight (optional column, default 1) multiplies a row's contribution: the Bested day-1 profile is 4 repeated
# measures on the same subjects per analyte, so each of its rows carries 1/4 and the profile counts as about one
# observation per analyte. The objective sum(obj_weight * ((log pred - log obs) / CV)^2) is -2 log-likelihood (up to a
# constant) with CV as the SD of the log value, so likelihood-ratio thresholds (chi-square) apply to WSS differences.
pk_obs_cv <- function(obs, cv_floor = PK_CV_FLOOR, tmax_cv_floor = PK_TMAX_CV_FLOOR) {
  cv <- obs$sd / obs$mean
  if (!is.null(obs$cv_assumed)) cv <- ifelse(is.finite(cv), cv, obs$cv_assumed)
  cv[!is.finite(cv)] <- 0
  pmax(cv, ifelse(obs$stat == "Tmax", tmax_cv_floor, cv_floor))
}
pk_obs_weight <- function(obs) {
  if (is.null(obs$obj_weight)) return(rep(1, nrow(obs)))
  ifelse(is.finite(obs$obj_weight), obs$obj_weight, 1)
}
pk_row_wss <- function(pred, obs, cv = pk_obs_cv(obs), w = pk_obs_weight(obs)) w * ((log(pred) - log(obs$mean)) / cv)^2

# Weighted least squares on log summary statistics. Tmax is IN the objective (review C3: it is the only datum that
# identifies ka); its CV is floored at PK_TMAX_CV_FLOOR. Parameters not in `pars` are held at their `start` values.
fit_pk <- function(obs, start = pk_default(), pars = PK_FIT_PARS, cv_floor = PK_CV_FLOOR, tmax_cv_floor = PK_TMAX_CV_FLOOR) {
  cv <- pk_obs_cv(obs, cv_floor, tmax_cv_floor); w <- pk_obs_weight(obs)
  th0 <- pk_to_theta(start, pars)
  obj <- function(th) {
    names(th) <- pars
    pred <- tryCatch(pk_predict_stats(pk_from_theta(th, start, pars), obs), error = function(e) rep(NA_real_, nrow(obs)))
    if (any(!is.finite(pred) | pred <= 0)) return(1e10)
    sum(pk_row_wss(pred, obs, cv, w))
  }
  if (length(pars) == 1) {
    o <- stats::optimize(obj, th0 + c(-6, 6), tol = 1e-10)
    o2 <- list(par = stats::setNames(o$minimum, pars), value = o$objective, convergence = 0L)
  } else {
    o1 <- stats::optim(th0, obj, method = "Nelder-Mead", control = list(maxit = 6000, reltol = 1e-12))
    o2 <- stats::optim(o1$par, obj, method = "BFGS", control = list(maxit = 1000, reltol = 1e-14))
  }
  th <- stats::setNames(o2$par, pars)
  H <- stats::optimHess(th, obj)
  dimnames(H) <- list(pars, pars)
  par <- pk_from_theta(th, start, pars)
  list(par = par, theta = th, value = o2$value, convergence = o2$convergence,
       hessian = H, pars = pars, start = start, obs = obs, cv = cv, obj_weight = w,
       row_wss = pk_row_wss(pk_predict_stats(par, obs), obs, cv, w), objective_stats = sort(unique(obs$stat)))
}

# Curvature at the fit. The objective is -2 log-likelihood (Gaussian on log statistics with known CV), so the
# covariance of the transformed parameters is solve(hessian / 2). SEs are returned only if that matrix is positive definite.
pk_hessian_summary <- function(fit) {
  H <- fit$hessian / 2
  ev <- eigen(H, symmetric = TRUE)
  pd <- all(is.finite(ev$values)) && min(ev$values) > 0
  cov <- if (pd) solve(H) else matrix(NA_real_, nrow(H), ncol(H))
  dimnames(cov) <- list(fit$pars, fit$pars)
  list(eigen = ev$values, vectors = ev$vectors, pd = pd, cov = cov,
       se = stats::setNames(if (pd) sqrt(diag(cov)) else rep(NA_real_, length(fit$pars)), fit$pars))
}

# Parent disposition half-lives (h) from the micro-constants: one-compartment gives ln2 / (CL / V1) for both; two-compartment
# gives the alpha (distribution) and beta (terminal) phases, the eigenvalues of the central-peripheral rate matrix.
pk_parent_half_lives <- function(p) {
  k10 <- p$CL / p$V1
  if (!isTRUE(p$Q > 1e-4)) return(c(distribution = log(2) / k10, terminal = log(2) / k10))
  k12 <- p$Q / p$V1; k21 <- p$Q / p$V2
  s <- k10 + k12 + k21; d <- sqrt(s^2 - 4 * k10 * k21)
  c(distribution = log(2) / ((s + d) / 2), terminal = log(2) / ((s - d) / 2))
}

# One- vs two-compartment choice (two-compartment only if the WSS drops by more than lr_threshold, chi-square(2) 95%).
# `free` lists the parameters fitted in the two-compartment model (the one-compartment model drops Q and V2). The
# two-compartment model is fitted from two starts (`start` and the one-compartment optimum plus a small peripheral
# compartment) and the lower WSS is kept, so that a poor local optimum cannot decide the structure.
# One fit described entirely by its job list, so that a job can be shipped to a PSOCK worker: the function is top-level (its
# environment is the global environment, which every worker has after sourcing R/) and every input travels inside `j`.
pk_fit_job <- function(j) fit_pk(j$obs, j$start, pars = j$pars)

# `map_fits(X, FUN)` runs the two independent first-wave fits (one-compartment, and two-compartment from the package
# default); it defaults to sequential lapply, and script 01 passes a cluster-backed mapper. The fits are deterministic and
# use no RNG, so parallelising them is numerically neutral.
pk_select_structure <- function(obs, start, free = PK_FIT_PARS, lr_threshold = stats::qchisq(0.95, 2), map_fits = lapply) {
  force(obs); force(start); force(free)      # never ship unforced promises to a worker
  w1 <- map_fits(list(list(obs = obs, start = modifyList(start, list(Q = 1e-6)), pars = setdiff(free, c("Q", "V2"))),
                      list(obs = obs, start = start, pars = free)),
                 pk_fit_job)
  one <- w1[[1]]; two_a <- w1[[2]]
  two_b <- fit_pk(obs, modifyList(one$par, list(Q = 0.2 * one$par$CL, V2 = 2 * one$par$V1)), pars = free)
  two <- if (two_b$value < two_a$value) two_b else two_a
  two$start_used <- if (identical(two, two_b)) "from one-compartment optimum" else "package default"
  fit <- if (one$value - two$value > lr_threshold) two else one
  fit$structure <- if (identical(fit, two)) "two-compartment" else "one-compartment"
  list(one = one, two = two, fit = fit, lr_threshold = lr_threshold)
}

# Rows predicted outside a factor of `tol` of the observation (stated limitations in the fit output). `exclude` drops rows
# that the caller already reports in full, so they are not listed twice.
pk_misfit_notes <- function(obs, pred, tol = 1.5, exclude = rep(FALSE, nrow(obs))) {
  r <- pred / obs$mean; r[exclude] <- 1
  bad <- which(r > tol | r < 1 / tol)
  if (!length(bad)) return(character(0))
  sprintf("Limitation (stated; no structure added): %s %s %s (%s) observed %.3g, predicted %.3g (%.2f-fold).",
          obs$study[bad], obs$regimen[bad], obs$analyte[bad], obs$stat[bad], obs$mean[bad], pred[bad], r[bad])
}

# ---- R4 fix round 1: structure choice for the dose dependence, identifiability gate, and the PK outputs ----
# Adoption rule (review I7; stated before any fit was run): a one-extra-parameter structure (fp_dose or parent_loss) is
# adopted over the linear model only if (a) its WSS is lower by more than PK_ADOPT_DWSS = qchisq(0.95, 1) = 3.84 (a
# likelihood-ratio test with 1 df, valid because every weight is obj_weight / CV^2 with CV from the reported or a documented
# dispersion), AND (b) its Hessian is positive definite, AND (c) no Verebey 1976 row is predicted worse than by the linear
# model by more than a factor of PK_ADOPT_VEREBEY_FOLD = 1.5 (ratio of fold errors). If both alternatives qualify, the
# larger dWSS wins. `cand` has one row per structure: dose_model, delta_wss (linear WSS minus this WSS), hessian_pd,
# max_verebey_worsening.
PK_ADOPT_DWSS <- stats::qchisq(0.95, 1)
PK_ADOPT_VEREBEY_FOLD <- 1.5
pk_adoption_rule <- function(cand) {
  q <- cand$dose_model != "linear" & is.finite(cand$delta_wss) & cand$delta_wss > PK_ADOPT_DWSS &
    cand$hessian_pd %in% TRUE & is.finite(cand$max_verebey_worsening) & cand$max_verebey_worsening <= PK_ADOPT_VEREBEY_FOLD
  adopted <- if (any(q)) cand$dose_model[q][which.max(cand$delta_wss[q])] else "linear"
  why <- vapply(which(cand$dose_model != "linear"), function(i) sprintf("%s: dWSS %.2f (%s %.2f), Hessian PD %s, worst Verebey row %.2f-fold of linear (limit %.1f) -> %s",
    cand$dose_model[i], cand$delta_wss[i], if (isTRUE(cand$delta_wss[i] > PK_ADOPT_DWSS)) ">" else "<=", PK_ADOPT_DWSS,
    cand$hessian_pd[i], cand$max_verebey_worsening[i], PK_ADOPT_VEREBEY_FOLD, if (q[i]) "qualifies" else "does not qualify"), "")
  list(adopted = adopted, qualifies = q,
       reason = paste0(paste(why, collapse = "; "), ". Adopted: ", adopted,
                       if (sum(q) > 1) " (both qualify: larger dWSS)" else if (!any(q)) " (no alternative qualifies)" else ""))
}
# Largest ratio, over the reference rows, of the alternative fold error to the linear fold error (fold error =
# max(pred / obs, obs / pred)); > 1 means the alternative predicts that row worse.
pk_verebey_worsening <- function(obs_mean, pred_alt, pred_lin)
  max(exp(abs(log(pred_alt / obs_mean))) / exp(abs(log(pred_lin / obs_mean))))

# Identifiability gate (R1; a function since R4 fix round 1 so that the ADOPTED fit, whichever structure, goes through it):
# the Hessian must be positive definite and ka physiologically plausible (absorption half-life >= ka_thalf_min minutes).
# Otherwise the least-identified parameter (or ka, at the floor) is fixed and the fit repeated; notes say why.
PK_KA_THALF_MIN_FLOOR <- 10
pk_identifiability_gate <- function(fit, obs, ka_thalf_min = PK_KA_THALF_MIN_FLOOR, max_iter = 3) {
  ka_floor <- log(2) / (ka_thalf_min / 60)
  keep <- fit[setdiff(names(fit), c("par", "theta", "value", "convergence", "hessian", "pars", "start", "obs", "cv",
                                    "obj_weight", "row_wss", "objective_stats"))]
  notes <- character(0); fixed <- character(0)
  for (iter in seq_len(max_iter)) {
    hs <- pk_hessian_summary(fit)
    ka_ok <- !("ka" %in% fit$pars) || fit$par$ka <= ka_floor
    if (hs$pd && ka_ok) break
    if (!ka_ok) {
      fix_nm <- "ka"; fix_val <- ka_floor
      note <- sprintf("ka fitted to %.3f per h (absorption half-life %.1f min, below the %d min plausibility floor): fixed at %.3f per h, the fastest absorption judged plausible (half-life %d min; plainly-stated assumption, not a sourced value).",
                      fit$par$ka, 60 * log(2) / fit$par$ka, ka_thalf_min, ka_floor, ka_thalf_min)
    } else {
      load <- abs(hs$vectors[, which.min(hs$eigen)]); fix_nm <- fit$pars[which.max(load)]; fix_val <- fit$par[[fix_nm]]
      note <- sprintf("Hessian not positive definite (smallest eigenvalue %.3g); the flattest direction loads most on %s (|loading| %.2f), which is fixed at its current estimate %.4g (not identifiable from these summaries).",
                      min(hs$eigen), fix_nm, max(load), fix_val)
    }
    notes <- c(notes, note); fixed <- c(fixed, fix_nm); message(note)
    st <- modifyList(fit$par, stats::setNames(list(fix_val), fix_nm))
    fit <- c(fit_pk(obs, st, pars = setdiff(fit$pars, fix_nm)), keep)
  }
  list(fit = fit, notes = notes, fixed = fixed)
}

# Parameter table of a fit: fitted parameters with SEs on the transformed scale and back-transformed 95% limits (only if
# the Hessian is positive definite), then every fixed parameter.
pk_param_table <- function(fit, hs = pk_hessian_summary(fit)) {
  tr <- function(nm, x) if (nm %in% PK_LOGIT_PARS) inv_logit(x) else exp(x)
  fx <- unique(as.character(fit$fixed)); all_pars <- c(fit$pars, fx)
  data.frame(parameter = all_pars, estimate = vapply(all_pars, function(nm) as.numeric(fit$par[[nm]]), 0),
             se_transformed = c(unname(hs$se[fit$pars]), rep(NA_real_, length(fx))),
             lo95 = c(mapply(tr, fit$pars, fit$theta[fit$pars] - 1.96 * hs$se[fit$pars]), rep(NA_real_, length(fx))),
             hi95 = c(mapply(tr, fit$pars, fit$theta[fit$pars] + 1.96 * hs$se[fit$pars]), rep(NA_real_, length(fx))),
             transform = ifelse(all_pars %in% PK_LOGIT_PARS, "logit", "log"),
             convergence = fit$convergence,
             source_type = c(rep("fitted (training PK summaries)", length(fit$pars)), rep("fixed (see pk_fit_notes.txt)", length(fx))),
             row.names = NULL, stringsAsFactors = FALSE)
}

# One summary line naming the adopted model and every parameter, written first in pk_fit_notes.txt.
pk_model_line <- function(fit, hs = pk_hessian_summary(fit)) {
  dm <- if (is.null(fit$dose_model)) "linear" else fit$dose_model
  fx <- unique(as.character(fit$fixed))
  sprintf("Adopted PK model: %s, dose model [%s] (%s). Fitted: %s. Fixed: %s. WSS %.4f; Hessian positive definite %s.",
          fit$structure, dm, PK_DOSE_MODELS[[dm]]$label,
          paste(sprintf("%s = %.4g (SE %s on the %s scale)", fit$pars, vapply(fit$pars, function(nm) fit$par[[nm]], 0),
                        ifelse(is.finite(hs$se[fit$pars]), sprintf("%.3g", hs$se[fit$pars]), "NA"),
                        ifelse(fit$pars %in% PK_LOGIT_PARS, "logit", "log")), collapse = "; "),
          paste(sprintf("%s = %.4g", fx, vapply(fx, function(nm) as.numeric(fit$par[[nm]]), 0)), collapse = "; "),
          fit$value, hs$pd)
}

# Every PK output that downstream scripts read is written from ONE fit object here (review C3): pk_fit.rds,
# pk_fit_table.csv, pk_fit_vs_obs.csv, pk_hessian.csv and pk_fit_notes.txt (model line first, then `notes`).
pk_write_outputs <- function(fit, obs, dir, notes = character(0)) {
  hs <- if (is.null(fit$hessian_summary)) pk_hessian_summary(fit) else fit$hessian_summary
  pred <- pk_predict_stats(fit$par, obs)
  vo <- obs; vo$pred <- pred; vo$fold <- pred / obs$mean; vo$in_objective <- TRUE
  vo$cv_used <- pk_obs_cv(obs); vo$obj_weight_used <- pk_obs_weight(obs)
  vo$wss_contribution <- pk_row_wss(pred, obs, vo$cv_used, vo$obj_weight_used)
  hsi <- fit$hessian_summary_initial
  eig <- rbind(if (!is.null(hsi)) data.frame(fit = "initial (before the identifiability gate)", eigenvalue = hsi$eigen, positive_definite = hsi$pd),
               data.frame(fit = "final", eigenvalue = hs$eigen, positive_definite = hs$pd))
  saveRDS(fit, file.path(dir, "pk_fit.rds"))
  utils::write.csv(pk_param_table(fit, hs), file.path(dir, "pk_fit_table.csv"), row.names = FALSE)
  utils::write.csv(vo, file.path(dir, "pk_fit_vs_obs.csv"), row.names = FALSE)
  utils::write.csv(eig, file.path(dir, "pk_hessian.csv"), row.names = FALSE)
  writeLines(c(pk_model_line(fit, hs), notes), file.path(dir, "pk_fit_notes.txt"))
  invisible(vo)
}

# Scripted check (review C3) that every PK output in `dir` describes `fit`: returns the problems found (none = OK).
pk_output_check <- function(dir, fit, obs) {
  bad <- character(0); f <- function(x) file.path(dir, x)
  dm <- if (is.null(fit$dose_model)) "linear" else fit$dose_model
  r <- readRDS(f("pk_fit.rds"))
  if (!identical(r$pars, fit$pars) || !isTRUE(all.equal(unlist(r$par[fit$pars]), unlist(fit$par[fit$pars]))) ||
      !identical(if (is.null(r$dose_model)) "linear" else r$dose_model, dm)) bad <- c(bad, "pk_fit.rds is not the adopted fit")
  tab <- utils::read.csv(f("pk_fit_table.csv"), stringsAsFactors = FALSE)
  want <- unique(c(fit$pars, as.character(fit$fixed)))
  if (!all(want %in% tab$parameter)) {
    bad <- c(bad, sprintf("pk_fit_table.csv lacks %s", paste(setdiff(want, tab$parameter), collapse = ", ")))
  } else if (!isTRUE(all.equal(tab$estimate[match(fit$pars, tab$parameter)], unname(unlist(fit$par[fit$pars])), tolerance = 1e-8)))
    bad <- c(bad, "pk_fit_table.csv estimates differ from the adopted fit")
  vo <- utils::read.csv(f("pk_fit_vs_obs.csv"), stringsAsFactors = FALSE)
  pr <- pk_predict_stats(fit$par, obs)
  if (nrow(vo) != nrow(obs) || !isTRUE(all.equal(vo$pred, pr, tolerance = 1e-6)))
    bad <- c(bad, "pk_fit_vs_obs.csv predictions are not those of the adopted fit")
  he <- utils::read.csv(f("pk_hessian.csv"), stringsAsFactors = FALSE)
  if (sum(he$fit == "final") != length(fit$pars)) bad <- c(bad, "pk_hessian.csv final eigenvalues do not match the adopted fit parameters")
  nl <- readLines(f("pk_fit_notes.txt"), warn = FALSE)
  if (!length(nl) || !grepl(sprintf("dose model [%s]", dm), nl[1], fixed = TRUE) ||
      !all(vapply(fit$pars, function(p) grepl(paste0(p, " = "), nl[1], fixed = TRUE), TRUE)))
    bad <- c(bad, "pk_fit_notes.txt does not open with the adopted model and its parameters")
  ss <- utils::read.csv(f("pk_structure_summary.csv"), stringsAsFactors = FALSE)
  sel <- ss[as.logical(ss$selected) %in% TRUE, , drop = FALSE]
  if (nrow(sel) != 1 || !identical(sel$dose_model, dm) || (!is.null(fit$structure) && !identical(sel$structure, fit$structure)))
    bad <- c(bad, "pk_structure_summary.csv does not mark the adopted structure and dose model")
  sc <- utils::read.csv(f("pk_structure_comparison.csv"), stringsAsFactors = FALSE)
  if (is.null(sc$pred_selected) || !isTRUE(all.equal(sc$pred_selected, pr, tolerance = 1e-6)))
    bad <- c(bad, "pk_structure_comparison.csv pred_selected is not the adopted fit prediction")
  bad
}

# Stable keys of the PK sensitivity fits (file pk_fit_sens_<key>.rds, column `key` of occupancy_sensitivity_J.csv); the
# consumers (scripts 10, 11) select rows by these keys, never by label text (review C1). The D50 keys exist only when a
# dose-dependent structure is adopted.
PK_SENS_KEYS <- c(fp090 = "fp090", fp099 = "fp099", terminal = "terminal", D50_10 = "D50_10", D50_100 = "D50_100")
OCC_SENS_KEYS <- c(primary = "primary", koff_eq = "koff_eq")
