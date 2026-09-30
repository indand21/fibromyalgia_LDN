binding_params <- function(Kd_NTX, Kd_BN, koff_NTX, koff_BN, keq) {
  list(Kd_NTX = Kd_NTX, Kd_BN = Kd_BN, koff_NTX = koff_NTX, koff_BN = koff_BN,
       kon_NTX = koff_NTX / Kd_NTX, kon_BN = koff_BN / Kd_BN, keq = keq)
}

# Brain unbound conc equilibrates with plasma unbound conc (rate keq, 1/h); competitive MOR binding.
occupancy_ode <- function(t, y, p) {
  free <- 1 - y[["R_NTX"]] - y[["R_BN"]]
  list(c(p$keq * (p$f_ntx(t) - y[["Cb_NTX"]]),
         p$keq * (p$f_bn(t) - y[["Cb_BN"]]),
         p$kon_NTX * y[["Cb_NTX"]] * free - p$koff_NTX * y[["R_NTX"]],
         p$kon_BN * y[["Cb_BN"]] * free - p$koff_BN * y[["R_BN"]]))
}

simulate_occupancy <- function(pk_df, bp) {
  p <- c(bp, list(f_ntx = stats::approxfun(pk_df$time, pk_df$Cu_NTX, rule = 2),
                  f_bn = stats::approxfun(pk_df$time, pk_df$Cu_BN, rule = 2)))
  out <- as.data.frame(deSolve::lsoda(c(Cb_NTX = 0, Cb_BN = 0, R_NTX = 0, R_BN = 0), pk_df$time,
                                      occupancy_ode, p, rtol = 1e-8, atol = 1e-12, hmax = 0.25))
  out$RO <- out$R_NTX + out$R_BN
  out
}

equilibrium_occupancy <- function(c_ntx, c_bn, Kd_NTX, Kd_BN) {
  a <- c_ntx / Kd_NTX; b <- c_bn / Kd_BN
  (a + b) / (1 + a + b)
}
