options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
obs <- read_dataset("data/training/pk_summary.csv")
bc <- read_dataset("data/params/binding_constants.csv")
start <- pk_default()
fu_NTX <- pool_binding(bc, "NTX", "plasma_protein", "fu")$value
fu_BN <- pool_binding(bc, "BN", "plasma_protein", "fu")$value
start$fu_NTX <- fu_NTX
start$fu_BN <- if (is.finite(fu_BN)) fu_BN else start$fu_NTX

two <- fit_pk(obs, start)
one <- fit_pk(obs, modifyList(start, list(Q = 1e-6)), pars = setdiff(PK_FIT_PARS, c("Q", "V2")))
# Heuristic: WSS difference vs chi-square(2) (2 extra parameters: Q, V2)
fit <- if (one$value - two$value > qchisq(0.95, 2)) two else one
fit$structure <- if (identical(fit, two)) "two-compartment" else "one-compartment"
if (fit$convergence != 0) stop("PK fit did not converge (optim convergence code ", fit$convergence, ")")
na_se <- rep(NA_real_, length(fit$pars))
se <- tryCatch(sqrt(diag(solve(fit$hessian / 2))), warning = function(w) na_se, error = function(e) na_se)
if (any(!is.finite(se))) { se <- na_se; warning("Hessian not positive definite: SEs unavailable") }
tab <- data.frame(parameter = fit$pars, estimate = unlist(fit$par[fit$pars]), se_transformed = se,
                  convergence = fit$convergence, source_type = "fitted (training PK summaries)")
obs$pred <- pk_predict_stats(fit$par, obs)
saveRDS(fit, "output/results/pk_fit.rds")
write.csv(tab, "output/results/pk_fit_table.csv", row.names = FALSE)
write.csv(obs, "output/results/pk_fit_vs_obs.csv", row.names = FALSE)
cat(fit$structure, "WSS =", fit$value, "\n")
