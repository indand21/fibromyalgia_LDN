# Indirect-response turnover per cytokine: dC/dt = kin*(1 - f_c * signal) - kout_c*C, with f_c = 1 (upper bound:
# production fully dependent on the naltrexone-sensitive TLR4/TRIF branch; Wang 2016). Cytokine half-lives are hours,
# so at day 56 each cytokine is at steady state: % change = -100 * signal * (1 - exp(-kout_c * 24 * day)).
# Half-lives (quantity thalf, or thalf_lower + thalf_upper) may be given in minutes ("min") or hours ("h"); any other unit is an error.
cytokine_prediction <- function(hyp, th, prof, cyt_params, day) {
  m <- profile_means(prof, gamma = 1,
                     IC50_TLR4 = if (is.null(th$IC50_TLR4)) Inf else th$IC50_TLR4,
                     Kd_OGFr = if (is.null(th$Kd_OGFr)) Inf else th$Kd_OGFr)
  signal <- switch(hyp, H0 = 0, H1 = 0, H2 = m$tlr4, H3 = th$Emax * m$ogfr, stop("unknown hypothesis ", hyp))
  hl <- cyt_params[cyt_params$quantity %in% c("thalf", "thalf_lower", "thalf_upper"), ]
  bad <- setdiff(unique(hl$unit), c("min", "h"))
  if (length(bad)) stop("unsupported half-life unit(s): ", paste(bad, collapse = ", "), " (use 'min' or 'h')")
  hl$hours <- hl$value / ifelse(hl$unit == "min", 60, 1)
  # Prefer a point estimate (thalf); otherwise use the arithmetic mean of thalf_lower and thalf_upper.
  rows <- lapply(unique(hl$cytokine), function(cy) {
    d <- hl[hl$cytokine == cy, ]
    if (any(d$quantity == "thalf")) return(data.frame(cytokine = cy, hours = d$hours[d$quantity == "thalf"][1], basis = "thalf"))
    if (all(c("thalf_lower", "thalf_upper") %in% d$quantity))
      return(data.frame(cytokine = cy, hours = mean(c(d$hours[d$quantity == "thalf_lower"][1], d$hours[d$quantity == "thalf_upper"][1])),
                        basis = "mean of bounds"))
    NULL
  })
  hh <- do.call(rbind, rows)
  if (is.null(hh) || nrow(hh) == 0) stop("no usable half-life rows (need quantity 'thalf', or both 'thalf_lower' and 'thalf_upper')")
  kout <- log(2) / hh$hours   # 1/h
  data.frame(cytokine = hh$cytokine, pct_change_pred = -100 * signal * (1 - exp(-kout * 24 * day)), halflife_basis = hh$basis)
}
