# Tables 1-4 from result CSVs (no held-out test-data file access; held-out results only via output/results).
options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
res <- "output/results"
dir.create(res, recursive = TRUE, showWarnings = FALSE)

need <- function(...) {
  f <- c(...); miss <- f[!file.exists(f)]
  if (length(miss)) stop("Table script cannot run; missing input file(s): ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(f)
}
rd <- function(f) { need(f); utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }
need_cols <- function(d, cols, f) {
  m <- setdiff(cols, names(d))
  if (length(m)) stop(f, " lacks required column(s): ", paste(m, collapse = ", "), call. = FALSE)
}
# resolve a column by regex, stop if absent
pick_col <- function(d, pattern, what, f) {
  hit <- grep(pattern, names(d), value = TRUE)
  if (!length(hit)) stop(f, ": no ", what, " column found among: ", paste(names(d), collapse = ", "), call. = FALSE)
  hit[1]
}

need("data/module_screen.csv", "data/params/binding_constants.csv",
     file.path(res, c("binding_stage.txt", "pk_fit_table.csv", "training_fit.csv", "pet_metrics.csv",
                      "test_predictions.csv", "test_elpd.csv", "cytokine_test.csv", "virtual_summary.csv")))

# ---- Table 1 ----
t1 <- rd("data/module_screen.csv")
write.csv(t1, file.path(res, "table1_module_screen.csv"), row.names = FALSE)

# ---- Table 2 ----
stage <- readLines(file.path(res, "binding_stage.txt"))[1]
bp_file <- file.path(res, sprintf("binding_params_%s.rds", stage)); need(bp_file)
bp <- readRDS(bp_file)
bc <- read_dataset("data/params/binding_constants.csv")
# Select on `analyte` (NTX or BN); `species` holds the organism (human/mouse). Several rows may match
# (e.g. three human MOR Ki rows for NTX); they are pooled by geometric mean in the scripts, and all
# contributing sources are listed here together with their organism.
src_of <- function(cmp, target, quantity) {
  r <- bc$analyte == cmp & bc$target == target & bc$quantity == quantity
  if (!any(r)) stop("binding_constants.csv has no row for ", cmp, " ", target, " ", quantity, call. = FALSE)
  hum <- r & bc$species == "human"
  use <- if (any(hum)) hum else r
  txt <- paste(unique(paste0(bc$source[use], " [organism: ", bc$species[use], "]")), collapse = "; ")
  if (sum(use) > 1) txt <- paste0("geometric mean of ", sum(use), " rows: ", txt)
  if (!any(hum)) txt <- paste0("non-human fallback (no human row): ", txt)
  txt
}
a2 <- identical(stage, "A2")
kd_type <- if (a2) "PET-calibrated (stage A2)" else "in vitro"
# A2 is the expected path if the bottom-up PET prediction fails; binding_params() does not return koff_source,
# so take it from the A2 object when present, otherwise from the stage-A object.
bpA_file <- file.path(res, "binding_params_A.rds"); need(bpA_file)
bpA <- readRDS(bpA_file)
koff_source <- if (!is.null(bp$koff_source)) bp$koff_source else bpA$koff_source
kp <- rd(file.path(res, "koff_provenance.csv"))   # written by script 02: which koff values are assumed
koff_ntx_assumed <- kp$koff_assumed[kp$analyte == "NTX"]; koff_bn_assumed <- kp$koff_assumed[kp$analyte == "BN"]
koff_ntx_src <- kp$source[kp$analyte == "NTX"]; koff_bn_src <- kp$source[kp$analyte == "BN"]
bind <- data.frame(
  parameter = c("Kd_NTX", "Kd_BN", "koff_NTX", "koff_BN", "keq"),
  estimate = c(bp$Kd_NTX, bp$Kd_BN, bp$koff_NTX, bp$koff_BN, bp$keq),
  lo95 = NA_real_, hi95 = NA_real_,
  source_type = c(kd_type, kd_type, if (koff_ntx_assumed) "assumed (justified)" else "in vitro",
                  if (koff_bn_assumed) "assumed (justified)" else "in vitro", "assumed (justified)"),
  source = c(src_of("NTX", "MOR", "Ki"), src_of("BN", "MOR", "Ki"),
             if (koff_ntx_assumed) koff_ntx_src else src_of("NTX", "MOR", "koff"),
             if (koff_bn_assumed) koff_bn_src else src_of("BN", "MOR", "koff"),
             "Fixed at 2 per hour on the basis of rapid passive blood-brain barrier permeation; not varied (no sensitivity analysis was run)"),
  stringsAsFactors = FALSE)
bpA2_file <- file.path(res, "binding_params_A2.rds")   # A2 is fitted whenever stage A fails the PET test, but adopted only if it improves held-out RMSE
if (file.exists(bpA2_file)) {
  bpA2 <- readRDS(bpA2_file)
  bind <- rbind(bind, data.frame(parameter = "kd_scale", estimate = bpA2$kd_scale, lo95 = NA_real_, hi95 = NA_real_,
    source_type = "PET-calibrated (stage A2)",
    source = paste0("Kd scale fitted to Rabiner 2011 PET occupancy only (all other PET data held out); ",
                    if (a2) "adopted for downstream analyses" else "NOT adopted (worse than stage A on the held-out rows; see stage_decision.csv)"),
    stringsAsFactors = FALSE))
}

pk <- rd(file.path(res, "pk_fit_table.csv"))
need_cols(pk, c("parameter", "estimate", "source_type"), "pk_fit_table.csv")
pk$source <- "This study: fitted to data/training/pk_summary.csv"
pk$lo95 <- NA_real_; pk$hi95 <- NA_real_

tf <- rd(file.path(res, "training_fit.csv"))
need_cols(tf, c("hypothesis", "variable"), "training_fit.csv")
qmed <- pick_col(tf, "^(X?50\\.?%?|q50|median)$", "posterior median", "training_fit.csv")
qlo  <- pick_col(tf, "^(X?2\\.5\\.?%?|q2\\.5|q025)$", "posterior 2.5% quantile", "training_fit.csv")
qhi  <- pick_col(tf, "^(X?97\\.5\\.?%?|q97\\.5|q975)$", "posterior 97.5% quantile", "training_fit.csv")
post <- data.frame(parameter = paste(tf$hypothesis, tf$variable), estimate = tf[[qmed]], lo95 = tf[[qlo]], hi95 = tf[[qhi]],
                   source_type = "fitted (training)",
                   source = "Posterior median and 95% CrI; priors in data/priors/placebo_priors.csv",
                   stringsAsFactors = FALSE)
cols <- c("parameter", "estimate", "lo95", "hi95", "source_type", "source")
t2 <- rbind(pk[, cols], bind[, cols], post[, cols])
nosrc <- is.na(t2$source) | trimws(t2$source) == ""
if (any(nosrc)) stop("Table 2 has parameters without a source: ", paste(t2$parameter[nosrc], collapse = ", "), call. = FALSE)
write.csv(t2, file.path(res, "table2_parameters.csv"), row.names = FALSE)

# ---- Table 3: one block per held-out dataset, long format ----
long <- function(block, hypothesis, item, observed = NA, pred_median = NA, pred_lo90 = NA, pred_hi90 = NA,
                 inside90 = NA, metric = NA_character_, value = NA_real_, footnote = NA_character_)
  data.frame(block = block, hypothesis = hypothesis, item = item, observed = observed, pred_median = pred_median,
             pred_lo90 = pred_lo90, pred_hi90 = pred_hi90, inside90 = inside90, metric = metric, value = value,
             footnote = footnote, stringsAsFactors = FALSE)

pm <- rd(file.path(res, "pet_metrics.csv")); need_cols(pm, c("stage", "n", "RMSE", "MAE", "n_within_15"), "pet_metrics.csv")
# stage A and A2 (if present); metrics are hypothesis-independent
pm <- pm[grepl("^A", pm$stage), , drop = FALSE]
if (!nrow(pm)) stop("pet_metrics.csv has no stage A rows", call. = FALSE)
met <- c("n", "RMSE", "MAE", "n_within_15")
b_pet <- do.call(rbind, lapply(seq_len(nrow(pm)), function(i)
  long("PET occupancy (Lee 1988, Rabiner 2011)", "all (hypothesis-independent)", paste("stage", pm$stage[i]),
       metric = met, value = as.numeric(pm[i, met]))))

tp <- rd(file.path(res, "test_predictions.csv"))
need_cols(tp, c("hypothesis", "trial", "endpoint", "day", "observed", "pred_median", "pred_lo90", "pred_hi90", "lpd", "inside90"), "test_predictions.csv")
dd <- tp[tp$endpoint == "diff_nrs", , drop = FALSE]   # unadjusted difference in change (primary + secondary days); adjusted-level rows are sensitivity only
if (!nrow(dd)) stop("test_predictions.csv has no diff_nrs rows", call. = FALSE)
b_tr <- long("Held-out trials (FINAL, INNOVA): diff in NRS change", dd$hypothesis, paste0(dd$trial, " day ", dd$day, if ("role" %in% names(dd)) paste0(" (", dd$role, ")") else ""),
             dd$observed, dd$pred_median, dd$pred_lo90, dd$pred_hi90, dd$inside90, metric = "lpd", value = dd$lpd)
el <- rd(file.path(res, "test_elpd.csv"))
need_cols(el, c("hypothesis", "elpd_primary", "primary_endpoints", "elpd_all_days_secondary_nonindependent"), "test_elpd.csv")
# Pairwise elpd contrasts on the PRIMARY observations (loo-compare style): diff = sum of per-observation lpd differences,
# SE = sd(per-observation differences) * sqrt(n). With n = 2 the SE is very poorly estimated.
pr <- tp[tp$endpoint == "diff_nrs" & tp$role == "primary", c("hypothesis", "trial", "day", "lpd")]
obs_key <- unique(pr[c("trial", "day")]); n_obs <- nrow(obs_key)
lp <- function(h) { x <- pr[pr$hypothesis == h, ]; x$lpd[match(paste(obs_key$trial, obs_key$day), paste(x$trial, x$day))] }
hh <- sort(unique(pr$hypothesis)); contr <- list()
for (i in seq_along(hh)) for (j in seq_along(hh)) if (i < j) {
  dv <- lp(hh[i]) - lp(hh[j]); est <- sum(dv); se <- stats::sd(dv) * sqrt(n_obs)
  contr[[length(contr) + 1]] <- data.frame(hypothesis_a = hh[i], hypothesis_b = hh[j], elpd_diff_a_minus_b = est,
    se_diff = se, n_obs = n_obs, ratio_diff_over_se = est / se,
    note = "PRIMARY endpoints only (FINAL d84, INNOVA d90). With n = 2 the SE is poorly estimated; the comparison is indicative only, not decisive.")
}
contr <- do.call(rbind, contr)
write.csv(contr, file.path(res, "test_elpd_contrasts.csv"), row.names = FALSE)
b_ct2 <- long("Held-out trials: pairwise PRIMARY elpd contrasts (a minus b; n = 2, SE poorly estimated, indicative only)",
              contr$hypothesis_a, paste(contr$hypothesis_a, "vs", contr$hypothesis_b),
              metric = "elpd_diff_a_minus_b", value = contr$elpd_diff_a_minus_b,
              footnote = sprintf("SE %.3f; diff/SE %.2f; n_obs %d. %s", contr$se_diff, contr$ratio_diff_over_se, contr$n_obs, contr$note))
b_el <- rbind(
  long(paste0("Held-out trials: PRIMARY test log predictive density, one endpoint per trial (", el$primary_endpoints[1], ")"),
       el$hypothesis, "primary endpoints", metric = "elpd_primary", value = el$elpd_primary),
  long("Held-out trials: SECONDARY, non-independent sum over all scored days (INNOVA days correlated)",
       el$hypothesis, "all days", metric = "elpd_all_days_secondary_nonindependent",
       value = el$elpd_all_days_secondary_nonindependent))

ct <- rd(file.path(res, "cytokine_test.csv"))
need_cols(ct, c("cytokine", "pred_median", "pred_lo90", "pred_hi90", "hypothesis", "pct_change", "predicted"), "cytokine_test.csv")
ct$predicted <- as.logical(ct$predicted)
cp <- ct[ct$predicted %in% TRUE, , drop = FALSE]
if (!nrow(cp)) stop("cytokine_test.csv has no predicted cytokines", call. = FALSE)
unpred <- sort(setdiff(unique(ct$cytokine), unique(cp$cytokine)))
fn <- if (length(unpred)) paste("No prediction (no curated half-life):", paste(unpred, collapse = ", ")) else NA_character_
b_ct <- long("Cytokines (Parkitny 2017), % change", cp$hypothesis, cp$cytokine, cp$pct_change, cp$pred_median, cp$pred_lo90, cp$pred_hi90,
             cp$pct_change >= cp$pred_lo90 & cp$pct_change <= cp$pred_hi90, footnote = fn)
t3 <- rbind(b_pet, b_tr, b_el, b_ct2, b_ct)
write.csv(t3, file.path(res, "table3_heldout.csv"), row.names = FALSE)

# ---- Table 4 (omega multiplier 1) and supplement (multiplier 2) ----
vs <- rd(file.path(res, "virtual_summary.csv")); need_cols(vs, "omega_mult", "virtual_summary.csv")
if (!any(vs$omega_mult == 1)) stop("virtual_summary.csv has no omega_mult == 1 rows", call. = FALSE)
write.csv(vs[vs$omega_mult == 1, ], file.path(res, "table4_virtual.csv"), row.names = FALSE)
if (any(vs$omega_mult == 2))
  write.csv(vs[vs$omega_mult == 2, ], file.path(res, "tableS_virtual_omega2.csv"), row.names = FALSE)
cat("Tables written to", res, "\n")
