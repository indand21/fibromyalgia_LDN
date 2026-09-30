# Figures 1-6. Reads result CSVs/RDS (plus training data and, post-freeze, the test trials for observed CIs).
options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
res <- "output/results"; figdir <- "output/figures"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

need <- function(...) {
  f <- c(...); miss <- f[!file.exists(f)]
  if (length(miss)) stop("Figure script cannot run; missing input file(s): ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(f)
}
rd <- function(name) { f <- file.path(res, name); need(f); utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) }
need_cols <- function(d, cols, name) {
  m <- setdiff(cols, names(d))
  if (length(m)) stop(name, " lacks required column(s): ", paste(m, collapse = ", "), call. = FALSE)
}

# Okabe-Ito hypothesis palette; observed data always black
pal_h <- c(H0 = "#999999", H1 = "#0072B2", H2 = "#E69F00", H3 = "#009E73")
col_train <- "#56B4E9"; col_test <- "#D55E00"
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
                                      legend.position = "bottom")
theme_set(th)
scale_h <- scale_colour_manual(values = pal_h, name = "Hypothesis")

save_fig <- function(p, name, w = 180, h = 120) {
  ggplot2::ggsave(file.path(figdir, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300)
  ggplot2::ggsave(file.path(figdir, paste0(name, ".tiff")), p, width = w, height = h, units = "mm", dpi = 300, compression = "lzw")
}

# ============ Figure 1: schematic (ggplot geoms only, no external images) ============
box <- function(id, x0, x1, y0, y1, label, fill) data.frame(id = id, x0 = x0, x1 = x1, y0 = y0, y1 = y1, label = label, fill = fill)
model_fill <- "grey95"
bx <- rbind(
  box("pk",   0.0, 2.4, 3.6, 5.4, "Oral PK\nNTX + 6-beta-naltrexol\n(2-cpt, first-pass)", model_fill),
  box("bind", 3.0, 5.4, 3.6, 5.4, "Brain binding\nunbound conc. -> MOR\noccupancy (competitive)", model_fill),
  box("H0",   6.0, 8.4, 6.0, 7.0, "H0  placebo only", model_fill),
  box("H1",   6.0, 8.4, 4.8, 5.8, "H1  MOR occupancy", model_fill),
  box("H2",   6.0, 8.4, 3.6, 4.6, "H2  TLR4 / cytokines", model_fill),
  box("H3",   6.0, 8.4, 2.4, 3.4, "H3  OGFr", model_fill),
  box("sens", 9.0, 11.2, 3.6, 5.4, "Sensitization\n+ placebo response\n(trial-level Pmax)", model_fill),
  box("pain", 11.8, 13.8, 3.6, 5.4, "Pain NRS\n(change, diff,\nresponders)", model_fill),
  box("tr_pk",  0.0, 2.4, 6.6, 7.6, "Training: PK summaries", col_train),
  box("tr_eff", 6.0, 8.4, 8.0, 9.0, "Training: Younger 2009, 2013", col_train),
  box("te_pet", 3.0, 5.4, 1.4, 2.4, "Test: PET occupancy\n(Lee, Rabiner)", col_test),
  box("te_tr",  11.8, 13.8, 1.4, 2.4, "Test: FINAL, INNOVA", col_test),
  box("te_cyt", 6.0, 8.4, 0.8, 1.8, "Test: cytokines (Parkitny)", col_test))
bx$xm <- (bx$x0 + bx$x1) / 2; bx$ym <- (bx$y0 + bx$y1) / 2
ar <- function(x, y, xend, yend) data.frame(x = x, y = y, xend = xend, yend = yend)
arrows <- rbind(
  ar(2.4, 4.5, 3.0, 4.5), ar(5.4, 4.5, 6.0, 5.3), ar(5.4, 4.5, 6.0, 4.1), ar(5.4, 4.5, 6.0, 2.9),
  ar(8.4, 6.5, 9.0, 5.1), ar(8.4, 5.3, 9.0, 4.8), ar(8.4, 4.1, 9.0, 4.3), ar(8.4, 2.9, 9.0, 3.9),
  ar(11.2, 4.5, 11.8, 4.5),
  ar(1.2, 6.6, 1.2, 5.4), ar(7.2, 8.0, 7.2, 7.0),
  ar(4.2, 2.4, 4.2, 3.6), ar(12.8, 2.4, 12.8, 3.6), ar(7.2, 1.8, 7.2, 2.4))
p1 <- ggplot() +
  geom_rect(data = bx, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill), colour = "black", linewidth = 0.3) +
  geom_text(data = bx, aes(x = xm, y = ym, label = label), size = 2.6, lineheight = 0.9) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = 0.4) +
  annotate("rect", xmin = 9.0, xmax = 9.4, ymin = 8.3, ymax = 8.7, fill = col_train, colour = "black", linewidth = 0.3) +
  annotate("text", x = 9.55, y = 8.5, label = "Training input (fitted, then frozen)", hjust = 0, size = 2.6) +
  annotate("rect", xmin = 9.0, xmax = 9.4, ymin = 7.5, ymax = 7.9, fill = col_test, colour = "black", linewidth = 0.3) +
  annotate("text", x = 9.55, y = 7.7, label = "Held-out test input (never fitted)", hjust = 0, size = 2.6) +
  scale_fill_identity() + coord_cartesian(xlim = c(-0.2, 14.2), ylim = c(0.5, 9.3), expand = FALSE) + theme_void()
save_fig(p1, "fig1_schematic", 180, 110)

# ============ Figure 2: PK and PET ============
pkv <- rd("pk_fit_vs_obs.csv"); need_cols(pkv, c("mean", "pred"), "pk_fit_vs_obs.csv")
pkv <- pkv[is.finite(pkv$mean) & is.finite(pkv$pred) & pkv$mean > 0 & pkv$pred > 0, ]
aes_pk <- if (all(c("analyte", "stat") %in% names(pkv))) aes(mean, pred, colour = analyte, shape = stat) else aes(mean, pred)
lim <- range(c(pkv$mean, pkv$pred))
pa <- ggplot(pkv, aes_pk) + geom_abline(slope = 1, intercept = 0, colour = "grey40") + geom_point(size = 1.8) +
  (if (all(c("analyte", "stat") %in% names(pkv))) scale_shape_manual(values = rep(c(16, 17, 15, 3, 4, 8, 18, 7, 9, 10), length.out = length(unique(pkv$stat)))) else NULL) +
  scale_x_log10(limits = lim) + scale_y_log10(limits = lim) + coord_equal() +
  labs(x = "Observed summary statistic", y = "Predicted", title = "A  PK fit (training)") + theme(legend.box = "vertical")

# Stages: the ADOPTED stage (binding_stage.txt, chosen by the held-out rule in script 03) and the alternative stage, both
# drawn so that a reader can see the failure of each. Panel B plots the quantity that was SCORED: predicted occupancy for the
# actual regimen and time of each Rabiner observation (single dose, scan time from pet_test.csv), computed with the same
# predict_pet() as script 03, against a dose grid; it is NOT the steady-state peak.
need(file.path(res, "binding_stage.txt")); stage <- readLines(file.path(res, "binding_stage.txt"))[1]
alt <- setdiff(c("A", "A2"), stage)
lab_of <- function(sg) sprintf("Stage %s (%s)", sg, if (sg == stage) "adopted" else "not adopted")
stage_files <- c(A = file.path(res, "binding_params_A.rds"), A2 = file.path(res, "binding_params_A2.rds"))
need(stage_files[[stage]]); stages_shown <- c(stage, if (file.exists(stage_files[[alt]])) alt)
labs_shown <- vapply(stages_shown, lab_of, ""); col_vals <- stats::setNames(c(pal_h[["H1"]], "#D55E00")[seq_along(stages_shown)], labs_shown)
lt_vals <- stats::setNames(c("solid", "dashed")[seq_along(stages_shown)], labs_shown)
pet <- rd("pet_test.csv"); need_cols(pet, c("study", "dose_mg", "regimen", "time_h", "occupancy_pct", "stage", "predicted", "reference"), "pet_test.csv")
petA <- pet[grepl("^A \\(", pet$stage), ]
if (!nrow(petA)) stop("pet_test.csv has no stage A rows (stage label 'A (bottom-up)')", call. = FALSE)
pk <- readRDS(file.path(res, "pk_fit.rds"))
rab <- petA[grepl("Rabiner", petA$study) & petA$regimen == "single", ]
if (!nrow(rab) || length(unique(rab$time_h)) != 1) stop("expected Rabiner single-dose rows sharing one scan time in pet_test.csv", call. = FALSE)
grid <- data.frame(dose_mg = sort(unique(c(seq(1, 50, length.out = 40), rab$dose_mg))), regimen = "single", time_h = rab$time_h[1], reference = "absolute")
curveB <- do.call(rbind, lapply(stages_shown, function(sg)
  data.frame(dose_mg = grid$dose_mg, occ = predict_pet(readRDS(stage_files[[sg]]), grid, pk$par), curve = lab_of(sg))))
# the row-level predictions that were actually scored in script 03 (hollow markers on the curves)
scoredB <- do.call(rbind, lapply(stages_shown, function(sg) {
  r <- pet[grepl("Rabiner", pet$study) & pet$regimen == "single" & grepl(if (sg == "A") "^A \\(" else "^A2 ", pet$stage), ]
  data.frame(dose_mg = r$dose_mg, occ = r$predicted, curve = lab_of(sg)) }))
pb <- ggplot() + geom_line(data = curveB, aes(dose_mg, occ, colour = curve, linetype = curve)) +
  geom_point(data = scoredB, aes(dose_mg, occ, colour = curve), shape = 1, size = 2.2) +
  geom_point(data = rab, aes(dose_mg, occupancy_pct), colour = "black", size = 1.8) +
  scale_colour_manual(values = col_vals, name = NULL) + scale_linetype_manual(values = lt_vals, name = NULL) +
  labs(x = "Dose (mg, single dose)", y = sprintf("MOR occupancy at %g h (%%)", rab$time_h[1]),
       title = "B  Occupancy vs dose") + coord_cartesian(ylim = c(0, 100)) + theme(legend.position = "inside", legend.position.inside = c(0.98, 0.03), legend.justification = c(1, 0), legend.background = element_rect(fill = "white", colour = NA), legend.text = element_text(size = 6), legend.key.size = unit(3.5, "mm"))

# Panel C plots the SCORED quantity. Lee 1988 is percent blockade normalised to that study's own 1 h scan, so the curve is
# 100 * RO(t) / RO(1 h) for the same dose and regimen (same convention as predict_pet(), PET_NORM_REF, script 03).
curve_for <- function(bpf, lab) {
  bp <- readRDS(bpf); s <- simulate_pk(50, seq(0, max(tt) + 1, 0.25), pk$par, 24, 1); o <- simulate_occupancy(s, bp)
  ro <- stats::approx(o$time, o$RO, tt)$y; ro1 <- stats::approx(o$time, o$RO, 1)$y
  data.frame(time_h = tt, occ = 100 * ro / ro1, curve = lab)
}
tt <- seq(1, 168, 1)
cv <- do.call(rbind, lapply(stages_shown, function(sg) curve_for(stage_files[[sg]], lab_of(sg))))
lee <- petA[grepl("Lee", petA$study) & petA$dose_mg == 50 & petA$reference == PET_NORM_REF, ]
if (!nrow(lee)) stop("expected normalised Lee 1988 50 mg rows in pet_test.csv", call. = FALSE)
pc <- ggplot() + geom_line(data = cv, aes(time_h, occ, colour = curve, linetype = curve)) +
  geom_point(data = lee, aes(time_h, occupancy_pct), colour = "black", size = 1.8) +
  scale_colour_manual(values = col_vals, name = NULL) + scale_linetype_manual(values = lt_vals, name = NULL) +
  labs(x = "Time after 50 mg (h)", y = "Occupancy relative to 1 h scan (%)", title = "C  Time course, 50 mg") +
  coord_cartesian(ylim = c(0, 105)) + theme(legend.position = "none")
save_fig(pa + pb + pc + plot_layout(nrow = 1), "fig2_pk_pet", 180, 100)

# ============ Figure 3: training fits (posterior-predicted arm outcomes) ============
arms <- read_dataset("data/training/efficacy_training.csv")
# Only arms in the fitted likelihood (finite se; see script 04) have a trial effect in the posterior.
arms <- arms[is.finite(arms$se), , drop = FALSE]
need(file.path(res, sprintf("exposure_cache_%s.rds", stage)), file.path(res, sprintf("posterior_%s.rds", c("H0", "H1", "H2", "H3"))))
cache <- readRDS(file.path(res, sprintf("exposure_cache_%s.rds", stage)))
pred_rows <- do.call(rbind, lapply(c("H0", "H1", "H2", "H3"), function(h) {
  f <- readRDS(file.path(res, sprintf("posterior_%s.rds", h)))
  idx <- round(seq(1, nrow(f$draws), length.out = 300))
  do.call(rbind, lapply(seq_len(nrow(arms)), function(i) {
    a <- arms[i, ]
    v <- vapply(idx, function(k) { par <- f$draws[k, ]
      sim <- simulate_arm(h, theta_from(par, h, f$fixed), parse_schedule(a$schedule), a$day, cache, pmax_trial(par, a$trial))
      arm_outcome(sim, a$day, a$outcome, a$baseline_nrs) }, 0)
    data.frame(hypothesis = h, trial = a$trial, arm = a$arm, med = stats::median(v),
               lo = unname(stats::quantile(v, 0.05)), hi = unname(stats::quantile(v, 0.95)))
  }))
}))
arms$arm_id <- paste(arms$trial, arms$arm, sep = "\n")
pred_rows$arm_id <- paste(pred_rows$trial, pred_rows$arm, sep = "\n")
arms$se_plot <- ifelse(is.na(arms$se), NA_real_, arms$se)
pd <- position_dodge(width = 0.6)
p3 <- ggplot() +
  geom_pointrange(data = pred_rows, aes(arm_id, med, ymin = lo, ymax = hi, colour = hypothesis), position = pd, size = 0.35) +
  geom_point(data = arms, aes(arm_id, mean), colour = "black", shape = 15, size = 1.8, position = position_nudge(x = -0.42)) +
  geom_linerange(data = arms[is.finite(arms$se), ], aes(arm_id, mean, ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                 colour = "black", position = position_nudge(x = -0.42)) +
  scale_h + labs(x = NULL, y = "Pain reduction (%)", title = "Training fit: observed mean (black square; 95% CI only where SE reported) vs posterior prediction (median, 90% interval)") +
  theme(plot.title = element_text(size = 8))
save_fig(p3, "fig3_training_fits", 180, 100)

# ============ Figure 4: held-out trial predictions (forest) ============
tp <- rd("test_predictions.csv")
need_cols(tp, c("hypothesis", "trial", "endpoint", "day", "observed", "pred_median", "pred_lo90", "pred_hi90", "obs_se"), "test_predictions.csv")
tp <- tp[tp$endpoint != "diff_nrs_adjusted_level", , drop = FALSE]   # sensitivity rows (adjusted level estimand) are not plotted
tp$label <- paste0(ifelse(tp$endpoint == "diff_nrs", "diff (LDN - placebo)", sub("^change_nrs:", "change ", tp$endpoint)), ", d", tp$day)
# observed 95% CI from the held-out trial table (read after the stage-B freeze, via the guarded reader)
to <- read_test_data("data/test/trials_test.csv", "B")
to$label <- paste0(ifelse(to$endpoint == "diff_nrs", "diff (LDN - placebo)", paste0("change ", to$arm)), ", d", to$day)
to <- to[to$endpoint %in% c("diff_nrs", "change_nrs"), c("trial", "label", "lo95", "hi95")]
# Row order (top to bottom): trials FINAL then INNOVA (facet order); within a trial, follow-up day ascending; within a day,
# the LDN - placebo difference, then the LDN arm, then the placebo arm. Hypotheses dodge in the order H0, H1, H2, H3.
row_key <- unique(tp[, c("trial", "label", "day", "endpoint")])
row_key$kind <- ifelse(row_key$endpoint == "diff_nrs", 1, ifelse(grepl("placebo", row_key$endpoint), 3, 2))
row_key <- row_key[order(factor(row_key$trial, c("final", "innova")), row_key$day, row_key$kind), ]
lab_levels <- rev(unique(row_key$label))   # discrete y axis draws the first level at the bottom
tp$label <- factor(tp$label, levels = lab_levels)
tp$hypothesis <- factor(tp$hypothesis, levels = c("H0", "H1", "H2", "H3"))
tp$trial <- factor(tp$trial, levels = c("final", "innova"))
obs4 <- unique(tp[tp$hypothesis == levels(tp$hypothesis)[1], c("trial", "label", "observed")])
obs4 <- merge(obs4, to, by = c("trial", "label"), all.x = TRUE)
# Diff rows: the scored observation may be the unadjusted change difference (INNOVA), so its interval is observed +/- 1.96 x the
# SE used in the scoring (borrowed from the adjusted CI for INNOVA), not the adjusted-level CI from the trial table.
se4 <- unique(tp[is.finite(tp$obs_se) & tp$endpoint == "diff_nrs", c("trial", "label", "obs_se")])
obs4 <- merge(obs4, se4, by = c("trial", "label"), all.x = TRUE)
i4 <- is.finite(obs4$obs_se); obs4$lo95[i4] <- obs4$observed[i4] - 1.96 * obs4$obs_se[i4]; obs4$hi95[i4] <- obs4$observed[i4] + 1.96 * obs4$obs_se[i4]
p4 <- ggplot() +
  geom_pointrange(data = tp, aes(pred_median, label, xmin = pred_lo90, xmax = pred_hi90, colour = hypothesis),
                  position = position_dodge(width = 0.7), size = 0.3) +
  geom_point(data = obs4, aes(observed, label), colour = "black", shape = 15, size = 1.8, position = position_nudge(y = 0.42)) +
  geom_linerange(data = obs4[is.finite(obs4$lo95) & is.finite(obs4$hi95), ], aes(y = label, xmin = lo95, xmax = hi95),
                 colour = "black", position = position_nudge(y = 0.42), orientation = "y") +
  geom_vline(xintercept = 0, linetype = "dotted") + facet_grid(trial ~ ., scales = "free_y", space = "free_y") +
  scale_h + labs(x = "Change in NRS (points)", y = NULL, title = "Held-out trials: prediction (median, 90% interval) vs observed (black square; 95% CI where available)") +
  theme(plot.title = element_text(size = 8))
save_fig(p4, "fig4_heldout", 180, 140)

# ============ Figure 5: cytokines ============
ct <- rd("cytokine_test.csv"); need_cols(ct, c("cytokine", "pred_median", "pred_lo90", "pred_hi90", "hypothesis", "pct_change", "se", "predicted"), "cytokine_test.csv")
ct$predicted <- as.logical(ct$predicted); cp <- ct[ct$predicted %in% TRUE, ]
if (!nrow(cp)) stop("cytokine_test.csv has no predicted cytokines", call. = FALSE)
obs5 <- unique(cp[, c("cytokine", "pct_change", "se")])
p5 <- ggplot() +
  geom_pointrange(data = cp, aes(cytokine, pred_median, ymin = pred_lo90, ymax = pred_hi90, colour = hypothesis),
                  position = position_dodge(width = 0.6), size = 0.35) +
  geom_point(data = obs5, aes(cytokine, pct_change), colour = "black", shape = 15, size = 1.8, position = position_nudge(x = -0.4)) +
  geom_linerange(data = obs5[is.finite(obs5$se), ], aes(cytokine, pct_change, ymin = pct_change - 1.96 * se, ymax = pct_change + 1.96 * se),
                 colour = "black", position = position_nudge(x = -0.4)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_h + labs(x = NULL, y = "Change from baseline (%)",
                 title = "Cytokines: predicted (median, 90% interval) vs observed (black square; 95% CI bars only where SE reported)") +
  theme(plot.title = element_text(size = 8))
save_fig(p5, "fig5_cytokines", 180, 100)

# ============ Figure 6: virtual trials and power ============
vt <- rd("virtual_trials.csv"); need_cols(vt, c("diff", "hypothesis", "trial", "omega_mult"), "virtual_trials.csv")
vs <- rd("virtual_summary.csv"); need_cols(vs, c("hypothesis", "trial", "omega_mult", "obs_diff"), "virtual_summary.csv")
vt <- vt[vt$omega_mult == 1, ]; ob6 <- unique(vs[vs$omega_mult == 1, c("trial", "obs_diff")])
need_cols(vs, "obs_diff_adjusted_level_sensitivity", "virtual_summary.csv")
ob6s <- unique(vs[vs$omega_mult == 1 & is.finite(vs$obs_diff_adjusted_level_sensitivity), c("trial", "obs_diff_adjusted_level_sensitivity")])
p6a <- ggplot(vt, aes(hypothesis, diff, fill = hypothesis)) + geom_violin(colour = NA, alpha = 0.8, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_hline(data = ob6, aes(yintercept = obs_diff), colour = "black", linewidth = 0.6) +
  geom_hline(data = ob6s, aes(yintercept = obs_diff_adjusted_level_sensitivity), colour = "black", linewidth = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~trial) + scale_fill_manual(values = pal_h, guide = "none") +
  labs(x = NULL, y = "Simulated diff in NRS change (LDN - placebo)", title = "A  Simulated treatment difference") + theme(plot.title = element_text(size = 8))
pw <- rd("power_curve.csv")
need_cols(pw, c("hypothesis", "trial", "n_per_arm", "rej_rate_with_model_unc", "rej_rate_sampling_only"), "power_curve.csv")
best <- unique(pw$hypothesis)[1]
pwl <- rbind(data.frame(pw[c("trial", "n_per_arm")], variant = "with model-structure uncertainty", rate = pw$rej_rate_with_model_unc),
             data.frame(pw[c("trial", "n_per_arm")], variant = "sampling variability only", rate = pw$rej_rate_sampling_only))
p6b <- ggplot(pwl, aes(n_per_arm, rate, linetype = variant, shape = trial)) + geom_line(colour = pal_h[[best]]) + geom_point(colour = pal_h[[best]], size = 1.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "black", linewidth = 0.4) +
  annotate("text", x = min(pwl$n_per_arm), y = 0.05, label = "nominal level 0.05", hjust = 0, vjust = -0.6, size = 2.4) +
  coord_cartesian(ylim = range(c(pwl$rate, 0.05)) + c(-0.004, 0.004)) +
  labs(x = "Participants per arm", y = "One-sided 5% rejection rate", linetype = NULL, shape = "Design",
       title = paste0("B  Type I error check (", best, ")")) + guides(linetype = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) + theme(legend.box = "vertical", legend.text = element_text(size = 7))
save_fig(p6a + p6b + plot_layout(widths = c(1.4, 1)), "fig6_virtual_trials", 180, 110)
cat("Figures written to", figdir, "\n")
