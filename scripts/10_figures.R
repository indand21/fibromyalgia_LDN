# Figures (R3b). Reads result files only (no held-out test-data file: observed values and their SEs come from the
# result tables). Every figure is written as PNG and as LZW-compressed TIFF at 300 dpi; caption notes, including the
# estimator statements, go to output/results/figure_notes.csv for the manuscript.
#   fig1_schematic        model structure; training, calibration and held-out inputs
#   fig2_pk_pet           A PK fit (all training rows); B PET dose-response; C Lee time course + Weerts (stages A, A2, J)
#   fig3_training_fits    posterior-predicted training arms
#   fig4_heldout          held-out differences: observed (95% CI), predictive and true-difference intervals
#   fig5_cytokines        cytokine predictions vs Parkitny
#   fig6_virtual_trials   simulated LDN - placebo differences vs observed
#   fig7_occupancy        stage-J steady-state occupancy vs dose with its band; structural ranges at 4.5 mg
#   fig8_power            assurance (H1, H3) and the H0 type-I-error check
#   figS1_sobol           Sobol indices, eligible hypotheses, days 14 and 84
#   figS2_power_sd_sensitivity  primary vs SD-sensitivity calibration (consistency check)
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
notes <- list()
note <- function(figure, panel, text, sources)
  notes[[length(notes) + 1]] <<- data.frame(figure = figure, panel = panel, note = text, source_files = paste(sources, collapse = "; "))

# Okabe-Ito hypothesis palette; observed data always black. Binding stages get their own Okabe-Ito colours.
pal_h <- c(H0 = "#999999", H1 = "#0072B2", H2 = "#E69F00", H3 = "#009E73")
col_train <- "#56B4E9"; col_test <- "#D55E00"; col_calib <- "#F0E442"
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
                                      legend.position = "bottom")
theme_set(th)
scale_h <- scale_colour_manual(values = pal_h, name = "Hypothesis")

save_fig <- function(p, name, w = 180, h = 120) {
  ggplot2::ggsave(file.path(figdir, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300)
  ggplot2::ggsave(file.path(figdir, paste0(name, ".tiff")), p, width = w, height = h, units = "mm", dpi = 300, compression = "lzw")
}

# ---- adopted binding stage (never a fallback to another stage) ----
need(file.path(res, "binding_stage.txt")); stage <- readLines(file.path(res, "binding_stage.txt"))[1]
dec <- rd("stage_decision.csv"); need_cols(dec, "stage_adopted", "stage_decision.csv")
if (!identical(dec$stage_adopted[1], stage))
  stop("binding_stage.txt says ", stage, " but stage_decision.csv says ", dec$stage_adopted[1], call. = FALSE)
STAGES <- c("A", "A2", "J")
stage_files <- stats::setNames(file.path(res, sprintf("binding_params_%s.rds", STAGES)), STAGES)
need(stage_files, file.path(res, sprintf("exposure_cache_%s.rds", stage)), file.path(res, "pk_fit.rds"))
role_of <- function(sg) switch(sg,
  A = "A: bottom-up, in vitro (held-out prediction)",
  A2 = "A2: Kd scaled to Rabiner only (sensitivity)",
  J = "J: Kd and koff calibrated on all PET (calibration, not a validation)")
lab_of <- function(sg) paste0(role_of(sg), if (sg == stage) " [adopted]" else "")
col_stage <- stats::setNames(c("#D55E00", "#CC79A7", "#0072B2"), vapply(STAGES, lab_of, ""))
lt_stage <- stats::setNames(c("dashed", "dotdash", "solid"), vapply(STAGES, lab_of, ""))
scale_stage <- list(scale_colour_manual(values = col_stage, name = NULL), scale_linetype_manual(values = lt_stage, name = NULL))
pk <- readRDS(file.path(res, "pk_fit.rds"))
bps <- lapply(stage_files, readRDS)

# ============ Figure 1: schematic (ggplot geoms only, no external images) ============
box <- function(id, x0, x1, y0, y1, label, fill) data.frame(id = id, x0 = x0, x1 = x1, y0 = y0, y1 = y1, label = label, fill = fill)
model_fill <- "grey95"
bx <- rbind(
  box("pk",   0.0, 2.4, 3.6, 5.4, sprintf("Oral PK\nNTX + 6-beta-naltrexol\n(%s, first-pass)", sub("-compartment", "-cpt", pk$structure)), model_fill),
  box("bind", 3.0, 5.4, 3.6, 5.4, "Brain binding\nunbound conc. -> MOR\noccupancy (competitive)", model_fill),
  box("H0",   6.0, 8.4, 6.0, 7.0, "H0  placebo only", model_fill),
  box("H1",   6.0, 8.4, 4.8, 5.8, "H1  MOR occupancy", model_fill),
  box("H2",   6.0, 8.4, 3.6, 4.6, "H2  TLR4 / cytokines", model_fill),
  box("H3",   6.0, 8.4, 2.4, 3.4, "H3  OGFr", model_fill),
  box("sens", 9.0, 11.2, 3.6, 5.4, "Sensitisation\n+ placebo response\n(trial-level Pmax)", model_fill),
  box("pain", 11.8, 13.8, 3.6, 5.4, "Pain NRS\n(change, diff,\nresponders)", model_fill),
  box("tr_pk",  0.0, 2.4, 6.6, 7.6, "Training: PK summaries\n(Verebey 100 mg;\nBested 4.5 mg)", col_train),
  box("tr_eff", 6.0, 8.4, 8.0, 9.0, "Training:
Younger 2013
(Younger 2009 excluded)", col_train),
  box("te_pet", 2.7, 5.7, 1.2, 2.6, "PET occupancy\n(Lee, Rabiner, Weerts)", col_calib),
  box("te_tr",  11.8, 13.8, 1.4, 2.4, "Test:
FINAL, INNOVA", col_test),
  box("te_cyt", 6.0, 8.4, 0.8, 1.8, "Test: cytokines
(Parkitny)", col_test))
bx$xm <- (bx$x0 + bx$x1) / 2; bx$ym <- (bx$y0 + bx$y1) / 2
ar <- function(x, y, xend, yend) data.frame(x = x, y = y, xend = xend, yend = yend)
arrows <- rbind(
  ar(2.4, 4.5, 3.0, 4.5), ar(5.4, 4.5, 6.0, 5.3), ar(5.6, 3.25, 6.0, 3.9), ar(5.6, 3.25, 6.0, 3.0),
  ar(8.4, 6.5, 9.0, 5.1), ar(8.4, 5.3, 9.0, 4.8), ar(8.4, 4.1, 9.0, 4.3), ar(8.4, 2.9, 9.0, 3.9),
  ar(11.2, 4.5, 11.8, 4.5),
  ar(1.2, 6.6, 1.2, 5.4),
  ar(4.2, 2.6, 4.2, 3.6), ar(12.8, 2.4, 12.8, 3.6))
h23_path <- rbind(ar(1.2, 3.6, 1.2, 3.25), ar(1.2, 3.25, 5.6, 3.25))
# Training flow (Younger 2013 -> all four hypotheses): a bracket in the training colour on the LEFT of the H boxes.
# Test flow (cytokines -> H2 and H3): dashed, in the test colour, on the RIGHT, so the two flows never share a path
# (R4 fix round 1, review M10). Short arrowheads keep the 0.25-unit stubs from being all head.
yr_path <- rbind(ar(6.0, 8.5, 5.8, 8.5), ar(5.8, 8.5, 5.8, 2.9))
yr_arrows <- rbind(ar(5.8, 6.5, 6.0, 6.5), ar(5.8, 5.3, 6.0, 5.3), ar(5.8, 4.1, 6.0, 4.1), ar(5.8, 2.9, 6.0, 2.9))
cyt_path <- rbind(ar(8.4, 1.3, 8.7, 1.3), ar(8.7, 1.3, 8.7, 3.75))
cyt_arrows <- rbind(ar(7.8, 1.8, 7.8, 2.4), ar(8.7, 3.75, 8.4, 3.75))
leg1 <- data.frame(y = c(8.5, 7.7, 6.9), fill = c(col_train, col_calib, col_test),
  label = c("Training input (fitted, then frozen)",
            sprintf("PET: calibrates binding (stage %s);
held-out test only for stage A", stage),
            "Held-out test input (never fitted)"))
p1 <- ggplot() +
  geom_rect(data = bx, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill), colour = "black", linewidth = 0.3) +
  geom_text(data = bx, aes(x = xm, y = ym, label = label), size = 2.6, lineheight = 0.9) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = 0.4) +
  # H2 (TLR4) and H3 (OGFr) are driven by unbound plasma naltrexone from the PK model, not by MOR occupancy
  geom_segment(data = h23_path, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.4) +
  annotate("text", x = 1.3, y = 3.05, label = "unbound plasma NTX (H2, H3)", hjust = 0, size = 2.1) +
  geom_segment(data = yr_path, aes(x = x, y = y, xend = xend, yend = yend), colour = "#0072B2", linewidth = 0.5) +
  geom_segment(data = yr_arrows, aes(x = x, y = y, xend = xend, yend = yend), colour = "#0072B2", linewidth = 0.5,
               arrow = arrow(length = unit(1.1, "mm"), type = "closed")) +
  geom_segment(data = cyt_path, aes(x = x, y = y, xend = xend, yend = yend), colour = col_test, linewidth = 0.5, linetype = "22") +
  geom_segment(data = cyt_arrows, aes(x = x, y = y, xend = xend, yend = yend), colour = col_test, linewidth = 0.5, linetype = "22",
               arrow = arrow(length = unit(1.3, "mm"), type = "closed")) +
  geom_rect(data = leg1, aes(xmin = 9.0, xmax = 9.4, ymin = y - 0.2, ymax = y + 0.2, fill = fill), colour = "black", linewidth = 0.3) +
  geom_text(data = leg1, aes(x = 9.55, y = y, label = label), hjust = 0, size = 2.4) +
  scale_fill_identity() + coord_cartesian(xlim = c(-0.2, 14.4), ylim = c(0.5, 9.3), expand = FALSE) + theme_void()
save_fig(p1, "fig1_schematic", 180, 110)
note("fig1_schematic", "", sprintf("PET occupancy is calibration data for the adopted binding stage (%s) and a held-out test only for the bottom-up stage A.", stage),
     c("binding_stage.txt", "stage_decision.csv"))

# ============ Figure 2: PK and PET ============
pkv <- rd("pk_fit_vs_obs.csv"); need_cols(pkv, c("analyte", "stat", "regimen", "mean", "pred"), "pk_fit_vs_obs.csv")
pkv <- pkv[is.finite(pkv$mean) & is.finite(pkv$pred) & pkv$mean > 0 & pkv$pred > 0, ]
pkv$analyte <- ifelse(pkv$analyte == "BN", "6-beta-naltrexol", "naltrexone")
need_cols(pkv, c("study", "dose_mg"), "pk_fit_vs_obs.csv")
# Low-dose rows (below the 100 mg reference dose) get a square; once-daily rows a ring. Both legends are computed from the
# rows themselves (review M2), e.g. "square: 4.5 mg (Bested 2023)".
pkv$low_dose <- pkv$dose_mg < FP_DREF_MG
low_lab <- paste(unique(sprintf("%g mg (%s)", pkv$dose_mg[pkv$low_dose], pkv$study[pkv$low_dose])), collapse = "; ")
ring_lab <- paste(unique(sprintf("%g mg", pkv$dose_mg[pkv$regimen == "daily"])), collapse = " and ")
lim <- range(c(pkv$mean, pkv$pred)) * c(0.7, 1.4)
# every PK statistic has a shape (R4 fix round 1 adds the Bested day-1 profile times)
stat_shapes <- c(Cmax = 16, Tmax = 17, thalf = 15, C24 = 18, C72 = 8, AUC = 3, "C0.25" = 6, "C0.5" = 7, "C0.75" = 9, C1 = 10, Css13h = 2)
if (length(setdiff(unique(pkv$stat), names(stat_shapes)))) stop("Fig 2A: no shape for stat(s) ", paste(setdiff(unique(pkv$stat), names(stat_shapes)), collapse = ", "), call. = FALSE)
dm_pk <- pk_dose_model_of(pk$par)
fp_txt <- if (dm_pk == "linear") sprintf("fp %.2f at every dose", pk$par$fp) else
  PK_DOSE_MODELS[[dm_pk]]$label
pa <- ggplot(pkv, aes(mean, pred)) + geom_abline(slope = 1, intercept = 0, colour = "grey40") +
  geom_point(data = pkv[pkv$regimen == "daily", ], colour = "black", shape = 1, size = 3.6, stroke = 0.4) +
  # R4: low-dose Bested rows get an outer square marker to distinguish from 100 mg training rows
  geom_point(data = pkv[pkv$low_dose %in% TRUE, ], colour = "black", shape = 0, size = 3.8, stroke = 0.4) +
  geom_point(aes(colour = analyte, shape = stat), size = 1.9) +
  scale_colour_manual(values = c("naltrexone" = "#0072B2", "6-beta-naltrexol" = "#009E73"), name = NULL) +
  scale_shape_manual(values = stat_shapes[sort(unique(pkv$stat))], name = NULL) +
  scale_x_log10(limits = lim) + scale_y_log10(limits = lim) + coord_equal() +
  guides(colour = guide_legend(order = 1, ncol = 1), shape = guide_legend(order = 2, ncol = 2)) +
  labs(x = "Observed (ng/mL or h)", y = "Predicted", title = sprintf("A  PK fit, %s", sub("-compartment", "-cpt", pk$structure)),
       subtitle = sprintf("ring: once daily, steady state (%s)\nsquare: %s", ring_lab, low_lab)) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1), legend.box = "vertical",
        legend.text = element_text(size = 5.5), legend.key.size = unit(2.4, "mm"), legend.spacing.y = unit(0.5, "mm"),
        legend.background = element_rect(fill = "white", colour = NA), plot.subtitle = element_text(size = 6.5),
        legend.margin = margin(1, 1, 1, 1))
note("fig2_pk_pet", "A", sprintf("PK structure %s, dose model %s (implied systemic bioavailability %.3f at 4.5 mg and %.3f at 100 mg; fp = %.3f at 4.5 mg and %.3f at 100 mg). Every training row, including the 24 h concentrations (C24), the once-daily steady-state rows and the low-dose rows (C0.25-C1 = concentration 15-60 min after the first dose; Css13h = mean over the 1-3 h window at steady state), is plotted against its prediction. Square outline: %s; ring: once daily at steady state (%s). Fit notes: pk_fit_notes.txt; dose-model comparison: pk_dose_model_comparison.csv.",
                                 pk$structure, fp_txt, pk_bioavailability(4.5, pk$par), pk_bioavailability(100, pk$par),
                                 fp_at_dose(4.5, pk$par), fp_at_dose(100, pk$par), low_lab, ring_lab),
     c("pk_fit_vs_obs.csv", "pk_fit.rds", "pk_fit_notes.txt", "pk_dose_model_comparison.csv"))

pet <- rd("pet_test.csv"); need_cols(pet, c("study", "dose_mg", "regimen", "time_h", "occupancy_pct", "sd", "stage", "predicted", "reference"), "pet_test.csv")
stage_rows <- function(sg) pet[grepl(switch(sg, A = "^A \\(", A2 = "^A2 ", J = "^J "), pet$stage), ]
for (sg in STAGES) if (!nrow(stage_rows(sg))) stop("pet_test.csv has no stage ", sg, " rows", call. = FALSE)
obsP <- stage_rows("A")                                            # each observation once
rab <- obsP[grepl("Rabiner", obsP$study) & obsP$regimen == "single", ]
if (!nrow(rab) || length(unique(rab$time_h)) != 1) stop("expected Rabiner single-dose rows sharing one scan time in pet_test.csv", call. = FALSE)
gridB <- data.frame(dose_mg = sort(unique(c(exp(seq(log(0.5), log(60), length.out = 50)), rab$dose_mg))), regimen = "single",
                    time_h = rab$time_h[1], reference = "absolute", n_doses = 1)
curveB <- do.call(rbind, lapply(STAGES, function(sg) data.frame(dose_mg = gridB$dose_mg, occ = predict_pet(bps[[sg]], gridB, pk$par), curve = lab_of(sg))))
scored <- function(keep) do.call(rbind, lapply(STAGES, function(sg) { r <- stage_rows(sg); r <- r[keep(r), ]
  data.frame(dose_mg = r$dose_mg, time_h = r$time_h, occ = r$predicted, curve = lab_of(sg)) }))
scoredB <- scored(function(r) grepl("Rabiner", r$study) & r$regimen == "single")
pb <- ggplot() + geom_line(data = curveB, aes(dose_mg, occ, colour = curve, linetype = curve)) +
  geom_point(data = scoredB, aes(dose_mg, occ, colour = curve), shape = 1, size = 2, show.legend = FALSE) +
  geom_point(data = rab, aes(dose_mg, occupancy_pct), colour = "black", size = 1.8) + scale_stage +
  scale_x_log10() + coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Single dose (mg, log scale)", y = sprintf("MOR occupancy at %g h (%%)", rab$time_h[1]), title = "B  Rabiner dose-response")
note("fig2_pk_pet", "B", sprintf("Black: Rabiner 2011 single-dose occupancy at %g h. Lines: model occupancy for the same regimen and scan time; hollow markers: the row-level predictions scored in script 03. Stage A (bottom-up) fits no PET row and is the only genuine held-out prediction; stage A2 is calibrated on these Rabiner rows; stage J is calibrated on all PET rows, so its agreement is a calibration, not a validation.", rab$time_h[1]),
     c("pet_test.csv", "binding_params_A.rds", "binding_params_A2.rds", "binding_params_J.rds", "pk_fit.rds"))

lee <- obsP[grepl("Lee", obsP$study) & obsP$reference == PET_NORM_REF, ]
if (!nrow(lee)) stop("expected normalised Lee 1988 rows in pet_test.csv", call. = FALSE)
if (length(unique(lee$dose_mg)) != 1) stop("expected one Lee 1988 dose", call. = FALSE)
wee <- obsP[grepl("Weerts", obsP$study), ]
if (nrow(wee) != 1) stop("expected one Weerts 2008 row in pet_test.csv", call. = FALSE)
gridC <- data.frame(dose_mg = lee$dose_mg[1], regimen = "single", time_h = seq(1, 168, 1), reference = PET_NORM_REF, n_doses = 1)
curveC <- do.call(rbind, lapply(STAGES, function(sg) data.frame(time_h = gridC$time_h, occ = predict_pet(bps[[sg]], gridC, pk$par), curve = lab_of(sg))))
scoredW <- scored(function(r) grepl("Weerts", r$study))
scoredL <- scored(function(r) grepl("Lee", r$study) & r$reference == PET_NORM_REF)
ymaxC <- max(115, max(curveC$occ) + 8)
pc <- ggplot() + geom_line(data = curveC, aes(time_h, occ, colour = curve, linetype = curve)) +
  geom_point(data = scoredL, aes(time_h, occ, colour = curve), shape = 1, size = 2, show.legend = FALSE) +
  geom_point(data = lee, aes(time_h, occupancy_pct), colour = "black", size = 1.8) +
  geom_errorbar(data = wee, aes(x = time_h, ymin = occupancy_pct - sd, ymax = occupancy_pct + sd), colour = "black", width = 4) +
  geom_point(data = wee, aes(time_h, occupancy_pct), colour = "black", shape = 17, size = 2.2) +
  geom_point(data = scoredW, aes(time_h, occ, colour = curve), shape = 2, size = 2.2, show.legend = FALSE) + scale_stage +
  # Weerts label in the empty top-right corner with a leader line to its marker, so it never sits on a curve or a stage
  # marker (R4 fix round 1, review M11; checked rendered)
  annotate("segment", x = wee$time_h + 4, y = wee$occupancy_pct + 1, xend = wee$time_h + 44, yend = ymaxC - 4, colour = "grey40", linewidth = 0.3) +
  annotate("text", x = wee$time_h + 46, y = ymaxC - 4, label = "Weerts 2008", hjust = 0, size = 2.3) +
  coord_cartesian(ylim = c(0, ymaxC)) +
  labs(x = sprintf("Time after %g mg (h)", lee$dose_mg[1]), y = "Occupancy (% of 1 h scan; Weerts absolute)", title = "C  Lee time course")
note("fig2_pk_pet", "C", sprintf("Black circles: Lee 1988, %g mg single dose, percent blockade normalised to that study's own 1 h scan (curves on the same estimand, 100 x RO(t) / RO(1 h)). Black triangle: Weerts 2008, absolute occupancy (mean +/- SD) after %d once-daily %g mg doses, at an ASSUMED %g h after the last dose (scan interval not reported; pet_weerts_scan_time_sensitivity.csv); hollow triangles: each stage's prediction for that row; hollow circles: the scored Lee predictions. A normalised curve can exceed 100%%: under a slow association rate, occupancy is still rising at the 1 h reference scan, so later occupancy exceeds the reference.",
                                 lee$dose_mg[1], as.integer(wee$n_doses), wee$dose_mg, wee$time_h),
     c("pet_test.csv", "pet_weerts_scan_time_sensitivity.csv"))
p2 <- pa + (pb + pc + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.direction = "vertical",
                                                             legend.text = element_text(size = 6.5))) +
  plot_layout(widths = c(1, 2))
save_fig(p2, "fig2_pk_pet", 180, 105)

# ============ Figure 3: training fits (posterior-predicted arm outcomes) ============
arms <- read_dataset("data/training/efficacy_training.csv")
# Only arms in the fitted likelihood (finite se; see script 04) have a trial effect in the posterior.
arms <- arms[is.finite(arms$se), , drop = FALSE]
need(file.path(res, sprintf("posterior_%s.rds", c("H0", "H1", "H2", "H3"))))
cache <- readRDS(file.path(res, sprintf("exposure_cache_%s.rds", stage)))
pred_rows <- do.call(rbind, lapply(c("H0", "H1", "H2", "H3"), function(h) {
  f <- readRDS(file.path(res, sprintf("posterior_%s.rds", h)))
  idx <- round(seq(1, nrow(f$draws), length.out = 300))
  do.call(rbind, lapply(seq_len(nrow(arms)), function(i) {
    a <- arms[i, ]
    v <- vapply(idx, function(k) { par <- f$draws[k, ]
      sim <- simulate_arm(h, theta_from(par, h, f$fixed), parse_schedule(a$schedule), a$day, cache, pmax_trial(par, a$trial))
      arm_outcome(sim, a$day, a$outcome, a$baseline_nrs) }, 0)
    data.frame(hypothesis = h, trial = a$trial, arm = a$arm, day = a$day, med = stats::median(v),
               lo = unname(stats::quantile(v, 0.05)), hi = unname(stats::quantile(v, 0.95)))
  }))
}))
# Readable arm labels from the data (R/reporting.R arm_label), the SAME function and columns for the observed arms and the
# predictions, so both layers share their x categories (review I4), e.g. "Younger 2013\nLDN 4.5 mg (day 84)".
arms$arm_id <- arm_label(arms$trial, arms$arm, arms$day)
pred_rows$arm_id <- arm_label(pred_rows$trial, pred_rows$arm, pred_rows$day)
if (!identical(sort(unique(pred_rows$arm_id)), sort(unique(arms$arm_id))))
  stop("Fig 3: observed and predicted arms do not share their categories", call. = FALSE)
pd <- position_dodge(width = 0.6)
p3 <- ggplot() +
  geom_pointrange(data = pred_rows, aes(arm_id, med, ymin = lo, ymax = hi, colour = hypothesis), position = pd, size = 0.35) +
  geom_point(data = arms, aes(arm_id, mean), colour = "black", shape = 15, size = 1.8, position = position_nudge(x = -0.42)) +
  geom_linerange(data = arms, aes(arm_id, mean, ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                 colour = "black", position = position_nudge(x = -0.42)) +
  scale_h + labs(x = NULL, y = "Pain reduction (%)", title = "Training fit: observed mean (black square, 95% CI) vs posterior prediction (median, 90% interval)") +
  theme(plot.title = element_text(size = 8))
save_fig(p3, "fig3_training_fits", 180, 100)
note("fig3_training_fits", "", sprintf("Posterior predictions use the adopted stage-%s exposure cache; 300 posterior draws per hypothesis.", stage),
     c("posterior_H0-H3.rds", sprintf("exposure_cache_%s.rds", stage), "data/training/efficacy_training.csv"))

# ============ Figure 4: held-out trial predictions ============
tp <- rd("test_predictions.csv")
need_cols(tp, c("hypothesis", "trial", "endpoint", "day", "role", "observed", "obs_se", "pred_median", "pred_lo90", "pred_hi90",
                "true_lo90", "true_hi90"), "test_predictions.csv")
dd <- tp[tp$endpoint == "diff_nrs", , drop = FALSE]
if (!nrow(dd)) stop("test_predictions.csv has no diff_nrs rows", call. = FALSE)
trial_lab <- function(t) toupper(t)
dd$trial_f <- factor(trial_lab(dd$trial), levels = unique(trial_lab(dd$trial[order(match(tolower(dd$trial), c("final", "innova")), dd$day)])))
# rows: trial, then day, then hypothesis (observed first); top to bottom
keys <- unique(dd[order(dd$trial_f, dd$day), c("trial_f", "day", "role")])
lev <- unlist(lapply(seq_len(nrow(keys)), function(i) {
  dl <- sprintf("day %g%s", keys$day[i], if (keys$role[i] == "primary") " (primary)" else "")
  paste(dl, c("observed", "H0", "H1", "H2", "H3"), sep = " | ") }))
dd$row <- sprintf("day %g%s | %s", dd$day, ifelse(dd$role == "primary", " (primary)", ""), dd$hypothesis)
ob <- unique(dd[, c("trial_f", "day", "role", "observed", "obs_se")])
ob$row <- sprintf("day %g%s | observed", ob$day, ifelse(ob$role == "primary", " (primary)", ""))
ob$lo95 <- ob$observed - 1.96 * ob$obs_se; ob$hi95 <- ob$observed + 1.96 * ob$obs_se
# numeric row positions (a discrete axis would regroup rows that are absent from the first layer)
lev <- unique(lev); ypos <- function(r) length(lev) - match(r, lev) + 1
if (anyNA(ypos(dd$row)) || anyNA(ypos(ob$row))) stop("held-out figure: unmatched row labels", call. = FALSE)
dd$row <- ypos(dd$row); ob$row <- ypos(ob$row)
p4 <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_linerange(data = dd, aes(y = row, xmin = pred_lo90, xmax = pred_hi90, colour = hypothesis, linewidth = "posterior predictive 90% (observed difference)"),
                 alpha = 0.35, orientation = "y") +
  geom_linerange(data = dd, aes(y = row, xmin = true_lo90, xmax = true_hi90, colour = hypothesis, linewidth = "90% interval of the true difference"),
                 orientation = "y") +
  geom_point(data = dd, aes(pred_median, row, colour = hypothesis), size = 1.4) +
  geom_linerange(data = ob, aes(y = row, xmin = lo95, xmax = hi95), colour = "black", linewidth = 0.5, orientation = "y") +
  geom_point(data = ob, aes(observed, row), colour = "black", shape = 15, size = 1.9) +
  facet_grid(trial_f ~ ., scales = "free_y", space = "free_y") +
  scale_y_continuous(breaks = length(lev):1, labels = lev, expand = expansion(add = 0.6), minor_breaks = NULL) +
  scale_h + scale_linewidth_manual(values = c("posterior predictive 90% (observed difference)" = 2.6, "90% interval of the true difference" = 0.6), name = NULL) +
  guides(colour = guide_legend(order = 1), linewidth = guide_legend(order = 2, nrow = 2)) +
  labs(x = "LDN minus placebo, change in pain NRS (points)", y = NULL,
       title = "Held-out trials: observed difference (black, 95% CI) against each hypothesis's predictions") +
  theme(plot.title = element_text(size = 8), axis.text.y = element_text(size = 6.5), legend.box = "vertical")
save_fig(p4, "fig4_heldout", 180, 150)
note("fig4_heldout", "", paste("Rows ordered by trial, then follow-up day, then hypothesis. Black: observed LDN minus placebo difference in change from baseline with its 95% CI",
  "(observed +/- 1.96 x the SE used for scoring; for INNOVA the unadjusted difference of the arm changes, SE borrowed from the adjusted 95% CI). Wide pale bar: posterior",
  "predictive 90% interval for the OBSERVED difference (model draws plus observation error), against which inside90 is judged; thin bar: 90% interval of the TRUE",
  "difference (model draws only); point: predictive median. Primary endpoints: FINAL day 84, INNOVA day 90; other INNOVA days are secondary and correlated."),
  c("test_predictions.csv"))

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
note("fig5_cytokines", "", "Cytokines without a curated half-life have no prediction and are not plotted (listed in Table 3).", "cytokine_test.csv")

# ============ Figure 6: virtual trials, simulated difference ============
need(file.path(res, "virtual_trials.csv"))
hdr <- names(utils::read.csv(file.path(res, "virtual_trials.csv"), nrows = 1))
need_cols(stats::setNames(as.list(hdr), hdr), c("diff", "hypothesis", "trial", "omega_mult"), "virtual_trials.csv")
vt <- utils::read.csv(file.path(res, "virtual_trials.csv"), stringsAsFactors = FALSE,
                      colClasses = ifelse(hdr %in% c("diff", "omega_mult"), "numeric", ifelse(hdr %in% c("hypothesis", "trial"), "character", "NULL")))
vs <- rd("virtual_summary.csv"); need_cols(vs, c("hypothesis", "trial", "omega_mult", "obs_diff", "obs_diff_adjusted_level_sensitivity", "ppp_two_sided"), "virtual_summary.csv")
vt <- vt[vt$omega_mult == 1, ]; vt$trial <- trial_lab(vt$trial)
v1 <- vs[vs$omega_mult == 1, ]; v1$trial <- trial_lab(v1$trial)
ob6 <- unique(v1[, c("trial", "obs_diff")])
ob6s <- unique(v1[is.finite(v1$obs_diff_adjusted_level_sensitivity), c("trial", "obs_diff_adjusted_level_sensitivity")])
v1$lab <- sprintf("p = %.2f", v1$ppp_two_sided)
n_sim <- table(paste(vt$trial, vt$hypothesis))          # simulated differences per violin, computed (review I9)
n_sim_txt <- if (length(unique(n_sim)) == 1) format(unique(as.integer(n_sim)), big.mark = ",") else
  sprintf("%s-%s", format(min(n_sim), big.mark = ","), format(max(n_sim), big.mark = ","))
p6 <- ggplot(vt, aes(hypothesis, diff, fill = hypothesis)) + geom_violin(colour = NA, alpha = 0.8, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_hline(data = ob6, aes(yintercept = obs_diff, linetype = "Observed (unadjusted)"), colour = "black", linewidth = 0.6, key_glyph = "path") +
  geom_hline(data = ob6s, aes(yintercept = obs_diff_adjusted_level_sensitivity, linetype = "INNOVA adjusted (sensitivity)"), colour = "black", linewidth = 0.4, key_glyph = "path") +
  geom_hline(data = data.frame(y = 0), aes(yintercept = y, linetype = "Zero reference"), colour = "grey50", linewidth = 0.3, key_glyph = "path") +
  geom_text(data = v1, aes(hypothesis, Inf, label = lab), vjust = 1.5, size = 2.3, inherit.aes = FALSE) +
  facet_wrap(~trial) + scale_fill_manual(values = pal_h, guide = "none") +
  scale_linetype_manual(values = c("Observed (unadjusted)" = "solid", "INNOVA adjusted (sensitivity)" = "dashed", "Zero reference" = "dotted"),
                        name = NULL) +
  labs(x = NULL, y = "Simulated difference in NRS change (LDN - placebo)",
       title = "Virtual replicate trials at the actual designs",
       subtitle = sprintf("Violin and box: %s simulated LDN - placebo differences per hypothesis and trial; p: two-sided posterior predictive p-value", n_sim_txt)) +
  theme(plot.title = element_text(size = 8), plot.subtitle = element_text(size = 6.5),
        legend.position = "bottom", legend.text = element_text(size = 6.5), legend.key.width = unit(8, "mm"))
save_fig(p6, "fig6_virtual_trials", 180, 95)
note("fig6_virtual_trials", "", "Violins: simulated LDN - placebo differences at each trial's actual design (omega as calibrated). Solid line: observed (unadjusted) difference; dashed: INNOVA's LMM-adjusted level difference (sensitivity estimand). p: two-sided posterior predictive p-value of the observed difference (2 x the smaller tail).",
     c("virtual_trials.csv", "virtual_summary.csv"))

# ============ Figure 7: occupancy ============
ssJ <- rd(sprintf("steady_state_occupancy_%s.csv", stage)); need_cols(ssJ, c("dose_mg", "mean_RO"), sprintf("steady_state_occupancy_%s.csv", stage))
unc <- rd(sprintf("occupancy_uncertainty_%s.csv", stage))
need_cols(unc, c("dose_mg", "uncertainty_source", "mean_RO_point", "mean_RO_lo95", "mean_RO_hi95"), sprintf("occupancy_uncertainty_%s.csv", stage))
osn <- rd(sprintf("occupancy_sensitivity_%s.csv", stage))
need_cols(osn, c("key", "analysis", "fp", "mean_ss_occupancy_4.5mg", "wsse"), sprintf("occupancy_sensitivity_%s.csv", stage))
byst <- rd("steady_state_occupancy_by_stage.csv"); need_cols(byst, c("stage", "dose_mg", "mean_RO"), "steady_state_occupancy_by_stage.csv")
curve7 <- ssJ[ssJ$dose_mg >= 0.75 & ssJ$dose_mg <= 6, ]
band <- unc[unc$uncertainty_source %in% c("binding_only", "binding_and_pk"), ]
if (!nrow(band)) stop(sprintf("occupancy_uncertainty_%s.csv has no Laplace band rows", stage), call. = FALSE)
band$kind <- ifelse(band$uncertainty_source == "binding_only", "error bar: 95% band, binding (Laplace)", "error bar: 95% band, binding + PK (Laplace)")
band$x <- band$dose_mg + ifelse(band$uncertainty_source == "binding_only", -0.06, 0.06)
rab_ref <- rab[rab$dose_mg %in% c(2, 5), ]
p7a <- ggplot() +
  geom_line(data = curve7, aes(dose_mg, 100 * mean_RO), colour = col_stage[[lab_of(stage)]]) +
  geom_point(data = curve7, aes(dose_mg, 100 * mean_RO), colour = col_stage[[lab_of(stage)]], size = 1.2) +
  geom_errorbar(data = band, aes(x = x, ymin = 100 * mean_RO_lo95, ymax = 100 * mean_RO_hi95, linetype = kind), width = 0.12,
                colour = col_stage[[lab_of(stage)]]) +
  geom_point(data = rab_ref, aes(dose_mg, occupancy_pct, shape = "Rabiner 2011, single dose at 4 h
(reference, not a steady-state mean)"), colour = "black", size = 2) +
  scale_shape_manual(values = c(1), name = NULL) + scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  guides(linetype = guide_legend(nrow = 2), shape = guide_legend(nrow = 1)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Dose, once daily (mg)", y = "Mean steady-state MOR occupancy (%)", title = sprintf("A  Stage %s occupancy vs dose", stage)) +
  theme(legend.box = "vertical", legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.5, "mm"), legend.key.height = unit(3, "mm"))
at45 <- function(d) abs(d - 4.5) < 1e-9
pt45 <- unc$mean_RO_point[at45(unc$dose_mg)][1]
# Rows are selected by their stable key (PK_SENS_KEYS / OCC_SENS_KEYS in R/pk.R), never by label text (review C1).
okey <- function(k) osn[osn$key %in% k, , drop = FALSE]
fp_rows <- okey(PK_SENS_KEYS[c("fp090", "fp099")]); term <- okey(PK_SENS_KEYS[["terminal"]])
keq <- unc[unc$uncertainty_source == "keq_sensitivity" & at45(unc$dose_mg), ]
eqk <- okey(OCC_SENS_KEYS[["koff_eq"]]); prim <- okey(OCC_SENS_KEYS[["primary"]])
d50_rows <- okey(PK_SENS_KEYS[c("D50_10", "D50_100")]); pk_dd <- pk_dose_model_of(pk$par) != "linear"
d50_bd <- if (pk_dd) d50_boundary(res, pk_dose_model_of(pk$par)) else numeric(0)
for (x in list(list(term, "terminal-phase row"), list(keq, "keq_sensitivity row"), list(eqk, "koff_BN-equal row"), list(prim, "primary row")))
  if (nrow(x[[1]]) != 1) stop("occupancy sensitivity inputs: expected one ", x[[2]], call. = FALSE)
if (nrow(fp_rows) != 2) stop("occupancy_sensitivity: expected the fp 0.90 and 0.99 rows (keys ", paste(PK_SENS_KEYS[c("fp090", "fp099")], collapse = ", "), ")", call. = FALSE)
if (pk_dd && nrow(d50_rows) != 2) stop("occupancy_sensitivity: a dose-dependent PK is adopted but the D50 10 and 100 mg rows are missing", call. = FALSE)
if (!pk_dd && nrow(d50_rows)) stop("occupancy_sensitivity: D50 rows present although the linear PK is adopted", call. = FALSE)
st45 <- function(sg) { v <- byst$mean_RO[byst$stage == sg & at45(byst$dose_mg)]; if (length(v) != 1) stop("no 4.5 mg row for stage ", sg, call. = FALSE); v }
b45 <- band[at45(band$dose_mg), ]
rng <- rbind(
  data.frame(item = "Laplace 95% band, binding only", lo = b45$mean_RO_lo95[b45$uncertainty_source == "binding_only"], hi = b45$mean_RO_hi95[b45$uncertainty_source == "binding_only"], kind = "probability band"),
  data.frame(item = "Laplace 95% band, binding + PK", lo = b45$mean_RO_lo95[b45$uncertainty_source == "binding_and_pk"], hi = b45$mean_RO_hi95[b45$uncertainty_source == "binding_and_pk"], kind = "probability band"),
  data.frame(item = "keq 0.5 to 8 per h (refitted)", lo = keq$mean_RO_lo95, hi = keq$mean_RO_hi95, kind = "structural sensitivity"),
  data.frame(item = sprintf("fp at 100 mg %s (refitted)", paste(sprintf("%.2f", fp_rows$fp), collapse = " / ")), lo = min(fp_rows$mean_ss_occupancy_4.5mg), hi = max(fp_rows$mean_ss_occupancy_4.5mg), kind = "structural sensitivity"),
  if (pk_dd) data.frame(item = sprintf("D50 %s mg (refitted; primary %g)%s", paste(sprintf("%g%s", d50_rows$fp_D50_mg, ifelse(d50_rows$fp_D50_mg %in% d50_bd, "*", "")), collapse = " / "),
                                       pk$par$fp_D50, if (any(d50_rows$fp_D50_mg %in% d50_bd)) "; * at boundary" else ""),
                        lo = min(d50_rows$mean_ss_occupancy_4.5mg), hi = max(d50_rows$mean_ss_occupancy_4.5mg), kind = "structural sensitivity"),
  data.frame(item = "terminal-phase PK (refitted)", lo = term$mean_ss_occupancy_4.5mg, hi = term$mean_ss_occupancy_4.5mg, kind = "structural sensitivity"),
  data.frame(item = sprintf("koff_BN base = koff_NTX (PET wSSE %.1f vs %.1f)", eqk$wsse, prim$wsse), lo = eqk$mean_ss_occupancy_4.5mg, hi = eqk$mean_ss_occupancy_4.5mg, kind = "structural sensitivity"),
  data.frame(item = "stage A2 (Kd on Rabiner only)", lo = st45("A2"), hi = st45("A2"), kind = "outer bound (other binding stage)"),
  data.frame(item = "stage A (bottom-up, in vitro)", lo = st45("A"), hi = st45("A"), kind = "outer bound (other binding stage)"))
rng$item <- factor(rng$item, levels = rev(rng$item))
p7b <- ggplot(rng) + geom_vline(xintercept = 100 * pt45, colour = col_stage[[lab_of(stage)]], linetype = "dashed") +
  geom_linerange(aes(y = item, xmin = 100 * lo, xmax = 100 * hi, colour = kind), linewidth = 1.6, orientation = "y") +
  geom_point(aes(y = item, x = 100 * lo, colour = kind), size = 1.6) + geom_point(aes(y = item, x = 100 * hi, colour = kind), size = 1.6) +
  scale_colour_manual(values = c("probability band" = "#0072B2", "structural sensitivity" = "#E69F00", "outer bound (other binding stage)" = "#999999"),
                      breaks = c("probability band", "structural sensitivity", "outer bound (other binding stage)"), name = NULL) +
  guides(colour = guide_legend(nrow = 3)) + coord_cartesian(xlim = c(0, 100)) +
  labs(x = "Mean steady-state occupancy at 4.5 mg (%)", y = NULL, title = "B  4.5 mg: bands and sensitivities") +
  theme(axis.text.y = element_text(size = 6.5), legend.text = element_text(size = 6), plot.title = element_text(size = 9),
        legend.spacing.y = unit(0.5, "mm"), legend.key.height = unit(3, "mm"))
save_fig(p7a + p7b + plot_layout(widths = c(1, 1.1)), "fig7_occupancy", 180, 95)
note("fig7_occupancy", "A", sprintf("Line: stage-%s mean steady-state occupancy over one dosing interval (certified steady state). Error bars at 1.5, 4.5 and 6 mg: 95%% Laplace bands for the binding scales alone and jointly with independent PK draws (the joint band ignores the PK-binding correlation, so it is approximate, not a strict superset of the binding-only band). Open circles: Rabiner 2011 occupancy at 2 and 5 mg, single doses scanned at 4 h; they are reference points, not steady-state means.", stage),
     c(sprintf("steady_state_occupancy_%s.csv", stage), sprintf("occupancy_uncertainty_%s.csv", stage), "pet_test.csv"))
shift_pp <- function(v) 100 * (v - pt45)
note("fig7_occupancy", "B", sprintf("Dashed line: the adopted stage's point value. Probability bands (blue) and structural sensitivities (orange) are reported side by side; the structural rows are refits, not probability statements. Stages A2 and A bound the binding assumption from outside. Shifts from the point value: terminal-phase PK %+.1f pp; koff_BN base = koff_NTX %+.1f pp (and a worse PET fit, weighted SSE %.1f against %.1f); fp at 100 mg 0.90-0.99 %+.1f to %+.1f pp; keq 0.5-8 per h %+.1f to %+.1f pp; %s.",
  shift_pp(term$mean_ss_occupancy_4.5mg), shift_pp(eqk$mean_ss_occupancy_4.5mg), eqk$wsse, prim$wsse,
  min(shift_pp(fp_rows$mean_ss_occupancy_4.5mg)), max(shift_pp(fp_rows$mean_ss_occupancy_4.5mg)), shift_pp(keq$mean_RO_lo95), shift_pp(keq$mean_RO_hi95),
  if (pk_dd) sprintf("PK D50 %s mg (adopted %s, primary D50 %g mg) %+.1f to %+.1f pp", paste(sprintf("%g", d50_rows$fp_D50_mg), collapse = "/"),
                     PK_DOSE_MODELS[[pk_dose_model_of(pk$par)]]$label, pk$par$fp_D50,
                     min(shift_pp(d50_rows$mean_ss_occupancy_4.5mg)), max(shift_pp(d50_rows$mean_ss_occupancy_4.5mg)))
  else "no D50 item: the linear PK is adopted, so the dose-dependent D50 fits are PK outputs only (pk_dose_model_comparison.csv)"),
     c(sprintf("occupancy_uncertainty_%s.csv", stage), sprintf("occupancy_sensitivity_%s.csv", stage), "steady_state_occupancy_by_stage.csv", "binding_J_keq_sensitivity.csv"))

# ============ Figure 8: assurance and the type-I-error check ============
pw <- rd("power_curve.csv")
need_cols(pw, c("hypothesis", "role", "trial", "n_per_arm", "rej_rate_with_model_unc", "rej_rate_sampling_only"), "power_curve.csv")
pw <- pw[is.finite(pw$rej_rate_with_model_unc), ]                       # explanatory rows (e.g. an ineligible H2) carry no rate
is_type1 <- function(d) d$hypothesis == "H0" | grepl("type I", d$role, ignore.case = TRUE)
pw$trial <- trial_lab(pw$trial)
pa8 <- pw[!is_type1(pw), ]; p08 <- pw[is_type1(pw), ]
if (!nrow(p08)) stop("power_curve.csv has no H0 type-I-error rows", call. = FALSE)
role_lab <- function(d) sprintf("%s: %s", d$hypothesis, sub(" \\(net effect.*$", "", d$role))
p8a <- if (nrow(pa8)) {
  pa8$curve <- role_lab(pa8)
  ggplot(pa8, aes(n_per_arm, rej_rate_with_model_unc, colour = hypothesis, group = curve)) + geom_line() + geom_point(size = 1.3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.4) + facet_wrap(~trial) +
    scale_colour_manual(values = pal_h, labels = function(h) unique(pa8$curve[match(h, pa8$hypothesis)]), name = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Participants per arm", y = "Assurance (posterior-averaged power)", title = sprintf("A  %s: assurance of a one-sided 5%% test (dashed: 0.05)", paste(sort(unique(pa8$hypothesis)), collapse = " and ")))
} else ggplot() + annotate("text", 0, 0, label = "No hypothesis exceeds the effect threshold (effect_eligibility.csv)") + theme_void()
# R4: add linetype key for fig8B (solid = with model uncertainty; dotted = sampling variability only)
p8b <- ggplot(p08, aes(n_per_arm, rej_rate_with_model_unc,
                       linetype = "With model-structure uncertainty")) +
  geom_line(colour = pal_h[["H0"]]) + geom_point(colour = pal_h[["H0"]], size = 1.3) +
  geom_line(aes(y = rej_rate_sampling_only, linetype = "Sampling variability only"), colour = pal_h[["H0"]]) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.4) + facet_wrap(~trial) +
  scale_linetype_manual(values = c("With model-structure uncertainty" = "solid", "Sampling variability only" = "dotted"),
                        name = NULL) +
  coord_cartesian(ylim = c(0, max(0.1, max(p08$rej_rate_with_model_unc) + 0.01))) +
  labs(x = "Participants per arm", y = "Rejection rate under H0", title = "B  H0: type I error check (nominal 0.05, dashed)") +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.5), legend.key.width = unit(8, "mm"))
save_fig(p8a / p8b + plot_layout(heights = c(1.3, 1)), "fig8_power", 180, 140)
d_so <- max(abs(pw$rej_rate_with_model_unc - pw$rej_rate_sampling_only), na.rm = TRUE)
note("fig8_power", "", sprintf("Rejection rates of a one-sided 5%% test whose critical value comes from H0 replicates (independent seed). A: hypotheses whose net effect at day 84, evaluated at the posterior-median parameters, exceeds the threshold (effect_eligibility.csv); their rate averages power over the posterior, i.e. assurance. B: H0 as the type-I-error check; solid with model-structure uncertainty, dotted sampling variability only (the two variants differ by at most %.3f anywhere). Hypotheses below the threshold have no curve.", d_so),
     c("power_curve.csv", "effect_eligibility.csv"))

# ============ Figure S1: Sobol indices ============
sb <- rd("sobol.csv")
need_cols(sb, c("hypothesis", "day", "parameter", "S1", "S1_lo", "S1_hi", "ST", "ST_lo", "ST_hi", "eligible"), "sobol.csv")
sb <- sb[as.logical(sb$eligible) %in% TRUE & is.finite(sb$S1), ]
if (!nrow(sb)) stop("sobol.csv has no indices for an eligible hypothesis", call. = FALSE)
SOBOL_SEED_FILE <- "sobol_seed_agreement.csv"
seedf <- file.path(res, SOBOL_SEED_FILE); need(seedf)
sa <- utils::read.csv(seedf, stringsAsFactors = FALSE, check.names = FALSE)
need_cols(sa, c("hypothesis", "day", "parameter", "index", "disagree"), SOBOL_SEED_FILE)
sl <- rbind(data.frame(sb[, c("hypothesis", "day", "parameter")], index = "S1", est = sb$S1, lo = sb$S1_lo, hi = sb$S1_hi),
            data.frame(sb[, c("hypothesis", "day", "parameter")], index = "ST", est = sb$ST, lo = sb$ST_lo, hi = sb$ST_hi))
sl <- merge(sl, sa[, c("hypothesis", "day", "parameter", "index", "disagree")], by = c("hypothesis", "day", "parameter", "index"), all.x = TRUE)
if (anyNA(sl$disagree)) stop(SOBOL_SEED_FILE, " lacks rows for some plotted indices", call. = FALSE)
sl$seeds_agree <- !as.logical(sl$disagree)
sl$index <- factor(sl$index, levels = c("S1", "ST"), labels = c("first order S1 (Sobol 2007)", "total ST (Jansen)"))
sl$panel <- sprintf("%s, day %g", sl$hypothesis, sl$day)
# R4: map raw parameter codes to readable labels; unknown codes shown as-is with underscores replaced
par_readable <- c(ka = "ka (abs. rate, /h)", fp = "fp (presystemic conversion at 100 mg)", fp0 = "fp0 (presystemic conversion, dose -> 0)",
  floss0 = "floss0 (parent-only presystemic loss, dose -> 0)",
  # sampler-scale names as they appear in sobol.csv
  logit_Emax = "Emax (logit scale)", log_kpbo = "kpbo, placebo onset rate (log scale)", log_kout = "kout (log scale)",
  log_kA = "kA (log scale)", log_gamma = "gamma (log scale)", log10_Kd_OGFr = "Kd OGFr (log10 scale)",
  CL = "CL (apparent, L/h)", V1 = "V1 (apparent, L)", Q = "Q (intercompt., L/h)", V2 = "V2 (periph., L)",
  CLm = "CLm (metabolite, L/h)", Vm = "Vm (metabolite, L)", fm = "fm (fraction to metabolite)",
  Kd_NTX = "Kd NTX (MOR, nM)", Kd_BN = "Kd 6β-NTX (MOR, nM)", koff_NTX = "koff NTX (/h)", koff_BN = "koff 6β-NTX (/h)",
  keq = "keq (equilibration, /h)", kd_scale = "Kd scale (binding calib.)", koff_scale = "koff scale (binding calib.)",
  E0 = "E0 (baseline pain)", Emax = "Emax (max drug effect)", EC50 = "EC50 (occupancy at half-max)",
  Hill = "Hill coeff.", sigma_e = "σ_e (obs. noise)", omega_p = "ω_p (placebo variability)",
  tau_p = "τ_p (placebo autocorr.)", alpha = "α (temporal decay)", beta = "β (dose-response slope)",
  IC50_TLR4 = "IC50 TLR4 (NTX)", Emax_cyt = "Emax cytokine", kprod = "k_prod (cytokine prod.)",
  kdeg = "k_deg (cytokine deg.)", kOGF = "k_OGF (OGFr pathway)")
# Apply labels: known codes use the lookup, unknown codes get underscores → spaces for readability
sl$par_label <- ifelse(sl$parameter %in% names(par_readable), par_readable[sl$parameter],
                       gsub("_", " ", sl$parameter))
pS1 <- ggplot(sl, aes(est, par_label, colour = index)) + geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.55), size = 0.3) +
  geom_point(data = sl[!sl$seeds_agree, ], aes(x = hi, group = index), shape = 4, colour = "black", size = 2.4, position = position_dodge(width = 0.55)) +
  facet_wrap(~panel, ncol = 2, scales = "free_y") + scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
  labs(x = "Sobol index of the net effect (95% bootstrap CI)", y = NULL,
       title = "Sobol indices, hypotheses with an effect; x marks an index on which the two seeds disagree") +
  theme(plot.title = element_text(size = 8))
save_fig(pS1, "figS1_sobol", 180, 130)
note("figS1_sobol", "", paste("First-order indices by the Sobol 2007 estimator and total indices by Jansen's estimator, both from the same model runs (Saltelli et al. 2010).",
  "Jansen's first-order estimate (S1_jansen) is high-variance on these heavy-tailed outputs (its variance term comes from one sample) and is not shown. Parameters vary uniformly over their marginal 95% posterior",
  "intervals; only hypotheses whose net effect exceeds the threshold are analysed. Each index was computed with two seeds; x marks an index whose two estimates differ by more than the larger bootstrap CI half-width",
  sprintf("(%d of %d plotted indices; %s).", sum(!sl$seeds_agree), nrow(sl), SOBOL_SEED_FILE)),
  c("sobol.csv", SOBOL_SEED_FILE))

# ============ Figure S2: power under the SD-sensitivity calibration (consistency check) ============
pws <- rd("power_curve_sd_sensitivity.csv"); need_cols(pws, names(pw)[names(pw) %in% c("hypothesis", "role", "trial", "n_per_arm", "rej_rate_with_model_unc")], "power_curve_sd_sensitivity.csv")
pws <- pws[is.finite(pws$rej_rate_with_model_unc), ]; pws$trial <- trial_lab(pws$trial)
both <- rbind(data.frame(pw[, c("hypothesis", "trial", "n_per_arm", "rej_rate_with_model_unc")], calibration = "primary"),
              data.frame(pws[, c("hypothesis", "trial", "n_per_arm", "rej_rate_with_model_unc")], calibration = "SD sensitivity (consistency check, not a candidate calibration)"))
pS2 <- ggplot(both, aes(n_per_arm, rej_rate_with_model_unc, colour = hypothesis, linetype = calibration, group = interaction(hypothesis, calibration))) +
  geom_line() + geom_point(size = 1.1) + geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.3) + facet_wrap(~trial) +
  scale_h + scale_linetype_manual(values = c("solid", "22"), name = NULL) + guides(linetype = guide_legend(nrow = 2)) +
  coord_cartesian(ylim = c(0, 1)) + labs(x = "Participants per arm", y = "Rejection rate",
                                         title = "Rejection rates (H1, H3: assurance; H0: type I error) under the primary and the SD-sensitivity calibration") +
  theme(plot.title = element_text(size = 8), legend.box = "vertical")
save_fig(pS2, "figS2_power_sd_sensitivity", 180, 100)
note("figS2_power_sd_sensitivity", "", "SD sensitivity: variability recalibrated to the Younger 2013 SDs quoted by Bested et al. 2023 (FIQ pain reduction, sample-size section; 14.6 placebo, 12.5 LDN). It is a consistency check, not a candidate calibration; its responder rates are compared with the primary calibration's in tableS_responder_calibrations.csv.",
     c("power_curve.csv", "power_curve_sd_sensitivity.csv", "omegas_sd_sensitivity.csv"))

fn <- do.call(rbind, notes); fn$file <- paste0(fn$figure, c(".png / .tiff"))
write.csv(fn[, c("figure", "file", "panel", "note", "source_files")], file.path(res, "figure_notes.csv"), row.names = FALSE)
cat("Figures written to", figdir, "; caption notes in", file.path(res, "figure_notes.csv"), "\n")
