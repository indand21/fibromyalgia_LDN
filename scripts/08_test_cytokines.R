options(v5.root = normalizePath(".", winslash = "/")); source("R/utils.R"); source_all()
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
stage <- readLines("output/results/binding_stage.txt")[1]
cache <- readRDS(sprintf("output/results/exposure_cache_%s.rds", stage))
if (is.null(cache[["4.5"]])) stop("exposure cache has no entry for the 4.5 mg dose")
obs <- read_test_data("data/test/cytokines_parkitny.csv", "B")
cp <- read_dataset("data/params/cytokine_params.csv")
out <- do.call(rbind, lapply(c("H0", "H1", "H2", "H3"), function(h) {
  f <- readRDS(sprintf("output/results/posterior_%s.rds", h))
  idx <- round(seq(1, nrow(f$draws), length.out = 500))
  # each observation row is predicted at its OWN day (review M11; all six Parkitny rows are day 56 today)
  pr <- do.call(rbind, lapply(idx, function(i) cytokine_predict_rows(h, theta_from(f$draws[i, ], h, f$fixed), cache[["4.5"]], cp, obs)))
  agg <- stats::aggregate(pct_change_pred ~ cytokine + day, pr, function(x) c(med = stats::median(x), lo = stats::quantile(x, .05), hi = stats::quantile(x, .95)))
  agg <- data.frame(cytokine = agg$cytokine, day = agg$day, agg$pct_change_pred); names(agg)[3:5] <- c("pred_median", "pred_lo90", "pred_hi90")
  agg$hypothesis <- h; agg
}))
out <- merge(out, obs[, c("cytokine", "day", "pct_change", "se", "p_value", "n")], by = c("cytokine", "day"), all = TRUE)
out$predicted <- !is.na(out$pred_median)
no_pred <- setdiff(unique(obs$cytokine), unique(out$cytokine[out$predicted]))
if (length(no_pred)) warning("no curated half-life, so no prediction for: ", paste(no_pred, collapse = ", "))
write.csv(out, "output/results/cytokine_test.csv", row.names = FALSE)

# QST (FINAL exploratory paper, Due Bruun 2025): the source trial's own between-group differences, read through
# read_test_data so the file is hashed under stage T (R3a item 5, review I11). Descriptive: the model makes no QST
# prediction; the manuscript's QST paragraph and R3b's tables cite this file instead of a hand transcription.
qst <- qst_summary(read_test_data("data/test/qst_final.csv", "B"))
write.csv(qst, "output/results/qst_test.csv", row.names = FALSE)
print(qst[, c("measure", "diff", "lo95", "hi95", "unit", "point_favours", "ci_excludes_zero")])
