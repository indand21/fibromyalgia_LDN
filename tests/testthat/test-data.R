test_that("read_dataset rejects missing or vague provenance", {
  f <- tempfile(fileext = ".csv")
  write.csv(data.frame(x = 1, source = "Literature", location = "-", verified = "yes"), f, row.names = FALSE)
  expect_error(read_dataset(f), "not a source")
  write.csv(data.frame(x = 1, source = "A 2020", location = "T1", verified = "maybe"), f, row.names = FALSE)
  expect_error(read_dataset(f), "verified")
  write.csv(data.frame(x = 1), f, row.names = FALSE)
  expect_error(read_dataset(f), "provenance")
  # R4 fix round 1: values read off a figure are labelled `digitised` (Bested 2023 Fig 4)
  write.csv(data.frame(x = 1, source = "A 2020", location = "Fig 4", verified = "digitised"), f, row.names = FALSE)
  expect_s3_class(read_dataset(f), "data.frame")
})

test_that("parse_schedule handles single dose, placebo and titration", {
  expect_equal(parse_schedule("0:4.5"), data.frame(start_day = 0, dose_mg = 4.5))
  expect_equal(parse_schedule("0:0")$dose_mg, 0)
  s <- parse_schedule("0:1.5;7:3;14:4.5;21:6")
  expect_equal(s$start_day, c(0, 7, 14, 21)); expect_equal(s$dose_mg, c(1.5, 3, 4.5, 6))
  expect_error(parse_schedule("7:3;0:1.5"), "increasing")
})

test_that("all curated CSVs pass provenance validation", {
  files <- list.files(file.path("..", "..", "data"), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("module_screen", files)]
  skip_if(length(files) == 0, "data not curated yet")
  for (f in files) expect_s3_class(read_dataset(f), "data.frame")
})

rd_data <- function(...) {
  f <- file.path("..", "..", "data", ...)
  skip_if(!file.exists(f), paste("missing", f))
  read_dataset(f)
}

test_that("pk_summary.csv has finite doses and valid stat/analyte values", {
  d <- rd_data("training", "pk_summary.csv")
  expect_true(all(is.finite(d$dose_mg)))
  # R4: C0.5 (concentration at 30 min, Bested 2023 day 1) and Css13h (mean over 1-3 h window at steady
  # state, Bested 2023 day 21) added as new stat types; see R/pk.R PK_STATS.
  expect_true(all(d$stat %in% PK_STATS))
  expect_true(all(d$analyte %in% c("NTX", "BN")))
})

test_that("pet_occupancy.csv has finite times and a valid reference convention", {
  d <- rd_data("test", "pet_occupancy.csv")
  expect_true(all(is.finite(d$time_h)))
  expect_true(all(d$reference %in% c("absolute", "normalised to 1 h scan")))
})

test_that("binding_constants.csv is keyed on analyte and koff is per hour", {
  d <- rd_data("params", "binding_constants.csv")
  expect_true("analyte" %in% names(d))
  expect_true(all(d$analyte %in% c("NTX", "BN")))
  expect_gte(sum(d$analyte == "NTX" & d$target == "MOR" & d$quantity == "Ki"), 1)
  expect_gte(sum(d$analyte == "BN" & d$target == "MOR" & d$quantity == "Ki"), 1)
  expect_true(all(d$unit[d$quantity == "koff"] == "1/h"))
})

test_that("trials_test.csv has a change_nrs row for every arm and valid arm names", {
  d <- rd_data("test", "trials_test.csv")
  expect_true(all(d$arm == "placebo" | grepl("^ldn_", d$arm)))
  # Script 06 needs, for each (trial, day) with a diff_nrs row, exactly one placebo and one ldn_* change_nrs row.
  dd <- unique(d[d$endpoint == "diff_nrs", c("trial", "day")])
  expect_gt(nrow(dd), 0)
  for (k in seq_len(nrow(dd))) {
    x <- d[d$trial == dd$trial[k] & d$day == dd$day[k] & d$endpoint == "change_nrs", ]
    lab <- paste(dd$trial[k], "day", dd$day[k])
    expect_equal(sum(x$arm == "placebo"), 1, label = paste("placebo change_nrs rows,", lab))
    expect_equal(sum(grepl("^ldn_", x$arm)), 1, label = paste("ldn change_nrs rows,", lab))
  }
  # every arm of every trial also has at least one change_nrs row
  for (tr in unique(d$trial)) {
    x <- d[d$trial == tr, ]
    for (a in unique(x$arm))
      expect_gte(sum(x$arm == a & x$endpoint == "change_nrs"), 1, label = paste("change_nrs rows for", tr, a))
  }
})

test_that("efficacy_training.csv has a scale column and finite sd_pct for Younger 2013 rows", {
  d <- rd_data("training", "efficacy_training.csv")
  expect_true("scale" %in% names(d))
  y <- d[grepl("younger.*2013", d$trial, ignore.case = TRUE), ]
  expect_gt(nrow(y), 0)
  expect_true(all(is.finite(y$sd_pct)))
})

test_that("pool_binding converts uM rows to nM before pooling", {
  bc <- data.frame(analyte = "NTX", species = "mouse", target = "TLR4", quantity = "IC50", value = 105.5, unit = "uM",
                   source = "S", location = "L", verified = "yes")
  p <- pool_binding(bc, "NTX", "TLR4", "IC50")
  expect_equal(p$value, 105500); expect_equal(p$unit, "nM")
  bc2 <- rbind(bc, transform(bc, value = 50000, unit = "nM"), transform(bc, value = 1e8, unit = "pM"))
  expect_equal(pool_binding(bc2, "NTX", "TLR4", "IC50")$value, exp(mean(log(c(105500, 50000, 1e5)))))
  expect_equal(to_nM(c(1, 1, 1000), c("uM", "nM", "pM")), c(1000, 1, 1))
})

test_that("pool_binding stops on an unrecognised concentration unit", {
  bc <- data.frame(analyte = "NTX", species = "human", target = "MOR", quantity = "Ki", value = 1, unit = "mg/L",
                   source = "S", location = "L", verified = "yes")
  expect_error(pool_binding(bc, "NTX", "MOR", "Ki"), "unrecognised")
  expect_error(to_nM(1, "furlong"), "unrecognised")
})

test_that("the curated TLR4 IC50 is consumed as nM", {
  bc <- rd_data("params", "binding_constants.csv")
  expect_equal(pool_binding(bc, "NTX", "TLR4", c("IC50", "Kd"))$value, 105500)
  expect_equal(pool_binding(bc, "NTX", "MOR", "Ki")$unit, "nM")
})

test_that("pk_summary.csv carries a valid dosing regimen on every row", {
  d <- rd_data("training", "pk_summary.csv")
  expect_true("regimen" %in% names(d))
  expect_true(all(d$regimen %in% c("single", "daily")))
  expect_true(all(d$regimen[grepl("chronic", d$study)] == "daily"))
  # R4: Bested 2023 daily rows (Css13h, day-21 steady state) also have regimen "daily"
  # though the study name has no "chronic" substring; test updated to reflect this.
  expect_true(all(d$regimen[!grepl("chronic", d$study) & !grepl("Bested", d$study)] == "single"))
  expect_true(all(d$regimen[grepl("Bested", d$study) & d$stat == "Css13h"] == "daily"))
  expect_true(all(d$regimen[grepl("Bested", d$study) & d$stat == "C0.5"] == "single"))
})

test_that("pk_summary.csv carries the Verebey 1976 24 h concentrations after the acute dose", {
  d <- rd_data("training", "pk_summary.csv"); c24 <- d[d$stat == "C24", ]
  expect_equal(nrow(c24[c24$regimen == "single", ]), 2)
  expect_equal(c24$mean[c24$analyte == "NTX" & c24$regimen == "single"], 2.1)
  expect_equal(c24$mean[c24$analyte == "BN" & c24$regimen == "single"], 17.6)
  expect_true(all(c24$unit == "ng/mL"))
})

test_that("the terminal-phase SENSITIVITY rows live in their own file, never in the primary PK file", {
  s <- rd_data("training", "pk_summary_sensitivity_terminal.csv")
  expect_equal(nrow(s), 2)
  expect_equal(sort(s$stat), c("C24", "C72")); expect_equal(s$mean[s$stat == "C24"], 2.4); expect_equal(s$mean[s$stat == "C72"], 1.7)
  expect_true(all(s$analyte == "NTX" & s$regimen == "daily" & is.na(s$sd) & s$verified == "yes"))
  expect_true(all(grepl("which dose", s$location)))
  expect_false("C72" %in% rd_data("training", "pk_summary.csv")$stat)
  expect_identical(names(s), names(rd_data("training", "pk_summary.csv")))
})

# ---- R4 fix round 1 (review C2, I5): Bested 2023 rows weighted by their reported dispersion, curation corrected ----
test_that("the Bested 2023 day-1 profile carries digitised means and SDs at 15-60 min, down-weighted to one observation per analyte", {
  d <- rd_data("training", "pk_summary.csv"); b <- d[d$study == "Bested 2023", ]
  prof <- b[b$regimen == "single", ]
  expect_equal(nrow(prof), 8)
  for (an in c("NTX", "BN")) {
    x <- prof[prof$analyte == an, ]
    expect_setequal(x$stat, c("C0.25", "C0.5", "C0.75", "C1"))
    expect_true(all(is.finite(x$sd) & x$sd > 0))                  # dispersion read off the figure, never the CV floor
    expect_equal(sum(x$obj_weight), 1)                              # repeated measures on the same subjects: 4 rows count as 1
    expect_true(all(x$verified == "digitised"))
    expect_true(all(grepl("digitised", x$location)))
  }
  expect_equal(prof$mean[prof$analyte == "NTX" & prof$stat == "C0.5"], 0.66)   # stated in the text: used, not the plotted value
  expect_equal(prof$mean[prof$analyte == "BN" & prof$stat == "C0.5"], 4.61)
  expect_match(prof$location[prof$analyte == "BN" & prof$stat == "C0.5"], "conflict", ignore.case = TRUE)
  expect_equal(prof$n[prof$stat == "C0.5"], c(37, 37))             # n at 30 min is 37 (Fig 4), not the 52 per-protocol patients
  expect_true(all(prof$n[prof$stat == "C0.25"] == 34))
  expect_true(all(is.na(prof$cv_assumed)))
})

test_that("the Bested 2023 day-21 medians use an assumed CV of 0.8, n not reported, with the curation corrections", {
  d <- rd_data("training", "pk_summary.csv"); b <- d[d$study == "Bested 2023", ]
  ss <- b[b$regimen == "daily", ]
  expect_equal(nrow(ss), 2); expect_true(all(ss$stat == "Css13h"))
  expect_equal(sort(ss$mean), c(0.33, 5.29))
  expect_true(all(is.na(ss$sd))); expect_true(all(ss$cv_assumed == 0.8)); expect_true(all(ss$obj_weight == 1))
  expect_true(all(is.na(ss$n)))
  expect_true(all(grepl("median", ss$location)) && all(grepl("not reported", ss$location)))
  for (r in seq_len(nrow(b))) {
    expect_false(grepl("first treatment period|first-period", b$location[r]))   # it is a crossover
    expect_match(b$location[r], "0.015")                                     # LOQ, Supplemental Digital Content 4
    expect_match(b$location[r], "46 of 52")                                  # population note
  }
  v <- d[grepl("^Verebey", d$study), ]
  expect_true(all(is.na(v$cv_assumed)) && all(v$obj_weight == 1))           # Verebey rows unchanged: CV = SD / mean
})

test_that("pet_occupancy.csv records the number of doses given before each scan", {
  d <- rd_data("test", "pet_occupancy.csv")
  expect_true("n_doses" %in% names(d))
  expect_true(all(d$n_doses >= 1 & d$n_doses == round(d$n_doses)))
  expect_true(all(d$n_doses[d$regimen == "single"] == 1))
  w <- d[d$study == "Weerts 2008", ]
  # R4 fix: n_doses changed 4→3 (internally consistent: doses on days 15/16/17, scan on day 18 at 24 h trough after day 17)
  expect_equal(nrow(w), 1); expect_equal(w$n_doses, 3)
  expect_match(w$location, "days 15-19")
  expect_false(grepl("89.4|90.2|extchecks", w$location))   # R4 fix round 1 (M6): computed by script 03, not quoted from a scratch script
})

# ---- R3a item 7: the Younger floor is a modelling choice, and its provenance row says so ----
test_that("the Younger 2013 baseline floor row is marked secondary and states it is a modelling choice", {
  b <- read_dataset(file.path("..", "..", "data", "params", "baseline_bounds.csv"))
  y <- b[b$trial == "younger2013", ]
  expect_equal(nrow(y), 1)
  expect_equal(c(y$lower_nrs, y$upper_nrs), c(0.5, 10))                     # bound value unchanged
  expect_equal(y$verified, "secondary")
  expect_match(y$location, "modelling choice", ignore.case = TRUE)
  expect_match(y$location, "not an eligibility criterion", ignore.case = TRUE)
  expect_match(y$location, "05b_younger_floor_study.R", fixed = TRUE)
  expect_true(all(b$verified[b$trial != "younger2013"] == "yes"))          # the eligibility rows are unchanged
})

test_that("data/README.md marks pbo_emax_logit as reference-only", {
  rd <- readLines(file.path("..", "..", "data", "README.md"), encoding = "UTF-8")
  ln <- rd[grepl("pbo_emax_logit", rd, fixed = TRUE)]
  expect_true(length(ln) >= 1)
  expect_true(any(grepl("reference only", ln, ignore.case = TRUE)))
})
