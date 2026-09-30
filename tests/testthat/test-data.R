test_that("read_dataset rejects missing or vague provenance", {
  f <- tempfile(fileext = ".csv")
  write.csv(data.frame(x = 1, source = "Literature", location = "-", verified = "yes"), f, row.names = FALSE)
  expect_error(read_dataset(f), "not a source")
  write.csv(data.frame(x = 1, source = "A 2020", location = "T1", verified = "maybe"), f, row.names = FALSE)
  expect_error(read_dataset(f), "verified")
  write.csv(data.frame(x = 1), f, row.names = FALSE)
  expect_error(read_dataset(f), "provenance")
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
  expect_true(all(d$stat %in% c("Cmax", "Tmax", "AUCinf", "thalf")))
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
