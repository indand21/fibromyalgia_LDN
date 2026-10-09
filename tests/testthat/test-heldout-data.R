# R3a item 5 (review I11): the hand-derived INNOVA change rows and the QST table come from scripted, hashed reads.
dpath <- function(...) file.path("..", "..", "data", ...)

test_that("check_innova_change recomputes the curated INNOVA change rows (level minus baseline)", {
  lv <- read_dataset(dpath("test", "innova_nrs_levels.csv"))
  tt <- read_dataset(dpath("test", "trials_test.csv"))
  r <- check_innova_change(lv, tt)
  expect_equal(nrow(r), 6)                                         # 2 arms x days 90, 180, 365
  expect_setequal(r$day, c(90, 180, 365))
  expect_equal(r$change_derived[r$arm == "ldn_4.5mg" & r$day == 90], 7.33 - 7.66)
  expect_equal(r$change_derived[r$arm == "placebo" & r$day == 90], 6.83 - 7.47)
  expect_true(all(abs(r$change_derived - r$change_curated) <= 0.005 + 1e-9))
  # the unadjusted day-90 difference that decides the primary elpd
  d90 <- r$change_derived[r$arm == "ldn_4.5mg" & r$day == 90] - r$change_derived[r$arm == "placebo" & r$day == 90]
  expect_equal(d90, 0.31, tolerance = 1e-9)
})

test_that("check_innova_change stops when a curated change row disagrees beyond rounding", {
  lv <- read_dataset(dpath("test", "innova_nrs_levels.csv"))
  tt <- read_dataset(dpath("test", "trials_test.csv"))
  bad <- tt; i <- which(bad$trial == "innova" & bad$endpoint == "change_nrs" & bad$day == 180 & bad$arm == "placebo")
  bad$mean[i] <- bad$mean[i] + 0.02
  expect_error(check_innova_change(lv, bad), "disagree")
  bad2 <- tt; j <- which(bad2$trial == "innova" & bad2$endpoint == "change_nrs" & bad2$arm == "placebo")
  bad2$baseline_nrs[j] <- 7.5
  expect_error(check_innova_change(lv, bad2), "baseline")
  expect_error(check_innova_change(lv[lv$day != 90, ], tt), "level")
})

test_that("qst_summary reports each QST difference with its direction and whether the CI excludes zero", {
  q <- read_dataset(dpath("test", "qst_final.csv"))
  s <- qst_summary(q)
  expect_equal(nrow(s), nrow(q))
  expect_true(all(c("measure", "diff", "lo95", "hi95", "unit", "favours", "point_favours", "ci_excludes_zero") %in% names(s)))
  cpm <- s[s$measure == "conditioned_pain_modulation", ]
  expect_true(cpm$ci_excludes_zero); expect_equal(cpm$point_favours, "LDN")
  ts <- s[s$measure == "temporal_summation_of_pain", ]        # negative favours LDN for this measure
  expect_equal(ts$point_favours, "LDN"); expect_false(ts$ci_excludes_zero)
  ppt <- s[s$measure == "pressure_pain_tolerance_left_leg", ]
  expect_equal(ppt$point_favours, "placebo")
  expect_error(qst_summary(transform(q, favours = "unclear")), "favours")
})
