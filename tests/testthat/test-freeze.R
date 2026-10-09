make_csv <- function(dir, name = "d.csv") {
  f <- file.path(dir, name)
  write.csv(data.frame(x = 1, source = "Rabiner 2011", location = "Fig 1", verified = "yes"), f, row.names = FALSE)
  f
}

test_that("read_test_data refuses before the stage is frozen", {
  d <- tempfile(); dir.create(d); f <- make_csv(d); m <- file.path(d, "manifest.csv")
  expect_error(read_test_data(f, "A", manifest = m), "has not been frozen|No freeze manifest")
})

test_that("read_test_data works after freeze and refuses if a frozen file changes", {
  d <- tempfile(); dir.create(d); f <- make_csv(d); train <- make_csv(d, "train.csv")
  m <- file.path(d, "manifest.csv")
  freeze_stage("A", train, manifest = m)
  expect_equal(nrow(read_test_data(f, "A", manifest = m)), 1)
  cat("tamper\n", file = train, append = TRUE)
  expect_error(read_test_data(f, "A", manifest = m), "changed or missing")
})

test_that("a stage cannot be frozen twice", {
  d <- tempfile(); dir.create(d); train <- make_csv(d, "train.csv"); m <- file.path(d, "manifest.csv")
  freeze_stage("A", train, manifest = m)
  expect_error(freeze_stage("A", train, manifest = m), "already frozen")
})

test_that("training scripts never reference data/test", {
  scripts <- list.files(file.path("..", "..", "scripts"), pattern = "^0[1-5]_.*\\.R$", full.names = TRUE)
  scripts <- scripts[!grepl("^0[1-5]_test_", basename(scripts))]  # held-out test scripts legitimately read data/test (via read_test_data)
  for (s in scripts) expect_false(any(grepl("data/test|data\\\\test", readLines(s))), info = s)
})

test_that("PET script reads test data only through read_test_data", {
  s <- readLines(file.path("..", "..", "scripts", "03_test_pet.R"))
  expect_false(any(grepl("read_dataset\\(\"data/test", s)))
  expect_true(any(grepl("read_test_data\\(\"data/test/pet_occupancy.csv\", \"A\"\\)", s)))
})

test_that("trial test script reads test data only through read_test_data at stage B", {
  s <- readLines(file.path("..", "..", "scripts", "06_test_trials.R"))
  expect_false(any(grepl("read_dataset\\(\"data/test", s)))
  expect_false(any(grepl("read\\.csv|readRDS", s) & grepl("data/test", s)))
  expect_true(any(grepl("read_test_data\\(\"data/test/trials_test.csv\", \"B\"\\)", s)))
})

test_that("virtual-trial script reads test data only through read_test_data at stage B", {
  s <- readLines(file.path("..", "..", "scripts", "07_virtual_trials.R"))
  expect_false(any(grepl("read_dataset\\(\"data/test", s)))
  expect_false(any(grepl("read\\.csv|readRDS", s) & grepl("data/test", s)))
  expect_true(any(grepl("read_test_data\\(\"data/test/trials_test.csv\", \"B\"\\)", s)))
})

test_that("cytokine test script reads test data only through read_test_data at stage B", {
  s <- readLines(file.path("..", "..", "scripts", "08_test_cytokines.R"))
  expect_false(any(grepl("read_dataset\\(\"data/test", s)))
  expect_false(any(grepl("read\\.csv|readRDS", s) & grepl("data/test", s)))
  expect_true(any(grepl("read_test_data\\(\"data/test/cytokines_parkitny.csv\", \"B\"\\)", s)))
})

test_that("scripts 06-11 touch data/test only through read_test_data, and script 11 never does", {
  sd <- file.path("..", "..", "scripts")
  scripts <- list.files(sd, pattern = "^(0[6-9]|1[01])_.*[.]R$", full.names = TRUE)
  expect_gte(length(scripts), 6)
  for (s in scripts) {
    l <- readLines(s); hit <- grepl("data/test", l, fixed = TRUE) | grepl("data\\test", l, fixed = TRUE)
    if (grepl("^11_", basename(s))) expect_false(any(hit), info = paste(basename(s), "must not reference data/test at all"))
    else expect_true(all(grepl("read_test_data(", l[hit], fixed = TRUE)), info = paste(basename(s), "references data/test outside read_test_data("))
  }
})

test_that("read_test_data pins the test file hash on first read and refuses a changed file", {
  d <- tempfile(); dir.create(d); f <- make_csv(d); train <- make_csv(d, "train.csv"); m <- file.path(d, "manifest.csv")
  freeze_stage("A", train, manifest = m)
  expect_equal(nrow(read_test_data(f, "A", manifest = m)), 1)
  mm <- read.csv(m); t_rows <- mm[mm$stage == "T", ]
  expect_equal(nrow(t_rows), 1); expect_equal(t_rows$sha256, sha_file(f))
  expect_equal(nrow(read_test_data(f, "A", manifest = m)), 1)          # unchanged: fine, no duplicate row
  expect_equal(sum(read.csv(m)$stage == "T"), 1)
  write.csv(data.frame(x = 2, source = "Rabiner 2011", location = "Fig 1", verified = "yes"), f, row.names = FALSE)
  expect_error(read_test_data(f, "A", manifest = m), "Test data file changed")
})

test_that("pinning a test file does not block later freezing of other stages", {
  d <- tempfile(); dir.create(d); f <- make_csv(d); train <- make_csv(d, "train.csv"); t2 <- make_csv(d, "t2.csv")
  m <- file.path(d, "manifest.csv")
  freeze_stage("A", train, manifest = m); read_test_data(f, "A", manifest = m)
  expect_silent(freeze_stage("B", t2, manifest = m))
})

test_that("PET script reads the Rabiner ED50 only through read_test_data (so it is hashed under stage T)", {
  s <- readLines(file.path("..", "..", "scripts", "03_test_pet.R"))
  hit <- grepl("pet_ed50", s, fixed = TRUE)
  expect_true(any(hit))
  expect_true(all(grepl("read_test_data(\"data/test/pet_ed50.csv\", \"A\")", s[hit], fixed = TRUE)))
})

test_that("PET script adopts stage J and freezes it with its inputs and outputs", {
  s <- paste(readLines(file.path("..", "..", "scripts", "03_test_pet.R")), collapse = "\n")
  expect_match(s, "freeze_stage(\"J\"", fixed = TRUE)
  expect_match(s, "binding_params_J.rds", fixed = TRUE)
  expect_match(s, "steady_state_occupancy_by_stage.csv", fixed = TRUE)
  expect_match(s, "occupancy_uncertainty_J.csv", fixed = TRUE)
})

test_that("PET script writes the keq range and the PK sensitivity, and freezes the sensitivity inputs with stage J", {
  s <- paste(readLines(file.path("..", "..", "scripts", "03_test_pet.R")), collapse = "\n")
  expect_match(s, "keq_sensitivity", fixed = TRUE)
  expect_match(s, "occupancy_sensitivity_J.csv", fixed = TRUE)
  expect_match(s, "\"data/training/pk_summary_sensitivity_terminal.csv\", sens_files))", fixed = TRUE)
  p1 <- paste(readLines(file.path("..", "..", "scripts", "01_fit_pk.R")), collapse = "\n")
  expect_match(p1, "FP_FIXED <- 0.95", fixed = TRUE)
  expect_false(grepl("rbind(obs, read_dataset(\"data/training/pk_summary_sensitivity_terminal.csv\"))\nsel <-", p1, fixed = TRUE))
})

# ---- R3a item 5 (review I11): the INNOVA level file and the QST file are read through read_test_data ----
test_that("script 06 re-derives the INNOVA change rows from the level file read through read_test_data", {
  s <- readLines(file.path("..", "..", "scripts", "06_test_trials.R"))
  expect_true(any(grepl('read_test_data("data/test/innova_nrs_levels.csv", "B")', s, fixed = TRUE)))
  expect_true(any(grepl("check_innova_change(", s, fixed = TRUE)))
})

test_that("script 08 reports the QST table from qst_final.csv read through read_test_data", {
  s <- readLines(file.path("..", "..", "scripts", "08_test_cytokines.R"))
  expect_true(any(grepl('read_test_data("data/test/qst_final.csv", "B")', s, fixed = TRUE)))
  expect_true(any(grepl("qst_summary(", s, fixed = TRUE)))
})

# ---- R3a item 7: the Younger floor study is a pipeline script and reads no test data ----
test_that("the Younger floor study script exists, reads no test data and uses the pipeline's own functions", {
  f <- file.path("..", "..", "scripts", "05b_younger_floor_study.R")
  expect_true(file.exists(f))
  s <- readLines(f)
  expect_false(any(grepl("data/test|data\\\\test|read_test_data", s)))
  txt <- paste(s, collapse = "\n")
  for (fn in c("calibrate_omegas(", "match_trunc_parent(", "younger_floor_stability.csv"))
    expect_true(grepl(fn, txt, fixed = TRUE), info = fn)
  expect_false(grepl("set.seed(", txt, fixed = TRUE))           # no hand-duplicated RNG order
})


# ---- R3a fix round 1: post-hoc analyses are labelled exploratory (global constraints) ----
test_that("the responder-scoring definition text in script 07 labels it exploratory", {
  s <- paste(readLines(file.path("..", "..", "scripts", "07_virtual_trials.R")), collapse = "\n")
  expect_true(grepl("EXPLORATORY (added after the held-out results were seen", s, fixed = TRUE))
})
