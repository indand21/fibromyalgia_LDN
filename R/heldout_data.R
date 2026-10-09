# Scripted derivations from held-out test files (R3a item 5, review I11). Both files are read through
# read_test_data() by the scripts that use them, so they are hashed under stage T like every other test file.

# INNOVA reports absolute NRS levels (innova_nrs_levels.csv); trials_test.csv carries change_nrs rows that were
# derived by hand as level minus baseline. Recompute every such row and stop if any disagrees with the curated
# value beyond rounding (both levels are given to 2 decimals, so their difference is exact to 2 decimals and the
# curated change may differ from it by at most half a unit in the second decimal). The curated baseline_nrs and
# baseline_sd of each change row must also equal that arm's day-0 level row.
check_innova_change <- function(levels, tt, trial = "innova", tol = 0.005 + 1e-9) {
  lv <- levels[tolower(levels$trial) == tolower(trial), ]
  ch <- tt[tolower(tt$trial) == tolower(trial) & tt$endpoint == "change_nrs", ]
  if (nrow(ch) == 0) stop("check_innova_change: no change_nrs rows for trial '", trial, "' in trials_test.csv")
  get_level <- function(arm, day, col) {
    v <- lv[[col]][lv$arm == arm & lv$day == day]
    if (length(v) != 1 || !is.finite(v))
      stop("check_innova_change: need exactly one finite ", col, " level row for ", trial, " ", arm, " day ", day,
           " in innova_nrs_levels.csv (found ", length(v), ")")
    v
  }
  out <- do.call(rbind, lapply(seq_len(nrow(ch)), function(i) {
    r <- ch[i, ]
    b <- get_level(r$arm, 0, "nrs_mean"); bsd <- get_level(r$arm, 0, "nrs_sd"); l <- get_level(r$arm, r$day, "nrs_mean")
    if (abs(r$baseline_nrs - b) > tol || abs(r$baseline_sd - bsd) > tol)
      stop(sprintf("check_innova_change: %s %s day %g: curated baseline %g (SD %g) differs from the day-0 level row %g (SD %g)",
                   trial, r$arm, r$day, r$baseline_nrs, r$baseline_sd, b, bsd))
    data.frame(trial = r$trial, arm = r$arm, day = r$day, level = l, baseline = b,
               change_derived = round(l - b, 10), change_curated = r$mean)
  }))
  bad <- abs(out$change_derived - out$change_curated) > tol
  if (any(bad))
    stop("check_innova_change: curated change_nrs rows disagree with level minus baseline beyond rounding: ",
         paste(sprintf("%s day %g curated %g vs derived %g", out$arm[bad], out$day[bad], out$change_curated[bad],
                       out$change_derived[bad]), collapse = "; "))
  out$agrees_within_rounding <- TRUE
  out
}

# QST between-group differences (LDN minus placebo) from the FINAL exploratory paper. The `favours` column gives
# each measure's sign convention ("positive difference favours LDN" or "negative difference favours LDN"); the
# summary reports which arm the point estimate favours and whether the 95% CI excludes zero. Descriptive only.
qst_summary <- function(q) {
  need <- c("measure", "diff", "lo95", "hi95", "unit", "favours")
  miss <- setdiff(need, names(q))
  if (length(miss)) stop("qst_summary: missing column(s): ", paste(miss, collapse = ", "))
  sgn <- ifelse(grepl("^positive difference favours LDN", q$favours), 1,
                ifelse(grepl("^negative difference favours LDN", q$favours), -1, NA))
  if (anyNA(sgn)) stop("qst_summary: unrecognised favours text for: ", paste(q$measure[is.na(sgn)], collapse = ", "),
                       " (expected 'positive/negative difference favours LDN')")
  pf <- ifelse(q$diff * sgn > 0, "LDN", ifelse(q$diff * sgn < 0, "placebo", "neither"))
  data.frame(measure = q$measure, diff = q$diff, lo95 = q$lo95, hi95 = q$hi95, unit = q$unit, favours = q$favours,
             point_favours = pf, ci_excludes_zero = q$lo95 > 0 | q$hi95 < 0,
             note = "source-trial between-group difference (complete case), descriptive; not scored by the model",
             source = q$source, location = q$location)
}
