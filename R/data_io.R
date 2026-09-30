PROVENANCE_COLS <- c("source", "location", "verified")

read_dataset <- function(file) {
  d <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
  miss <- setdiff(PROVENANCE_COLS, names(d))
  if (length(miss)) stop(basename(file), " lacks provenance columns: ", paste(miss, collapse = ", "))
  if (any(is.na(d$source) | trimws(d$source) == "")) stop(basename(file), " has rows without a source")
  if (any(grepl("^\\s*literature\\s*$", d$source, ignore.case = TRUE)))
    stop(basename(file), ": 'Literature' is not a source; cite the reference")
  if (!all(d$verified %in% c("yes", "secondary", "no")))
    stop(basename(file), ": verified must be yes/secondary/no")
  d
}

parse_schedule <- function(s) {
  parts <- strsplit(strsplit(s, ";", fixed = TRUE)[[1]], ":", fixed = TRUE)
  out <- data.frame(start_day = as.numeric(vapply(parts, `[`, "", 1)),
                    dose_mg = as.numeric(vapply(parts, `[`, "", 2)))
  if (is.unsorted(out$start_day, strictly = TRUE) || out$start_day[1] != 0)
    stop("schedule start days must begin at 0 and be strictly increasing: ", s)
  out
}

# Convert a concentration constant to nM (uM x 1000, nM x 1, pM / 1000). Stops on an unrecognised unit.
CONC_QUANTITIES <- c("Ki", "Kd", "IC50", "EC50")
NON_CONC_UNITS <- c(koff = "1/h", fu = "fraction")
to_nM <- function(value, unit) {
  f <- c(uM = 1000, nM = 1, pM = 1e-3)
  bad <- setdiff(unique(unit), names(f))
  if (length(bad) || anyNA(unit))
    stop("to_nM: unrecognised concentration unit(s): ", paste(c(bad, if (anyNA(unit)) "NA"), collapse = ", "),
         " (expected uM, nM or pM)")
  value * unname(f[unit])
}

# Pool binding constants for one (analyte, target, quantity) by geometric mean.
# Selects on `analyte` (NTX/BN); prefers organism == "human" rows and falls back to non-human rows
# only when no human row exists (stated in `source`). Every contributing source string is recorded.
pool_binding <- function(bc, analyte, target, quantity) {
  if (!"analyte" %in% names(bc)) stop("binding_constants lacks an `analyte` column")
  r <- bc[bc$analyte == analyte & bc$target == target & bc$quantity %in% quantity, , drop = FALSE]
  if (nrow(r) == 0) return(list(value = NA_real_, unit = NA_character_, n = 0L, organism = NA_character_, source = NA_character_, fallback = FALSE))
  hum <- r$species == "human"
  fallback <- !any(hum)
  use <- if (fallback) r else r[hum, , drop = FALSE]
  # Unit-aware: concentration constants (Ki, Kd, IC50, EC50) are converted to nM before pooling, so a uM row
  # can never be consumed as nM. Non-concentration constants must carry their expected unit (koff 1/h, fu fraction).
  is_conc <- use$quantity %in% CONC_QUANTITIES
  if (any(is_conc) && !all(is_conc))
    stop("pool_binding: cannot pool concentration and non-concentration quantities together (", paste(quantity, collapse = "/"), ")")
  if (all(is_conc)) {
    use$value <- to_nM(use$value, use$unit); use$unit <- "nM"
  } else {
    exp_u <- NON_CONC_UNITS[use$quantity]
    if (anyNA(exp_u) || any(is.na(use$unit) | use$unit != exp_u))
      stop("pool_binding: unrecognised or unexpected unit for ", analyte, " ", target, " ", paste(quantity, collapse = "/"),
           " (found ", paste(unique(use$unit), collapse = ", "), "; expected ", paste(unique(exp_u), collapse = ", "), ")")
  }
  srcs <- paste0(use$source, " [", use$species, "]")
  src <- paste(srcs, collapse = "; ")
  if (nrow(use) > 1 && length(unique(use$value)) == 1)
    src <- paste0(nrow(use), " sources reporting the same value (", signif(use$value[1], 6), " ", use$unit[1], "): ", src)
  if (fallback) src <- paste0("non-human fallback (no human row): ", src)
  list(value = exp(mean(log(use$value))), unit = use$unit[1], n = nrow(use), organism = paste(unique(use$species), collapse = "/"),
       source = src, fallback = fallback)
}
