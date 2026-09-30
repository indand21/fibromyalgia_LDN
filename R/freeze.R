manifest_path <- function() file.path(v5_root(), "output", "results", "freeze_manifest.csv")

sha_file <- function(f) digest::digest(f, algo = "sha256", file = TRUE)

freeze_stage <- function(stage, files, manifest = manifest_path()) {
  stopifnot(length(files) > 0)
  if (!all(file.exists(files))) stop("cannot freeze missing files: ", paste(files[!file.exists(files)], collapse = ", "))
  rows <- data.frame(stage = stage, file = normalizePath(files, winslash = "/"),
                     sha256 = vapply(files, sha_file, ""),
                     frozen_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), row.names = NULL)
  old <- if (file.exists(manifest)) utils::read.csv(manifest, stringsAsFactors = FALSE) else rows[0, ]
  if (stage %in% old$stage)
    stop("Stage '", stage, "' is already frozen; log a CHANGELOG entry and use a new stage name")
  dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rbind(old, rows), manifest, row.names = FALSE)
  invisible(rows)
}

assert_frozen <- function(stage, manifest = manifest_path()) {
  if (!file.exists(manifest)) stop("No freeze manifest: stage '", stage, "' has not been frozen")
  m <- utils::read.csv(manifest, stringsAsFactors = FALSE)
  m <- m[m$stage == stage, ]
  if (nrow(m) == 0) stop("Stage '", stage, "' has not been frozen")
  now <- vapply(m$file, function(f) if (file.exists(f)) sha_file(f) else NA_character_, "")
  bad <- m$file[is.na(now) | now != m$sha256]
  if (length(bad)) stop("Frozen files changed or missing since freeze: ", paste(basename(bad), collapse = ", "))
  invisible(TRUE)
}

# Test files are hashed into the manifest under stage "T" on their FIRST read (after the guarding stage is
# verified) and re-verified on every later read; a changed test file stops the pipeline.
TEST_STAGE <- "T"
pin_test_file <- function(file, manifest = manifest_path()) {
  if (!file.exists(file)) stop("test data file not found: ", file)
  key <- normalizePath(file, winslash = "/"); h <- sha_file(file)
  m <- utils::read.csv(manifest, stringsAsFactors = FALSE)
  prev <- m[m$stage == TEST_STAGE & m$file == key, , drop = FALSE]
  if (nrow(prev) > 1) stop("freeze manifest has duplicate test-file rows for ", basename(file))
  if (nrow(prev) == 1) {
    if (prev$sha256 != h)
      stop("Test data file changed since its first read: ", basename(file), " (manifest stage '", TEST_STAGE, "')")
    return(invisible(TRUE))
  }
  utils::write.csv(rbind(m, data.frame(stage = TEST_STAGE, file = key, sha256 = h,
                                       frozen_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))),
                   manifest, row.names = FALSE)
  invisible(TRUE)
}

read_test_data <- function(file, stage, manifest = manifest_path()) {
  assert_frozen(stage, manifest)
  pin_test_file(file, manifest)
  read_dataset(file)
}
