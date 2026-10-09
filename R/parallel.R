# CPU parallelism helper (controller ruling 2026-10-01: PSOCK clusters with L'Ecuyer RNG streams, no GPU).
# with_cluster(n_jobs, expr_fun) starts a PSOCK cluster of n_jobs workers, points each worker at the model root and sources
# every R/ file there, optionally seeds the workers, calls expr_fun(cl), and ALWAYS stops the cluster on exit (also on error).
default_n_workers <- function() {
  n <- suppressWarnings(parallel::detectCores())
  if (!is.finite(n)) n <- 1
  as.integer(max(1, min(n - 2, 10)))
}

# Save and restore the master's RNG state, kind included. parallel::clusterSetRNGStream() sets
# RNGkind("L'Ecuyer-CMRG") and restores only `.Random.seed`; when `.Random.seed` did not exist (a session that has not yet
# drawn a random number) it REMOVES it, leaving the L'Ecuyer kind in force, so every later set.seed() silently draws from a
# different generator. R1 fix round 4: that is why the round-3 occupancy draws were generated under L'Ecuyer-CMRG.
rng_state <- function() list(kind = RNGkind(),
                             seed = if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
                               get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL)
rng_restore <- function(st) {
  do.call(RNGkind, as.list(st$kind))                 # resets .Random.seed, so restore it afterwards
  if (is.null(st$seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", st$seed, envir = globalenv())
  invisible(st)
}

with_cluster <- function(n_jobs = default_n_workers(), expr_fun, seed = NULL, root = v5_root()) {
  stopifnot(is.function(expr_fun), n_jobs >= 1)
  st <- rng_state()
  on.exit(rng_restore(st), add = TRUE)
  cl <- parallel::makePSOCKcluster(as.integer(n_jobs))
  on.exit(parallel::stopCluster(cl), add = TRUE, after = FALSE)
  parallel::clusterCall(cl, function(r) {
    options(v5.root = r)
    for (f in list.files(file.path(r, "R"), pattern = "\\.R$", full.names = TRUE)) source(f, local = globalenv())
    invisible(TRUE)
  }, normalizePath(root, winslash = "/"))
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
  expr_fun(cl)
}

# Per-TASK L'Ecuyer streams. clusterSetRNGStream() gives one stream per WORKER, so a job's results depend on how many
# workers ran it and on how tasks were scheduled onto them. task_seeds(seed, n) returns one independent stream per task;
# with_task_seed() installs one for the duration of a task. Results then depend only on `seed` and the task index.
task_seeds <- function(seed, n_tasks) {
  stopifnot(n_tasks >= 1)
  st <- rng_state()
  on.exit(rng_restore(st), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  out <- vector("list", n_tasks)
  out[[1]] <- .Random.seed
  for (i in seq_len(n_tasks - 1L)) out[[i + 1L]] <- parallel::nextRNGStream(out[[i]])
  out
}

with_task_seed <- function(stream, expr) {
  st <- rng_state()
  on.exit(rng_restore(st), add = TRUE)
  assign(".Random.seed", stream, envir = globalenv())
  force(expr)
}

# Per-task seeding for a task whose numbers must match an EXISTING sequential run produced with set.seed() under
# the default generator (R2 item 5: the parallel MCMC of script 04 must reproduce the sequential fit bit for bit).
# task_seeds()/with_task_seed() give independent L'Ecuyer streams, which is the right default for new work but
# would renumber such a run; with_mt_seed() installs an explicit Mersenne-Twister / Inversion / Rejection stream from `seed` instead, so
# the task depends only on `seed` -- never on the worker, the worker count or the ambient RNG state and kind.
with_mt_seed <- function(seed, expr) {
  st <- rng_state()
  on.exit(rng_restore(st), add = TRUE)
  # all three kinds are pinned, so neither normal.kind nor sample.kind can leak in from the caller (fix round 1)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  force(expr)
}
