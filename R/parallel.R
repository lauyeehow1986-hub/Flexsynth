# Internal: optional parallel generation of the `m` independent synthetic
# replicates. The replicates are statistically independent, so producing them on
# several workers is embarrassingly parallel. Serial synthesis (the default,
# `parallel = FALSE`) is byte-for-byte unchanged; the parallel path is opt-in via
# synth_control(parallel = ...) and used by both synth() and synth_linked().

# Normalise the `parallel` control to a worker count (1 == serial). TRUE means
# "use the machine", read from getOption("mc.cores") or detectCores(); a positive
# integer is taken verbatim (validated in synth_control()).
parallel_workers <- function(parallel) {
  if (is.null(parallel) || isFALSE(parallel)) return(1L)
  if (isTRUE(parallel)) {
    n <- getOption("mc.cores", parallel::detectCores(logical = TRUE))
    if (length(n) != 1L || is.na(n) || n < 1L) n <- 1L
    return(max(1L, as.integer(n)))
  }
  max(1L, as.integer(parallel))
}

# Run `gen_fun()` `m` times and return the list of results. When the control asks
# for more than one worker and there is more than one replicate to make, the work
# is spread over a PSOCK cluster; otherwise it runs serially in-process (the exact
# current code path). `seed` (when supplied) seeds independent L'Ecuyer streams so
# a parallel run is reproducible for a fixed (seed, worker count) — the streams
# differ from the single serial stream, which is inherent to parallel RNG, so a
# parallel result is not expected to equal the serial one.
run_replicates <- function(m, gen_fun, control, seed) {
  workers <- min(parallel_workers(control$parallel), m)
  if (m == 1L || workers <= 1L) {
    return(lapply(seq_len(m), function(i) gen_fun()))
  }

  cl <- tryCatch(parallel::makePSOCKcluster(workers), error = function(e) NULL)
  if (is.null(cl)) {
    warning("could not start a parallel cluster; running serially.",
            call. = FALSE)
    return(lapply(seq_len(m), function(i) gen_fun()))
  }
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Each PSOCK worker is a fresh R session; the synthesis closure resolves the
  # engine's internal functions through the flexsynth namespace, so the package
  # must be loadable on the workers. If it is not (e.g. an uninstalled dev tree),
  # fall back to serial rather than error.
  loadable <- tryCatch(
    all(unlist(parallel::clusterEvalQ(
      cl, requireNamespace("flexsynth", quietly = TRUE)))),
    error = function(e) FALSE
  )
  if (!isTRUE(loadable)) {
    warning("flexsynth is not loadable on the parallel workers; running serially.",
            call. = FALSE)
    return(lapply(seq_len(m), function(i) gen_fun()))
  }

  parallel::clusterSetRNGStream(cl, iseed = seed)
  parallel::parLapply(cl, seq_len(m), function(i) gen_fun())
}
