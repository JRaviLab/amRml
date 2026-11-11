# Shared CV execution for binary and multinomial models. Tuning, recipes,
# parameter grids and fold assignment are identical in both execution modes.
.tuneFold <- function(i, workflow, resamples, grid) {
  withr::local_options(list(cli.progress_handlers_only = character()))
  fold <- rsample::manual_rset(resamples$splits[i], ids = resamples$id[i])
  suppressMessages(suppressWarnings(tune::tune_grid(
    workflow, resamples = fold, grid = grid,
    metrics = yardstick::metric_set(yardstick::mcc),
    control = tune::control_grid(allow_par = FALSE, verbose = FALSE,
                                extract = .extractPenaltyPath))))
}

.extractPenaltyPath <- function(fit) {
  list(lambda = parsnip::extract_fit_engine(fit)$lambda)
}

.tuneModelFolds <- function(workflow, resamples, grid, workers = 1L) {
  .checkArgCount(workers, "cv_workers")
  workers <- min(workers, nrow(resamples))
  if (workers == 1L) {
    return(suppressMessages(suppressWarnings(tune::tune_grid(
      workflow, resamples = resamples, grid = grid,
      metrics = yardstick::metric_set(yardstick::mcc),
      control = tune::control_grid(allow_par = FALSE, verbose = FALSE,
                                  extract = .extractPenaltyPath)))))
  }
  withr::local_envvar(c(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1",
    R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)))
  backend <- BiocParallel::SnowParam(workers = workers, type = "SOCK",
                                     exportglobals = FALSE, progressbar = FALSE)
  on.exit(BiocParallel::bpstop(backend), add = TRUE)
  parts <- BiocParallel::bplapply(seq_len(nrow(resamples)), .tuneFold,
    workflow = workflow, resamples = resamples, grid = grid, BPPARAM = backend)
  # Identical grid/config IDs and outcome/metric attributes in every fold.
  out <- do.call(rbind, parts)
  attributes_to_keep <- attributes(parts[[1]])
  for (name in setdiff(names(attributes_to_keep), c("names", "row.names", "class"))) {
    attr(out, name) <- attributes_to_keep[[name]]
  }
  out
}

.penaltyCovered <- function(path, penalty) {
  if (!length(path) || any(!is.finite(path))) return(rep(FALSE, length(penalty)))
  penalty >= min(path) * (1 - 1e-10) & penalty <= max(path) * (1 + 1e-10)
}

# tune uses submodels to evaluate many penalties from one fitted path. Reject
# configurations whose requested penalty was never reached in any one fold,
# rather than allowing glmnet to silently substitute endpoint coefficients.
.checkTuningPaths <- function(tuned) {
  extracted <- tryCatch(tune::collect_extracts(tuned), error = function(e) NULL)
  notes <- .emptyModelNotes()
  if (is.null(extracted)) return(list(tuned = tuned, notes = notes))
  invalid <- character()
  ids <- grep("^id", names(extracted), value = TRUE)
  for (i in seq_len(nrow(extracted))) {
    path <- extracted$.extracts[[i]]$lambda
    if (!.penaltyCovered(path, extracted$penalty[i])) {
      invalid <- c(invalid, extracted$.config[i])
      resample <- if (length(ids)) paste(unlist(extracted[i, ids]), collapse = "/") else NA_character_
      notes <- dplyr::bind_rows(notes, tibble::tibble(
        stage = "tune", resample = resample,
        location = paste0("mixture=", extracted$mixture[i],
                          ", penalty=", signif(extracted$penalty[i], 8)),
        type = "warning",
        note = "Requested penalty is outside the fitted glmnet path; candidate excluded from selection.",
        n = 1L))
    }
  }
  invalid <- unique(invalid)
  if (length(invalid)) {
    tuned$.metrics <- lapply(tuned$.metrics, function(metrics) {
      metrics$.estimate[metrics$.config %in% invalid] <- NA_real_
      metrics
    })
  }
  list(tuned = tuned, notes = notes)
}

# Extend only a final path that misses the selected CV penalty. Dense bridge
# points preserve warm starts; this is not a different search grid.
.extendPenaltyPath <- function(path, penalty) {
  if (!length(path) || any(!is.finite(path)) || any(path <= 0)) {
    stop("Cannot extend an invalid glmnet penalty path.")
  }
  endpoint <- if (penalty < min(path)) min(path) else max(path)
  n <- max(2L, ceiling(abs(log10(endpoint / penalty)) * 50) + 1L)
  bridge <- 10^seq(log10(endpoint), log10(penalty), length.out = n)
  sort(unique(signif(c(path, bridge, penalty), 14)), decreasing = TRUE)
}
