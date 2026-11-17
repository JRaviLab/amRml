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

# Keeping legacy endpoint prediction behavior observable without changing grid
# or model selection
.tuningPathNotes <- function(tuned) {
  extracted <- tryCatch(tune::collect_extracts(tuned), error = function(e) NULL)
  notes <- .emptyModelNotes()
  if (is.null(extracted)) return(notes)
  ids <- grep("^id", names(extracted), value = TRUE)
  for (i in seq_len(nrow(extracted))) {
    path <- extracted$.extracts[[i]]$lambda
    if (!.penaltyCovered(path, extracted$penalty[i])) {
      resample <- if (length(ids)) paste(unlist(extracted[i, ids]), collapse = "/") else NA_character_
      notes <- dplyr::bind_rows(notes, tibble::tibble(
        stage = "tune", resample = resample,
        location = paste0("mixture=", extracted$mixture[i],
                          ", penalty=", signif(extracted$penalty[i], 8)),
        type = "warning",
        note = "Requested penalty is outside the fitted glmnet path; endpoint coefficients are used and the candidate remains eligible for compatibility.",
        n = 1L))
    }
  }
  notes
}
