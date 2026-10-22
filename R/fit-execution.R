# Validate and write resumable fits

# An explicit worker receives one group and reads its feature stores locally.
.runTaskGroup <- function(job, species, model_dir, save_fits, cache) {
  shared_fit <- NULL
  outcomes <- vector("list", nrow(job$tasks))
  for (j in seq_len(nrow(job$tasks))) {
    task <- job$tasks[j, ]
    keep_fit <- nrow(job$tasks) > 1L && is.null(shared_fit)
    outcome <- tryCatch(
      withr::with_seed(
        task$seed,
        .fitTask(
          task = task, path = job$paths[[j]],
          species = species,
          rows = job$rows[[j]], columns = job$columns[[j]],
          encoding = job$encoding[[j]], cache = cache,
          fitted = shared_fit, return_fitted = keep_fit
        ),
        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
        .rng_sample_kind = "Rejection"
      ),
      error = function(e) list(error = conditionMessage(e))
    )
    if (is.null(outcome[["error"]])) {
      if (keep_fit) shared_fit <- outcome$fitted
      outcome$fitted <- NULL
      if (save_fits) {
        .saveFitCheckpoint(outcome$model, file.path(model_dir, paste0(task$task_id, ".rds")))
      }
      outcome$model <- NULL
    }
    outcomes[[j]] <- outcome
  }
  # Commit the group only after any requested model files have been saved.
  .saveFitCheckpoint(outcomes, job$checkpoint)
  TRUE
}

# Only pending groups are sent to workers; no feature matrix is serialized.
.fitGroupJob <- function(group, listed, matrices, index, out_dir, checkpoint) {
  positions <- match(listed$matrix_id[group], matrices$matrix_id)
  list(
    tasks = listed[group, ],
    paths = file.path(out_dir, matrices$file[positions]),
    rows = index$rows[index$matrix_scopes[positions]],
    columns = index$columns[positions],
    encoding = matrices$encoding[positions],
    checkpoint = checkpoint
  )
}

# Describe the tasks and matrix files a checkpoint belongs to. Checksums prevent
# completed fits from being mixed with changed training inputs when resuming.
## This is draft functionality, and should be replaced with proper ORB manifest
## logging after verification, but modeling resume functionality is badly needed.
.fitCheckpointManifest <- function(tasks, groups, reuse_fits, save_fits) {
  files <- tasks$built$matrices$file[
    match(tasks$tasks$matrix_id, tasks$built$matrices$matrix_id)
  ]
  inputs <- unique(c("matrix_index.rds", files))
  checksums <- tools::md5sum(file.path(tasks$built$out_dir, inputs))
  names(checksums) <- inputs
  if (anyNA(checksums)) {
    stop("Cannot read the matrix files needed for checkpoints.", call. = FALSE)
  }
  list(
    version = 2L,
    contract = "saved-models",
    tasks = tasks$tasks,
    matrices = tasks$built$matrices,
    groups = groups,
    reuse_fits = reuse_fits,
    save_fits = save_fits,
    species = tasks$built$eligibility$profile$orb$dataset_id,
    input_checksums = checksums
  )
}

# Write one completed task group before its worker finishes. A temporary file
# prevents an interrupted write from appearing as a finished checkpoint.
.saveFitCheckpoint <- function(result, path) {
  temporary <- tempfile(".fit-", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(result, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) {
    stop("Could not save fit checkpoint: ", path, call. = FALSE)
  }
  invisible(path)
}

# Fail invalid checkpoints rather than treating them as model fits
.readFitCheckpoint <- function(path, task_ids, model_dir, save_fits) {
  result <- readRDS(path)
  valid <- is.list(result) && length(result) == length(task_ids)
  if (valid) for (j in seq_along(task_ids)) {
    x <- result[[j]]
    failure <- is.list(x) && identical(names(x), "error") &&
      is.character(x$error) && length(x$error) == 1L && !is.na(x$error)
    success <- is.list(x) &&
      all(c("performance", "importance", "predictions", "undefined") %in% names(x)) &&
      is.data.frame(x$performance) && is.data.frame(x$importance) &&
      is.data.frame(x$predictions) &&
      identical(as.character(x$performance$task_id), task_ids[[j]]) &&
      all(x$importance$task_id == task_ids[[j]]) &&
      all(x$predictions$task_id == task_ids[[j]]) &&
      is.character(x$undefined) && length(x$undefined) == 1L
    if (!failure && !success) valid <- FALSE
    if (success && save_fits &&
        !file.exists(file.path(model_dir, paste0(task_ids[[j]], ".rds")))) {
      stop("A saved model is missing for checkpoint task: ", task_ids[[j]], call. = FALSE)
    }
  }
  if (!valid) stop("Invalid fitted task checkpoint: ", path, call. = FALSE)
  result
}

# Empty schemas keep failed fits usable
.emptyFitResults <- function() {
  list(
    performance = tibble::tibble(
      task_id = character(), species = character(), num_obs = integer(), res_prop = numeric(),
      n_test = integer(), res_prop_test = numeric(), n_feat = integer(),
      model = character(), n_feats_returned = integer(), n_fold = integer(),
      fit_penalty = numeric(), fit_mixture = numeric(), mcc = numeric(), nmcc = numeric(),
      spec = numeric(), sens = numeric(), log2_apop = numeric(), f1 = numeric(),
      bal_acc = numeric(), auprc = numeric(), run_time_sec = numeric(),
      seed = integer(), date = character()
    ),
    importance = tibble::tibble(
      task_id = character(), species = character(), Variable = character(),
      Importance = numeric(), Sign = character()
    ),
    predictions = tibble::tibble(
      task_id = character(), species = character(), genome_id = character(), truth = character(),
      predicted = character(), prob_resistant = numeric()
    )
  )
}

# Restore previous files on failure
.restoreFitResults <- function(staging, out_dir, produced) {
  backup <- tempfile(".results-previous-", tmpdir = out_dir)
  dir.create(backup)
  old_moved <- character()
  new_moved <- character()
  tryCatch({
    for (name in produced) {
      if (file.exists(file.path(out_dir, name))) {
        if (!file.rename(file.path(out_dir, name), file.path(backup, name))) {
          stop("Could not move previous result: ", name)
        }
        old_moved <- c(old_moved, name)
      }
      if (!file.rename(file.path(staging, name), file.path(out_dir, name))) {
        stop("Could not install result: ", name)
      }
      new_moved <- c(new_moved, name)
    }
  }, error = function(e) {
    for (name in rev(new_moved)) unlink(file.path(out_dir, name), recursive = TRUE)
    restored <- vapply(rev(old_moved), function(name) {
      file.rename(file.path(backup, name), file.path(out_dir, name))
    }, logical(1))
    if (all(restored)) unlink(backup, recursive = TRUE)
    stop(conditionMessage(e), if (!all(restored)) paste0(
      " Previous results may be recoverable from: ", backup
    ), call. = FALSE)
  })
  unlink(backup, recursive = TRUE)
  invisible(out_dir)
}

# A completed fit can be reopened to finish an export without refitting
.readCompletedFit <- function(tasks, manifest, call) {
  out_dir <- tasks$built$out_dir
  path <- file.path(out_dir, "results", "fit-completion.rds")
  if (!file.exists(path)) {
    .amrAbort("invalid_argument",
              "No compatible checkpoints or completed fit were found; cannot resume.",
              call = call)
  }
  receipt <- readRDS(path)
  if (!is.list(receipt) || !identical(receipt$manifest, manifest) ||
      !is.character(receipt$checksums) || !length(receipt$checksums)) {
    .amrAbort("invalid_argument", "The completed fit does not match these inputs.", call = call)
  }
  observed <- tools::md5sum(file.path(out_dir, names(receipt$checksums)))
  names(observed) <- names(receipt$checksums)
  if (!identical(observed, receipt$checksums)) {
    .amrAbort("invalid_argument", "Completed fit outputs are missing or changed; cannot resume.",
              call = call)
  }
  results <- lapply(c("performance", "importance", "predictions"), function(name) {
    arrow::read_parquet(file.path(out_dir, "results", paste0(name, ".parquet")))
  })
  names(results) <- c("performance", "importance", "predictions")
  structure(c(list(tasks = arrow::read_parquet(file.path(out_dir, "tasks.parquet"))),
              results, list(built = tasks$built)), class = "amr_fit")
}

#' Fit the listed models
#'
#' Fits tasks from [expandTasks()] in parallel with old amRml's method:
#' glmnet logistic regression over 10 penalties and 6 mixtures, tuned by
#' stratified `n_fold` cross-validation and chosen by MCC, then scored on the
#' task's test genomes.
#'
#' @details
#' Cross-drug tasks with exactly matching training inputs reuse one fitted model;
#' set `reuse_fits = FALSE` to disable reuse for comparison. Reused tasks'
#' `run_time_sec` excludes shared training time for follower tasks; the first
#' task in a group includes the training cost.
#' Each completed task group is saved in `fit-checkpoints/` so it can be reused
#' after an interruption. Checkpoints are removed after results are written.
#' Use `resume = TRUE` with the same tasks and output directory to reuse them.
#' Each task runs under its own seed. A task that errors is recorded as
#' `fit_failed`, and one with undefined metrics as `metrics_undefined`; the rest
#' carry on. A shuffled task permutes the training matrix's labels before the
#' split.
#'
#' Completed results can also be reopened with `resume = TRUE` without fitting
#' again. Version 1 checkpoints require the original implementation or a fresh
#' run in another output directory.
#'
#' @param tasks An `amr_tasks` from [expandTasks()].
#' @param save_fits Save each fitted model to `results/models/<task_id>.rds`:
#'   the fitted workflow (`fit`) and the `features` its predictors `x1`, `x2`,
#'   ... stand for.
#'   Must match the interrupted run when resuming.
#' @param overwrite Replace results already in the output folder.
#' @param reuse_fits Reuse fitted models among compatible cross-drug tasks.
#' @param resume Reuse saved task groups from an interrupted fit or reopen an
#'   unchanged completed fit. The task plan, matrix files and `reuse_fits`
#'   setting must match the original run.
#'   `save_fits` must also match.
#'
#' @param BPPARAM A `BiocParallelParam`, e.g. `BiocParallel::MulticoreParam(4)`,
#'   or `BiocParallel::SerialParam()` to fit one task at a time.
#' @return An `amr_fit` list, also written to the output folder: `tasks` (with
#'   `status`, `rule_id` and `message`; `tasks.parquet`), `performance`,
#'   `importance` and `predictions` (`results/`), and the `built` matrices. In
#'   `importance`, `Sign` "POS" points to Susceptible and "NEG" to Resistant. The
#'   run's settings and task counts go to `run.json`.
#'
#' @section Errors:
#' Invalid arguments or incompatible resume inputs raise `amrml_invalid_argument`.
#' File and checkpoint errors stop the run without removing completed checkpoints.
#'
#' @examples
#' \dontrun{
#' fit <- fitModels(expandTasks(built))
#' fit$performance
#' }
#' @export
fitModels <- function(tasks, save_fits = TRUE,
                      BPPARAM = BiocParallel::bpparam(), overwrite = FALSE,
                      reuse_fits = TRUE, resume = FALSE) {
  call <- rlang::current_env()
  .checkArgClass(tasks, "tasks", "amr_tasks", "expandTasks()", call = call)
  .checkArgClass(BPPARAM, "BPPARAM", "BiocParallelParam", "BiocParallel", call = call)
  .checkArgFlag(overwrite, "overwrite", call = call)
  .checkArgFlag(reuse_fits, "reuse_fits", call = call)
  .checkArgFlag(resume, "resume", call = call)
  .checkArgFlag(save_fits, "save_fits", call = call)
  if (resume && overwrite) {
    .amrAbort("invalid_argument", "`resume` and `overwrite` cannot both be TRUE.", call = call)
  }
  out_dir <- tasks$built$out_dir
  produced <- c("tasks.parquet", "results", "run.json")
  if (any(file.exists(file.path(out_dir, produced))) && !overwrite && !resume) {
    .amrAbort("invalid_argument",
              "The output folder already has results; set `overwrite = TRUE` to replace them.",
              observed = list(out_dir = out_dir), call = call)
  }

  started <- Sys.time()
  species <- tasks$built$eligibility$profile$orb$dataset_id
  listed <- tasks$tasks
  listed$status <- rep("fitted", nrow(listed))
  listed$rule_id <- rep(NA_character_, nrow(listed))
  listed$message <- rep(NA_character_, nrow(listed))
  index <- readRDS(file.path(out_dir, "matrix_index.rds"))
  groups <- if (reuse_fits) {
    .fitTaskGroups(listed, tasks$built$matrices, index)
  } else {
    as.list(seq_len(nrow(listed)))
  }
  manifest <- .fitCheckpointManifest(tasks, groups, reuse_fits, save_fits)
  checkpoint_dir <- file.path(out_dir, "fit-checkpoints")
  manifest_path <- file.path(checkpoint_dir, "manifest.rds")
  if (resume) {
    if (!file.exists(manifest_path)) return(.readCompletedFit(tasks, manifest, call))
    if (!identical(readRDS(manifest_path), manifest)) {
      .amrAbort("invalid_argument",
                "Saved fits do not match these tasks, settings or matrix files; cannot resume.",
                call = call)
    }
  } else {
    if (dir.exists(checkpoint_dir) && !overwrite) {
      .amrAbort("invalid_argument",
                "Task checkpoints already exist; use `resume = TRUE` to reuse them.", call = call)
    }
    unlink(checkpoint_dir, recursive = TRUE)
    dir.create(checkpoint_dir, recursive = TRUE)
    .saveFitCheckpoint(manifest, manifest_path)
  }
  model_dir <- file.path(checkpoint_dir, "models")
  if (save_fits) dir.create(model_dir, showWarnings = FALSE)
  group_paths <- file.path(checkpoint_dir, sprintf("group-%06d.rds", seq_along(groups)))
  recovered <- which(file.exists(group_paths))
  # Validate recovered groups before scheduling more work.
  for (g in recovered) {
    .readFitCheckpoint(group_paths[[g]], listed$task_id[groups[[g]]], model_dir, save_fits)
  }
  pending <- which(!file.exists(group_paths))
  if (resume) message(length(recovered), " task groups recovered; ",
                      length(pending), " task groups remaining.")
  jobs <- lapply(pending, function(g) {
    .fitGroupJob(groups[[g]], listed, tasks$built$matrices, index, out_dir, group_paths[[g]])
  })
  # An empty cache is sent once; each worker populates its own copy across groups.
  cache <- new.env(parent = emptyenv())
  if (length(jobs)) BiocParallel::bplapply(
    jobs, .runTaskGroup, cache = cache, species = species, model_dir = model_dir, save_fits = save_fits,
    BPPARAM = BPPARAM
  )

  outcomes <- vector("list", nrow(listed))
  for (g in seq_along(groups)) {
    outcomes[groups[[g]]] <- .readFitCheckpoint(
      group_paths[[g]], listed$task_id[groups[[g]]], model_dir, save_fits
    )
  }
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[[i]]
    if (!is.null(outcome[["error"]])) {
      listed$status[[i]] <- "failed"
      listed$rule_id[[i]] <- "fit_failed"
      listed$message[[i]] <- outcome$error
      outcomes[i] <- list(NULL)
    } else if (!is.na(outcome$undefined)) {
      listed$rule_id[[i]] <- "metrics_undefined"
      listed$message[[i]] <- paste("Undefined:", outcome$undefined)
    }
  }
  empty <- .emptyFitResults()
  tables <- lapply(names(empty), function(name) {
    dplyr::bind_rows(c(list(empty[[name]]), lapply(outcomes, `[[`, name)))
  })
  names(tables) <- names(empty)
  fit <- structure(c(list(tasks = listed), tables, list(built = tasks$built)), class = "amr_fit")

  staging <- tempfile(".results-building-", tmpdir = out_dir)
  dir.create(file.path(staging, "results"), recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  arrow::write_parquet(fit$tasks, file.path(staging, "tasks.parquet"))
  for (name in names(empty)) {
    arrow::write_parquet(fit[[name]], file.path(staging, "results", paste0(name, ".parquet")))
  }
  if (save_fits) {
    destination <- file.path(staging, "results", "models")
    dir.create(destination)
    ids <- listed$task_id[listed$status == "fitted"]
    source <- file.path(model_dir, paste0(ids, ".rds"))
    if (length(ids) && !all(file.copy(source, destination))) {
      stop("Could not stage every saved model; checkpoints were preserved.", call. = FALSE)
    }
  }
  jsonlite::write_json(
    .runRecord(tasks, listed, started, save_fits, BPPARAM),
    file.path(staging, "run.json"), pretty = TRUE, auto_unbox = TRUE
  )
  files <- list.files(staging, recursive = TRUE, full.names = FALSE)
  checksums <- tools::md5sum(file.path(staging, files))
  names(checksums) <- files
  .saveFitCheckpoint(list(manifest = manifest, checksums = checksums),
                     file.path(staging, "results", "fit-completion.rds"))
  .restoreFitResults(staging, out_dir, produced)
  unlink(checkpoint_dir, recursive = TRUE)
  fit
}
