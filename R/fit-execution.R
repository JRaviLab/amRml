# Validate and write resumable fits

# An explicit worker receives one group and reads its feature stores locally.
.runTaskGroup <- function(job, species, model_dir, save_fits, cache) {
  shared_fit <- NULL
  outcomes <- vector("list", nrow(job$tasks))
  for (j in seq_len(nrow(job$tasks))) {
    task <- job$tasks[j, ]
    keep_fit <- nrow(job$tasks) > 1L && is.null(shared_fit)
    task_notes <- .emptyModelNotes()
    outcome <- tryCatch(
      withCallingHandlers(withr::with_seed(
        task$seed,
        .fitTask(
          task = task, path = job$paths[[j]],
          species = species,
          rows = job$rows[[j]], columns = job$columns[[j]],
          encoding = job$encoding[[j]], cache = cache,
          fitted = shared_fit, return_fitted = keep_fit, multi = job$multi[[j]]
        ),
        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
        .rng_sample_kind = "Rejection"
      ), warning = function(w) {
        task_notes <<- dplyr::bind_rows(task_notes, .conditionModelNote(w, "task"))
        invokeRestart("muffleWarning")
      }, message = function(m) {
        task_notes <<- dplyr::bind_rows(task_notes, .conditionModelNote(m, "task"))
        invokeRestart("muffleMessage")
      }),
      error = function(e) {
        task_notes <<- dplyr::bind_rows(task_notes, e$amr_notes)
        .fitDiagnostic(job$diagnostics, if (.isMemoryError(conditionMessage(e))) "oom_error" else "task_error",
                       list(task_id = task$task_id, group = basename(job$checkpoint),
                            error = conditionMessage(e)))
        list(error = conditionMessage(e))
      }
    )
    outcome$notes <- dplyr::bind_rows(outcome$notes, dplyr::bind_cols(
      tibble::tibble(task_id = rep(task$task_id, nrow(task_notes)),
                     species = rep(species, nrow(task_notes))), task_notes
    ))
    if (is.null(outcome[["error"]])) {
      if (keep_fit) shared_fit <- outcome$fitted
      outcome$fitted <- NULL
      if (save_fits) {
        .saveFitCheckpoint(
          outcome$model, file.path(model_dir, paste0(task$task_id, ".rds"))
        )
      }
      outcome$model <- NULL
    }
    outcomes[[j]] <- outcome
  }
  # Commit the group only after any requested model files have been saved.
  .saveFitCheckpoint(outcomes, job$checkpoint)
  TRUE
}

# How long is this thing taking? May need to revise format for long runs
.elapsedTime <- function(started) {
  seconds <- floor(as.numeric(difftime(Sys.time(), started, units = "secs")))
  sprintf("%02d:%02d", seconds %/% 60, seconds %% 60)
}

# Workers return completion counts; only the parent writes progress or manifests.
.fitProgressReporter <- function(started, manifest_path, total) {
  force(started)
  force(manifest_path)
  force(total)
  last_time <- NULL
  last_completed <- NULL
  function(completed, finish = FALSE) {
    now <- Sys.time()
    if (identical(completed, last_completed) && !finish) return(invisible(NULL))
    if (!finish && !is.null(last_time) &&
        as.numeric(difftime(now, last_time, units = "secs")) < 1) return(invisible(NULL))
    text <- sprintf("%d of %d models complete. %s elapsed.",
                    completed, total, .elapsedTime(started))
    if (!identical(completed, last_completed)) {
      if (interactive()) {
        cat("\r", text, "   ", sep = "")
        utils::flush.console()
      } else message(text)
      if (file.exists(manifest_path)) .manifestEvent(manifest_path, text, details = list(
        n_completed = completed, n_tasks = total,
        elapsed_seconds = as.numeric(difftime(now, started, units = "secs"))
      ))
    }
    last_time <<- now
    last_completed <<- completed
    if (finish && interactive()) cat("\n")
    invisible(NULL)
  }
}

.runFitGroup <- function(job, ...) {
  if (!is.null(job$lock)) {
    .acquireFitLock(job$lock, "Task group")
    on.exit(unlink(job$lock, recursive = TRUE), add = TRUE)
  }
  details <- list(group = basename(job$checkpoint), task_ids = job$tasks$task_id,
                  matrix_ids = job$tasks$matrix_id, matrix_files = unique(job$paths),
                  r_heap_used_mb = sum(gc(full = FALSE)[, 2]))
  .fitDiagnostic(job$diagnostics, "group_start", details)
  message("Running ", details$group, ": ", paste(details$matrix_ids, collapse = ", "),
          " (R heap ", details$r_heap_used_mb, " MiB).")
  tryCatch(.runTaskGroup(job, ...), error = function(e) {
    .fitDiagnostic(job$diagnostics, if (.isMemoryError(conditionMessage(e))) "oom_error" else "group_error",
                   list(group = details$group, error = conditionMessage(e)))
    stop(e)
  })
  .fitDiagnostic(job$diagnostics, "group_complete", list(group = details$group))
  nrow(job$tasks)
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
    multi = .isMulticlass(matrices$unit[positions], matrices$target[positions]),
    checkpoint = checkpoint,
    diagnostics = file.path(out_dir, "diagnostics")
  )
}

# Describe the tasks and matrix files a checkpoint belongs to. Checksums prevent
# completed fits from being mixed with changed training inputs when resuming.
.fitCheckpointManifest <- function(tasks, groups, reuse_fits, save_fits, save_state) {
  files <- tasks$built$matrices$file[
    match(tasks$tasks$matrix_id, tasks$built$matrices$matrix_id)
  ]
  inputs <- unique(c(file.path("state", "matrix_index.rds"), files))
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
    save_state = save_state,
    species = tasks$built$eligibility$profile$orb$dataset_id,
    input_checksums = checksums
  )
}

# Write one completed task group before its worker finishes. A temporary file
# prevents an interrupted write from appearing as a finished checkpoint.
.saveFitCheckpoint <- function(result, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
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
    failure <- is.list(x) && "error" %in% names(x) &&
      all(names(x) %in% c("error", "notes")) &&
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
    notes = dplyr::bind_cols(
      tibble::tibble(task_id = character(), species = character()), .emptyModelNotes()
    ),
    performance = tibble::tibble(
      task_id = character(), species = character(), num_obs = integer(), res_prop = numeric(),
      n_test = integer(), res_prop_test = numeric(), n_feat = integer(),
      n_labels = integer(), model = character(), n_feats_returned = integer(), n_fold = integer(),
      fit_penalty = numeric(), fit_mixture = numeric(), mcc = numeric(), nmcc = numeric(),
      spec = numeric(), sens = numeric(), log2_apop = numeric(), f1 = numeric(),
      bal_acc = numeric(), auprc = numeric(), run_time_sec = numeric(),
      seed = integer(), date = character(), metric_undefined = logical(),
      collapsed_model = logical(), collapse_class = character(),
      prediction_behavior = character(), tp = integer(), tn = integer(),
      fp = integer(), fn = integer(), n_pred_resistant = integer(),
      n_pred_susceptible = integer(), pred_resistant_fraction = numeric(),
      actual_resistant_fraction = numeric()
    ),
    importance = tibble::tibble(
      task_id = character(), species = character(), label = character(), Variable = character(),
      Importance = numeric(), Sign = character()
    ),
    predictions = tibble::tibble(
      task_id = character(), species = character(), genome_id = character(), truth = character(),
      predicted = character(), prob_resistant = numeric(), prob_predicted = numeric()
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
        dir.create(dirname(file.path(backup, name)), recursive = TRUE, showWarnings = FALSE)
        if (!file.rename(file.path(out_dir, name), file.path(backup, name))) {
          stop("Could not move previous result: ", name)
        }
        old_moved <- c(old_moved, name)
      }
      dir.create(dirname(file.path(out_dir, name)), recursive = TRUE, showWarnings = FALSE)
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
  paths <- .runPaths(out_dir)
  path <- paths$completion
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
    arrow::read_parquet(paths[[name]])
  })
  names(results) <- c("performance", "importance", "predictions")
  results$notes <- if (file.exists(paths$notes)) arrow::read_parquet(paths$notes) else
    .emptyFitResults()$notes
  structure(c(list(tasks = arrow::read_parquet(paths$tasks)),
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
#' Each completed task group is saved in `state/fit-checkpoints/` so it can be reused
#' after an interruption. Checkpoints are removed after results are written.
#' Use `resume = TRUE` with the same tasks and output directory to reuse them.
#' Each task runs under its own seed. A task that errors is recorded as
#' `fit_failed` or `fit_oom` for memory pressure errors, and one with undefined
#' metrics as `metrics_undefined`; the rest carry on. A shuffled task permutes
#' the training matrix's labels before the split.
#' Progress is reported by the parent process as checkpointed task groups return.
#' Counts include failed tasks and recovered checkpoints; the final summary lists
#' fitted and failed tasks separately. Elapsed time covers the current invocation.
#'
#' Completed results can also be reopened with `resume = TRUE` without fitting
#' again when `save_state = TRUE`.
#'
#' @param tasks An `amr_tasks` from [expandTasks()].
#' @param save_fits Save each fitted model to `results/models/<task_id>.rds`:
#'   the fitted workflow (`fit`) and the `features` its predictors `x1`, `x2`,
#'   ... stand for. Default `FALSE`.
#'   Must match the interrupted run when resuming.
#' @param save_state Retain the internal RDS matrix stores and resume state after
#'   successful fitting. Default `FALSE`; interrupted runs retain this state so
#'   they can still be resumed. Must match when resuming.
#' @param overwrite Replace results already in the output folder.
#' @param reuse_fits Reuse fitted models among cross-drug, cross-year, and
#'   cross-country tasks with identical training inputs, seeds, and fold settings.
#' @param resume Reuse saved task groups from an interrupted fit or reopen an
#'   unchanged completed fit. The task plan, matrix files and `reuse_fits`
#'   setting must match the original run.
#'   `save_fits` and `save_state` must also match.
#' @param execution Execution phase: `"fit"` (default), `"prepare"` to save an
#'   array plan, `"worker"` to checkpoint selected groups only, or `"collect"`
#'   to publish results after every group has finished. Workers and collection
#'   require `resume = TRUE`. See [runTaskChunk()] and [collectTaskArray()].
#' @param group_subset Group indices; required only for `execution = "worker"`.
#'   Use the saved chunk map rather than recomputing group numbers.
#' @param n_chunks Maximum number of array chunks when preparing.
#'   Fewer chunks assign more fit groups to each worker; see [runTaskChunk()].
#' @param BPPARAM A `BiocParallelParam`, e.g. `BiocParallel::MulticoreParam(4)`,
#'   or `BiocParallel::SerialParam()` to fit one task at a time. The default is
#'   `BiocParallel::bpparam()`, except that inside Positron or RStudio, which
#'   cannot fork the R session, it is a `SnowParam` with the same number of
#'   workers. Workers are fresh R processes, so amRml must be installed.
#' @return Preparation invisibly returns a plan (`tasks`, `chunks`, and fit
#'   settings); workers return a list of group IDs. Otherwise an `amr_fit` list,
#'   also written to the output folder: `tasks` (with
#'   `status`, `rule_id` and `message`; `results/tasks.parquet`), `performance`,
#'   `importance`, `predictions`, and `notes` (`results/`), and the `built` matrices.
#'   Notes retain tuning diagnostics and task warnings/messages, with task IDs,
#'   stage, resample, location, type, text, and occurrence counts. In
#'   `importance`, `Sign` "POS" points to Susceptible and "NEG" to Resistant. The
#'   Internal state is retained after interruption, or after success when
#'   `save_state = TRUE`. When called through [runModelingWorkflow()], settings,
#'   provenance and task counts are also recorded in `manifest.json`.
#'
#' MDR combination models retain per-label coefficients in `importance$label`.
#' Their binary-only performance metrics and `prob_resistant` are `NA`;
#' `n_labels` records label counts and `prob_predicted` the predicted-label probability.
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
fitModels <- function(tasks, save_fits = FALSE, save_state = FALSE,
                      BPPARAM = .defaultBPPARAM(), overwrite = FALSE,
                      reuse_fits = TRUE, resume = FALSE,
                      group_subset = NULL,
                      execution = c("fit", "prepare", "worker", "collect"),
                      n_chunks = 100L) {
  call <- rlang::current_env()
  started <- Sys.time()
  .checkArgClass(tasks, "tasks", "amr_tasks", "expandTasks()", call = call)
  .checkArgClass(BPPARAM, "BPPARAM", "BiocParallelParam", "BiocParallel", call = call)
  .checkArgFlag(overwrite, "overwrite", call = call)
  .checkArgFlag(reuse_fits, "reuse_fits", call = call)
  .checkArgFlag(resume, "resume", call = call)
  .checkArgFlag(save_fits, "save_fits", call = call)
  .checkArgFlag(save_state, "save_state", call = call)
  .checkArgGroupSubset(group_subset, call)
  execution <- match.arg(execution)
  .checkArgCount(n_chunks, "n_chunks", call = call)
  if ((execution == "worker") != !is.null(group_subset)) {
    .amrAbort("invalid_argument",
              "`group_subset` is required only for `execution = \"worker\"`.", call = call)
  }
  if (execution %in% c("worker", "collect") && !resume) {
    .amrAbort("invalid_argument", "Workers and collection require `resume = TRUE`.", call = call)
  }
  if (execution == "prepare" && resume) {
    .amrAbort("invalid_argument", "Preparation requires `resume = FALSE`.", call = call)
  }
  if (resume && overwrite) {
    .amrAbort("invalid_argument", "`resume` and `overwrite` cannot both be TRUE.", call = call)
  }
  out_dir <- tasks$built$out_dir
  paths <- .runPaths(out_dir)
  produced <- c("results", file.path("state", "fit-completion.rds"))
  if (any(file.exists(file.path(out_dir, produced))) && !overwrite && !resume) {
    .amrAbort("invalid_argument",
              "The output folder already has results; set `overwrite = TRUE` to replace them.",
              observed = list(out_dir = out_dir), call = call)
  }

  species <- tasks$built$eligibility$profile$orb$dataset_id
  listed <- tasks$tasks
  listed$status <- rep("fitted", nrow(listed))
  listed$rule_id <- rep(NA_character_, nrow(listed))
  listed$message <- rep(NA_character_, nrow(listed))
  if (execution == "worker") message("Array worker: loading matrix index and checking task groups.")
  index <- readRDS(paths$matrix_index)
  groups <- if (reuse_fits) {
    .fitTaskGroups(listed, tasks$built$matrices, index)
  } else {
    as.list(seq_len(nrow(listed)))
  }
  checkpoint_dir <- paths$checkpoints
  manifest_path <- file.path(checkpoint_dir, "manifest.rds")
  if (resume) {
    saved_contract <- if (file.exists(manifest_path)) readRDS(manifest_path) else
      if (file.exists(paths$completion)) readRDS(paths$completion)$manifest else NULL
    saved_groups <- saved_contract$groups
    membership <- integer(nrow(listed))
    for (g in seq_along(groups)) membership[groups[[g]]] <- g
    flat <- unlist(saved_groups, use.names = FALSE)
    # Older runs may contain fewer shared fits. Keep their checkpoint numbering
    # when those groups remain compatible; all other inputs are validated below.
    if (is.list(saved_groups) && is.integer(flat) &&
        identical(sort(flat), seq_len(nrow(listed))) &&
        all(vapply(saved_groups, function(g) length(g) > 0L &&
          length(unique(membership[g])) == 1L, logical(1)))) {
      groups <- saved_groups
    }
  }
  if (!is.null(group_subset) && any(group_subset > length(groups))) {
    .amrAbort("invalid_argument", "`group_subset` contains invalid group indices.", call = call)
  }
  array_path <- file.path(checkpoint_dir, "array.rds")
  if (execution == "fit" && file.exists(array_path)) {
    .amrAbort("invalid_argument", "This is an array run; please use `collectTaskArray()` to finish!", call = call)
  }
  if (execution == "worker") message("Array worker: verifying shared input checksums.")
  manifest <- .fitCheckpointManifest(tasks, groups, reuse_fits, save_fits, save_state)
  if (resume) {
    if (!file.exists(manifest_path)) return(.readCompletedFit(tasks, manifest, call))
    if (!identical(readRDS(manifest_path), manifest)) {
      .amrAbort("invalid_argument",
                "Saved fits don't match these tasks, settings or matrices; can't resume.",
                call = call)
    }
  } else {
    if (dir.exists(checkpoint_dir) && !overwrite) {
      .amrAbort("invalid_argument",
                "Task checkpoints already exist! Use `resume = TRUE` to reuse.", call = call)
    }
    unlink(checkpoint_dir, recursive = TRUE)
    dir.create(checkpoint_dir, recursive = TRUE)
    .saveFitCheckpoint(manifest, manifest_path)
  }
  model_dir <- file.path(checkpoint_dir, "models")
  if (save_fits) dir.create(model_dir, showWarnings = FALSE)
  group_paths <- file.path(checkpoint_dir, sprintf("group-%06d.rds", seq_along(groups)))
  if (execution == "prepare") {
    plan <- list(tasks = tasks, chunks = .splitGroupIds(length(groups), n_chunks),
                 save_fits = save_fits, save_state = save_state, reuse_fits = reuse_fits)
    .saveFitCheckpoint(plan, array_path)
    return(invisible(plan))
  }
  selected <- if (execution == "worker") sort(unique(as.integer(group_subset))) else seq_along(groups)
  if (execution == "collect") {
    if (length(list.files(checkpoint_dir, pattern = "^(chunk|group)-.*[.]lock$"))) {
      stop("Array workers may be still active, or stale locks remain; ",
           "inspect lock owner.json files and confirm jobs have exited before recovery.", call. = FALSE)
    }
    missing <- which(!file.exists(group_paths))
    if (length(missing)) stop("Cannot collect: missing task groups: ",
                              paste(missing, collapse = ", "), call. = FALSE)
  }
  recovered <- selected[file.exists(group_paths[selected])]
  # Validate recovered groups before scheduling more work.
  for (g in recovered) {
    .readFitCheckpoint(group_paths[[g]], listed$task_id[groups[[g]]], model_dir, save_fits)
  }
  pending <- selected[!file.exists(group_paths[selected])]
  if (resume) message(length(recovered), " task groups recovered; ",
                      length(pending), " task groups remaining.")
  jobs <- lapply(pending, function(g) {
    job <- .fitGroupJob(groups[[g]], listed, tasks$built$matrices, index, out_dir, group_paths[[g]])
    if (execution == "worker") job$lock <- paste0(group_paths[[g]], ".lock")
    job
  })
  # An empty cache is sent once; each worker populates its own copy across groups.
  cache <- new.env(parent = emptyenv())
  n_recovered <- sum(lengths(groups[recovered]))
  report <- .fitProgressReporter(
    started, if (execution == "worker") "" else paths$manifest,
    sum(lengths(groups[selected]))
  )
  report(n_recovered)
  if (length(jobs)) {
    next_job <- local({
      i <- 0L
      function() {
        i <<- i + 1L
        if (i <= length(jobs)) jobs[[i]] else NULL
      }
    })
    completed <- BiocParallel::bpiterate(
      next_job, .runFitGroup, cache = cache, species = species, model_dir = model_dir,
      save_fits = save_fits, BPPARAM = BPPARAM,
      REDUCE = function(completed, n) {
        report(completed + n)
        completed + n
      }, init = n_recovered, reduce.in.order = FALSE
    )
    report(completed, finish = TRUE)
  } else report(n_recovered, finish = TRUE)

  if (execution == "worker") return(invisible(list(groups = selected)))

  outcomes <- vector("list", nrow(listed))
  for (g in seq_along(groups)) {
    outcomes[groups[[g]]] <- .readFitCheckpoint(
      group_paths[[g]], listed$task_id[groups[[g]]], model_dir, save_fits
    )
  }
  task_notes <- lapply(outcomes, `[[`, "notes")
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[[i]]
    if (!is.null(outcome[["error"]])) {
      listed$status[[i]] <- "failed"
      listed$rule_id[[i]] <- if (.isMemoryError(outcome$error)) "fit_oom" else "fit_failed"
      listed$message[[i]] <- outcome$error
      outcomes[i] <- list(NULL)
    } else if (!is.na(outcome$undefined)) {
      listed$rule_id[[i]] <- "metrics_undefined"
      listed$message[[i]] <- paste("Undefined:", outcome$undefined)
    }
  }
  empty <- .emptyFitResults()
  tables <- lapply(names(empty), function(name) {
    dplyr::bind_rows(c(list(empty[[name]]),
                       if (name == "notes") task_notes else lapply(outcomes, `[[`, name)))
  })
  names(tables) <- names(empty)
  fit <- structure(c(list(tasks = listed), tables, list(built = tasks$built)), class = "amr_fit")

  staging <- tempfile(".results-building-", tmpdir = out_dir)
  stage_paths <- .runPaths(staging)
  dir.create(stage_paths$results, recursive = TRUE)
  dir.create(stage_paths$state, recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  .writeParquet(fit$tasks, stage_paths$tasks)
  for (name in names(empty)) {
    .writeParquet(fit[[name]], stage_paths[[name]])
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
  files <- list.files(staging, recursive = TRUE, full.names = FALSE)
  checksums <- tools::md5sum(file.path(staging, files))
  names(checksums) <- files
  .saveFitCheckpoint(list(manifest = manifest, checksums = checksums),
                     stage_paths$completion)
  .restoreFitResults(staging, out_dir, produced)
  unlink(checkpoint_dir, recursive = TRUE)
  if (!save_state) unlink(paths$state, recursive = TRUE)

  n_collapsed <- sum(fit$performance$collapsed_model, na.rm = TRUE)
  
  n_all_susceptible <- sum(fit$performance$prediction_behavior ==
                             "all_susceptible",
                           na.rm = TRUE)
  
  n_all_resistant <- sum(fit$performance$prediction_behavior ==
                           "all_resistant",
                         na.rm = TRUE)
  
  message(
    sprintf(
      "Model fitting complete: %d of %d fitted; %d failed. %s elapsed.",
      sum(listed$status == "fitted"),
      nrow(listed),
      sum(listed$status == "failed"),
      .elapsedTime(started)
    ),
    "\nDiagnostics:",
    "\n  Collapsed models: ", n_collapsed,
    "\n    All susceptible: ", n_all_susceptible,
    "\n    All resistant: ", n_all_resistant
  )
  message("  Diagnostic notes: ", nrow(fit$notes), " records in results/notes.parquet")
  fit
}
