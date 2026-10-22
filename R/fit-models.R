# Fitting a model for each task, with old amRml's method: glmnet logistic regression,
# tuned by stratified cross-validation and chosen by MCC.

# The phenotypes, with Resistant first so it is the event yardstick scores.
.PHENOTYPES <- c("Resistant", "Susceptible")
# Setting an initial cutoff for when matrices are treated as sparse in tidymodels
# but we should really benchmark this to see where the right threshold is
.SPARSE_THRESHOLD <- 0.8


# Read one scope from the shared sparse feature store and its index
.readMatrixWide <- function(path, rows, columns, encoding, cache) {
  if (!exists(path, envir = cache, inherits = FALSE)) {
    assign(path, readRDS(path), envir = cache)
  }
  store <- get(path, envir = cache, inherits = FALSE)
  genome_idx <- match(rows$genome_id, store$genomes)
  if (anyNA(genome_idx)) stop("Scope genome is missing from feature store.")
  values <- store$values[genome_idx, columns, drop = FALSE]
  if (encoding == "binary") {
    values@x[] <- as.numeric(values@x > 0)
    values <- Matrix::drop0(values)
  }

  list(
    rows = rows,
    features = store$features[columns],
    values = values
  )
}

# Ported from old amRml's shuffleLabels(), with the task's seed: the training matrix's labels
# are permuted before the split. Test rows (cross, leave-one-out) keep their real labels.
.shuffleLabels <- function(rows) {
  train <- which(rows$role == "train")
  rows$phenotype[train] <- rows$phenotype[train][sample.int(length(train))]
  rows
}

# Training and test rows. A matrix with test rows uses them; otherwise each phenotype keeps
# floor((1 - holdout) x n) genomes for training, as rsample's stratified split does.
# Training rows come grouped by phenotype, in rsample's order: with duplicate features,
# glmnet's coefficients depend on row order.
.splitTask <- function(rows, holdout) {
  if (any(rows$role == "test")) {
    return(list(train = which(rows$role == "train"), test = which(rows$role == "test")))
  }

  train <- unlist(lapply(split(seq_len(nrow(rows)), rows$phenotype), function(i) {
    sort(i[sample.int(length(i), floor((1 - holdout) * length(i)))])
  }), use.names = FALSE)
  list(train = train, test = setdiff(seq_len(nrow(rows)), train))
}

# Each training row's fold, assigned per phenotype as rsample's stratified folds do, but
# without pooling rare phenotypes (`pool = 0.1`), which can leave a fold without them.
.stratifiedFolds <- function(phenotype, n_fold) {
  fold <- integer(length(phenotype))
  for (level in intersect(.PHENOTYPES, phenotype)) {
    i <- which(phenotype == level)
    fold[i] <- sample(rep_len(seq_len(n_fold), length(i)))
  }
  fold
}

# The training data as a data frame: the label, then features renamed x1, x2, ..., so any
# feature ID works in a formula.
.modelData <- function(values, phenotype, sparse = FALSE) {
  if (sparse) {
    colnames(values) <- paste0("x", seq_len(ncol(values)))
    data <- sparsevctrs::coerce_to_sparse_tibble(values)
  } else {
    data <- as.data.frame(as.matrix(values))
  }
  names(data) <- paste0("x", seq_len(ncol(data)))
  data$phenotype <- factor(phenotype, levels = .PHENOTYPES)
  data
}

# Ported from old amRml's buildRecipe(): drop constant features, then normalize.
.buildRecipe <- function(data, sparse = FALSE) {
  recipe <- recipes::recipe(phenotype ~ ., data = data) |>
    recipes::step_zv(recipes::all_predictors())
  if (sparse) {
    # step_zv() alone seemingly misses constant nonzero sparse columns, so we're
    # filtering these out at prep to resolve those errors
    recipe <- recipes::step_rm(
      recipe,
      tidyselect::where(function(x) {
        sparsevctrs::is_sparse_vector(x) &&
          length(sparsevctrs::sparse_values(x)) == length(x) &&
          length(unique(sparsevctrs::sparse_values(x))) == 1L
      }),
      id = "zv_sparse_nonzero"
    )
    recipes::step_scale(recipe, recipes::all_predictors())
  } else {
    recipes::step_normalize(recipe, recipes::all_predictors())
  }
}

# Ported from old amRml's buildLRModel().
.modelSpec <- function() {
  parsnip::logistic_reg(penalty = tune::tune(), mixture = tune::tune()) |>
    parsnip::set_engine("glmnet")
}

# Ported from old amRml's buildTuningGrid(): 10 penalties x 6 mixtures.
.tuningGrid <- function() {
  penalty_vec <- 10^seq(-4, -1, length.out = 10)
  mix_vec <- 0:5 / 5
  tibble::tibble(
    penalty = rep(penalty_vec, each = length(mix_vec)),
    mixture = rep(mix_vec, length(penalty_vec))
  )
}

# Ported from old amRml's tuneGrid(), selectBestModel() and fitBestModel(): tune over
# stratified folds, choose by MCC, and fit on all training rows.
.fitModel <- function(data, n_fold, sparse = FALSE) {
  workflow <- workflows::workflow() |>
    workflows::add_recipe(.buildRecipe(data, sparse = sparse)) |>
    workflows::add_model(.modelSpec())

  # After the recipe, whose step IDs use random numbers, so the folds match rsample's.
  fold <- .stratifiedFolds(data$phenotype, n_fold)

  splits <- lapply(sort(unique(fold)), function(k) {
    rsample::make_splits(
      list(analysis = which(fold != k), assessment = which(fold == k)),
      data = data
    )
  })
  resamples <- rsample::manual_rset(splits, ids = paste0("Fold", sort(unique(fold))))

  tuned <- tune::tune_grid(
    workflow,
    resamples = resamples,
    grid = .tuningGrid(),
    metrics = yardstick::metric_set(
      yardstick::f_meas, yardstick::pr_auc, yardstick::spec, yardstick::sens,
      yardstick::bal_accuracy, yardstick::mcc
    ),
    control = tune::control_grid(allow_par = FALSE)
  )
  best <- tune::select_best(tuned, metric = "mcc")

  list(
    fit = parsnip::fit(tune::finalize_workflow(workflow, best), data = data),
    penalty = best$penalty,
    mixture = best$mixture
  )
}

# Ported from old amRml's .calculate*() helpers, keeping AUPRC. Resistant is the event.
# Each metric is rounded to 2, and nMCC and log2(AUPRC / prior) use the rounded values.
.taskMetrics <- function(scored) {
  estimate <- function(metric, ...) {
    round(as.numeric(metric(scored, truth = "phenotype", ...)$.estimate), 2)
  }
  mcc <- estimate(yardstick::mcc, estimate = ".pred_class")
  auprc <- estimate(yardstick::pr_auc, ".pred_Resistant")

  tibble::tibble(
    mcc = mcc,
    nmcc = round((mcc + 1) / 2, 2),
    spec = estimate(yardstick::spec, estimate = ".pred_class"),
    sens = estimate(yardstick::sens, estimate = ".pred_class"),
    log2_apop = round(log2(auprc / mean(scored$phenotype == "Resistant")), 2),
    f1 = estimate(yardstick::f_meas, estimate = ".pred_class"),
    bal_acc = estimate(yardstick::bal_accuracy, estimate = ".pred_class"),
    auprc = auprc
  )
}

# Ported from old amRml's .viGlmnet(): the non-zero coefficients at the chosen penalty, on
# normalized features. glmnet models the second level, so "POS" points to Susceptible.
.taskImportance <- function(fit, penalty, features) {
  coefs <- stats::coef(parsnip::extract_fit_engine(fit), s = penalty)[, 1]
  coefs <- coefs[names(coefs) != "(Intercept)" & coefs != 0]
  coefs <- coefs[order(-abs(coefs))]

  tibble::tibble(
    Variable = features[as.integer(sub("^x", "", names(coefs)))],
    Importance = abs(unname(coefs)),
    Sign = as.character(ifelse(unname(coefs) > 0, "POS", "NEG"))
  )
}

# Cross-test tasks can share a fit when the training inputs exactly match
.fitTaskGroups <- function(listed, matrices, index) {
  positions <- match(listed$matrix_id, matrices$matrix_id)
  signatures <- vector("list", nrow(listed))
  signature_set <- rep(NA_character_, nrow(listed))

  for (i in which(matrices$mode_id[positions] == "cross_drug")) {
    p <- positions[[i]]
    rows <- index$rows[[index$matrix_scopes[[p]]]]
    if (!any(rows$role == "test")) next

    train <- rows$role == "train"
    signatures[[i]] <- list(
      file = matrices$file[[p]],
      columns = index$columns[[p]],
      encoding = matrices$encoding[[p]],
      genome_id = rows$genome_id[train],
      phenotype = rows$phenotype[train],
      labels = listed$labels[[i]],
      seed = listed$seed[[i]],
      n_fold = listed$n_fold[[i]],
      holdout = listed$holdout[[i]]
    )
    # Collect all signature details
    signature_set[[i]] <- paste(
      signatures[[i]]$file, signatures[[i]]$encoding,
      signatures[[i]]$labels, signatures[[i]]$seed,
      signatures[[i]]$n_fold, signatures[[i]]$holdout,
      length(signatures[[i]]$columns), length(signatures[[i]]$genome_id),
      sep = "|"
    )
  }

  leader <- seq_len(nrow(listed))
  candidates <- which(!is.na(signature_set))
  for (part in split(candidates, signature_set[candidates])) {
    representatives <- integer()
    for (i in part) {
      found <- representatives[vapply(representatives, function(j) {
        identical(signatures[[i]], signatures[[j]])
      }, logical(1))]
      if (length(found)) {
        leader[[i]] <- found[[1L]]
      } else {
        representatives <- c(representatives, i)
      }
    }
  }

  split(seq_len(nrow(listed)), factor(leader, levels = unique(leader)))
}

# One task: shuffle if asked, split, tune, fit, and score on the test genomes.
.fitTask <- function(task, path, rows, columns, encoding, cache,
                     fitted = NULL, return_fitted = FALSE) {
  started <- Sys.time()
  matrix <- .readMatrixWide(path, rows, columns, encoding, cache)
  rows <- matrix$rows
  if (task$labels == "shuffled") {
    rows <- .shuffleLabels(rows)
  }

  split <- .splitTask(rows, task$holdout)
  train_values <- matrix$values[split$train, , drop = FALSE]
  # Not sure if there's a function to calculate how sparsely populated a matrix
  # is, but here finding it by non-zero values over total values in the matrix
  # and taking the inverse by subtracting the populated fraction from 1
  sparsity <- 1 - Matrix::nnzero(train_values) / (as.double(nrow(train_values)) * ncol(train_values))
  # Is that fraction greater than or equal to our cutoff? If so, that's sparse
  use_sparse <- sparsity >= .SPARSE_THRESHOLD
  train <- .modelData(train_values, rows$phenotype[split$train], sparse = use_sparse)
  test <- .modelData(matrix$values[split$test, , drop = FALSE],
                     rows$phenotype[split$test], sparse = use_sparse)
  if (is.null(fitted)) fitted <- .fitModel(train, task$n_fold, sparse = use_sparse)

  scored <- parsnip::augment(fitted$fit, test)
  metrics <- .taskMetrics(scored)
  importance <- .taskImportance(fitted$fit, fitted$penalty, matrix$features)
  input <- rows$role == "train"

  performance <- tibble::tibble(
    task_id = task$task_id,
    num_obs = sum(input),
    res_prop = round(mean(rows$phenotype[input] == "Resistant"), 2),
    n_test = length(split$test),
    res_prop_test = round(mean(test$phenotype == "Resistant"), 2),
    n_feat = length(matrix$features),
    model = "LR",
    n_feats_returned = nrow(importance),
    n_fold = task$n_fold,
    fit_penalty = fitted$penalty,
    fit_mixture = fitted$mixture,
    metrics,
    run_time_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2),
    seed = task$seed,
    date = as.character(Sys.Date())
  )
  predictions <- tibble::tibble(
    task_id = task$task_id,
    genome_id = rows$genome_id[split$test],
    truth = as.character(scored$phenotype),
    predicted = as.character(scored$.pred_class),
    prob_resistant = scored$.pred_Resistant
  )

  # Metrics that can't be computed, such as MCC when every test genome gets one prediction.
  undefined <- names(metrics)[vapply(metrics, is.na, logical(1))]
  if (length(undefined) && dplyr::n_distinct(scored$.pred_class) == 1) {
    undefined <- paste0(
      toString(undefined), "; every test genome was predicted ", scored$.pred_class[[1]]
    )
  }

  result <- list(
    performance = performance,
    importance = tibble::tibble(task_id = rep(task$task_id, nrow(importance)), importance),
    predictions = predictions,
    undefined = if (length(undefined)) toString(undefined) else NA_character_
  )
  if (return_fitted) result$fitted <- fitted
  result
}

# Describe the tasks and matrix files a checkpoint belongs to. Checksums prevent
# completed fits from being mixed with changed training inputs when resuming.
## This is draft functionality, and should be replaced with proper ORB manifest
## logging after verification, but modeling resume functionality is badly needed.
.fitCheckpointManifest <- function(tasks, groups, reuse_fits) {
  out_dir <- tasks$built$out_dir
  listed <- tasks$tasks
  files <- tasks$built$matrices$file[
    match(listed$matrix_id, tasks$built$matrices$matrix_id)
  ]
  inputs <- unique(c("matrix_index.rds", files))
  checksums <- tools::md5sum(file.path(out_dir, inputs))
  names(checksums) <- inputs
  if (anyNA(checksums)) {
    stop("Cannot read the matrix files needed for checkpoints.", call. = FALSE)
  }

  list(
    version = 1L,
    tasks = listed,
    matrices = tasks$built$matrices,
    groups = groups,
    reuse_fits = reuse_fits,
    input_checksums = checksums
  )
}

# Write one completed task group before its worker finishes. A temporary file
# prevents an interrupted write from appearing as a finished checkpoint.
.saveFitCheckpoint <- function(result, path) {
  temporary <- tempfile(".group-", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(result, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) {
    stop("Could not save fitted task group: ", path, call. = FALSE)
  }
  invisible(path)
}

# Fail invalid checkpoints rather than treating them as model fits
.readFitCheckpoint <- function(path, n_tasks) {
  result <- readRDS(path)
  valid <- is.list(result) && length(result) == n_tasks &&
    all(vapply(result, function(x) {
      inherits(x, "error") ||
        (is.list(x) && all(c("performance", "importance", "predictions", "undefined") %in% names(x)))
    }, logical(1)))
  if (!valid) {
    stop("Invalid fitted task checkpoint: ", path, call. = FALSE)
  }
  result
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
#' @param tasks An `amr_tasks` from [expandTasks()].
#' @param overwrite Replace results already in the output folder.
#' @param reuse_fits Reuse fitted models among compatible cross-drug tasks.
#' @param resume Reuse saved task groups from an interrupted fit. The task plan,
#'   matrix files and `reuse_fits` setting must match the original run.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] backend; by default,
#'   an eight-worker socket cluster with 128 scheduling chunks. Use
#'   [BiocParallel::SerialParam()] to run sequentially.
#'
#' @return An `amr_fit` list, also written to the output folder: `tasks` (with
#'   `status`, `rule_id` and `message`; `tasks.parquet`), and `performance`,
#'   `importance` and `predictions` (`results/`), plus the `built` matrices. `performance`
#'   keeps old amRml's columns and rounding, adding `auprc`;
#'   `n_feat` counts the matrix's features (the recipe may drop more).
#'   `importance` has old amRml's `Variable`, `Importance` and `Sign` ("POS"
#'   points to Susceptible).
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus `amrml_invalid_argument`: `tasks`
#' is not from [expandTasks()], `BPPARAM` is not a BiocParallel backend,
#' or results already exist and `overwrite` is `FALSE`.
#'
#' @examples
#' \dontrun{
#' fit <- fitModels(expandTasks(built))
#' fit$performance
#' }
#' @export
fitModels <- function(tasks, overwrite = FALSE,
                      BPPARAM = BiocParallel::SnowParam(workers = 8L, tasks = 128L),
                      reuse_fits = TRUE, resume = FALSE) {
  call <- rlang::current_env()
  .checkArgClass(tasks, "tasks", "amr_tasks", "expandTasks()", call = call)
  .checkArgClass(BPPARAM, "BPPARAM", "BiocParallelParam", "BiocParallel::SnowParam()", call = call)
  .checkArgFlag(reuse_fits, "reuse_fits", call = call)
  .checkArgFlag(resume, "resume", call = call)
  if (resume && overwrite) {
    .amrAbort("invalid_argument", "`resume` and `overwrite` cannot both be TRUE.", call = call)
  }
  out_dir <- tasks$built$out_dir
  results <- file.path(out_dir, c("tasks.parquet", "results"))
  .checkArgFlag(overwrite, "overwrite", call = call)
  if (any(file.exists(results)) && !overwrite && !resume) {
    .amrAbort(
      "invalid_argument",
      "The output folder already has results; set `overwrite = TRUE` to replace them.",
      observed = list(out_dir = out_dir),
      call = call
    )
  }

  listed <- tasks$tasks
  files <- tasks$built$matrices$file[match(listed$matrix_id, tasks$built$matrices$matrix_id)]
  listed$status <- "fitted"
  listed$rule_id <- NA_character_
  listed$message <- NA_character_
  fitted <- vector("list", nrow(listed))
  index <- readRDS(file.path(out_dir, "matrix_index.rds"))
  cache <- new.env(parent = emptyenv())
  matrix_positions <- match(listed$matrix_id, tasks$built$matrices$matrix_id)

  groups <- if (reuse_fits) {
    .fitTaskGroups(listed, tasks$built$matrices, index)
  } else {
    as.list(seq_len(nrow(listed)))
  }

  # Check existing fits against the tasks and matrices, only overwrite old checkpoints
  # on a new run when overwrite = TRUE
  checkpoint_dir <- file.path(out_dir, "fit-checkpoints")
  manifest_path <- file.path(checkpoint_dir, "manifest.rds")
  manifest <- .fitCheckpointManifest(tasks, groups, reuse_fits)
  if (resume) {
    if (!file.exists(manifest_path) ||
        !identical(readRDS(manifest_path), manifest)) {
      .amrAbort(
        "invalid_argument",
        "Saved fits do not match these tasks or matrix files; cannot resume.",
        call = call
      )
    }
  } else {
    if (dir.exists(checkpoint_dir) && !overwrite) {
      .amrAbort(
        "invalid_argument",
        "Task checkpoints already exist; use `resume = TRUE` to reuse them.",
        call = call
      )
    }
    unlink(checkpoint_dir, recursive = TRUE)
    dir.create(checkpoint_dir, recursive = TRUE)
    saveRDS(manifest, manifest_path)
  }

  group_paths <- file.path(checkpoint_dir, sprintf("group-%06d.rds", seq_along(groups)))
  pending <- which(!file.exists(group_paths))
  if (resume) {
    message(sum(file.exists(group_paths)), " task groups recovered; ",
            length(pending), " task groups remaining.")
  }

  # Workers save each finished group immediately, only fit missing groups on resumption.
  if (length(pending)) BiocParallel::bplapply(pending, function(g) {
    group <- groups[[g]]
    shared_fit <- NULL
    group_results <- vector("list", length(group))
    for (j in seq_along(group)) {
      i <- group[[j]]
      task <- listed[i, ]
      keep_fit <- length(group) > 1L && is.null(shared_fit)
      result <- tryCatch(
        withr::with_seed(
          task$seed,
          .fitTask(
            task,
            file.path(out_dir, files[[i]]),
            rows = index$rows[[index$matrix_scopes[[matrix_positions[[i]]]]]],
            columns = index$columns[[matrix_positions[[i]]]],
            encoding = tasks$built$matrices$encoding[[matrix_positions[[i]]]],
            cache = cache,
            fitted = shared_fit,
            return_fitted = keep_fit
          ),
          .rng_kind = "Mersenne-Twister",
          .rng_normal_kind = "Inversion",
          .rng_sample_kind = "Rejection"
        ),
        error = function(e)
          e
      )
      if (!inherits(result, "error") && keep_fit) {
        shared_fit <- result$fitted
        result$fitted <- NULL
      }
      group_results[[j]] <- result
    }
    .saveFitCheckpoint(group_results, group_paths[[g]])
    TRUE
  }, BPPARAM = BPPARAM)

  # Read every saved group, including groups completed by an earlier run
  results_by_group <- lapply(seq_along(groups), function(g) {
    .readFitCheckpoint(group_paths[[g]], length(groups[[g]]))
  })

  results_by_task <- vector("list", nrow(listed))
  for (g in seq_along(groups)) {
    results_by_task[groups[[g]]] <- results_by_group[[g]]
  }

  for (i in seq_len(nrow(listed))) {
    result <- results_by_task[[i]]
    if (inherits(result, "error")) {
      listed$status[[i]] <- "failed"
      listed$rule_id[[i]] <- "fit_failed"
      listed$message[[i]] <- conditionMessage(result)
      next
    }
    fitted[[i]] <- result
    if (!is.na(result$undefined)) {
      listed$rule_id[[i]] <- "metrics_undefined"
      listed$message[[i]] <- paste("Undefined:", result$undefined)
    }
  }

  fit <- list(
    tasks = listed,
    performance = dplyr::bind_rows(lapply(fitted, `[[`, "performance")),
    importance = dplyr::bind_rows(lapply(fitted, function(x) {
      if (is.null(x)) return(NULL)
      x$importance$Sign <- as.character(x$importance$Sign)
      x$importance
    })),
    predictions = dplyr::bind_rows(lapply(fitted, `[[`, "predictions")),
    built = tasks$built
  )

  # Formatting results table if no task fitted
  if (!any(listed$status == "fitted")) {
    fit$performance <- tibble::tibble(
      task_id = character(),
      num_obs = integer(),
      res_prop = numeric(),
      n_test = integer(),
      res_prop_test = numeric(),
      n_feat = integer(),
      model = character(),
      n_feats_returned = integer(),
      n_fold = integer(),
      fit_penalty = numeric(),
      fit_mixture = numeric(),
      mcc = numeric(),
      nmcc = numeric(),
      spec = numeric(),
      sens = numeric(),
      log2_apop = numeric(),
      f1 = numeric(),
      bal_acc = numeric(),
      auprc = numeric(),
      run_time_sec = numeric(),
      seed = integer(),
      date = character()
    )
    fit$importance <- tibble::tibble(
      task_id = character(),
      Variable = character(),
      Importance = numeric(),
      Sign = character()
    )
    fit$predictions <- tibble::tibble(
      task_id = character(),
      genome_id = character(),
      truth = character(),
      predicted = character(),
      prob_resistant = numeric()
    )
  }

  # Keep checkpoints until all Parquet files are finished
  staging <- tempfile(".results-building-", tmpdir = out_dir)
  dir.create(file.path(staging, "results"), recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  arrow::write_parquet(fit$tasks, file.path(staging, "tasks.parquet"))
  for (table in c("performance", "importance", "predictions")) {
    arrow::write_parquet(fit[[table]], file.path(staging, "results", paste0(table, ".parquet")))
  }

  backup <- tempfile(".results-previous-", tmpdir = out_dir)
  dir.create(backup)
  old_moved <- character()
  new_moved <- character()
  tryCatch({
    for (name in c("tasks.parquet", "results")) {
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
  unlink(checkpoint_dir, recursive = TRUE)

  structure(fit, class = "amr_fit")
}


#' @export
print.amr_fit <- function(x, ...) {
  orb <- x$built$eligibility$profile$orb
  t <- x$tasks
  modes <- .taskModes(t, x$built$matrices)
  by_mode <- split(t, factor(modes, unique(modes)))
  counts <- function(rows) {
    paste0(
      sum(rows$status == "fitted"), " fitted (", sum(rows$rule_id %in% "metrics_undefined"),
      " with undefined metrics), ", sum(rows$status == "failed"), " failed"
    )
  }

  cat("<amr_fit>", orb$dataset_id, "-", orb$dataset_label, "\n")
  cat("  in", x$built$out_dir, "\n")
  cat("  ", nrow(t), " tasks: ", counts(t), "\n", sep = "")
  cat(.modeLines(names(by_mode), vapply(by_mode, counts, character(1))), sep = "\n")
  invisible(x)
}
