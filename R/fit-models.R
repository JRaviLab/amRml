# Fitting a model for each task, with old amRml's method: glmnet logistic regression,
# tuned by stratified cross-validation and chosen by MCC.

# The phenotypes, with Resistant first so it is the event yardstick scores.
.PHENOTYPES <- c("Resistant", "Susceptible")


# A matrix file as its genomes (`rows`: genome_id, role, phenotype) and a genomes x
# features matrix. Files are complete grids sorted by genome then feature.
.readMatrixWide <- function(path) {
  long <- arrow::read_parquet(path)
  features <- long$feature_id[long$genome_id == long$genome_id[[1]]]
  rows <- long[seq(1, nrow(long), by = length(features)), c("genome_id", "role", "phenotype")]

  list(
    rows = rows,
    features = features,
    values = matrix(long$value, nrow = nrow(rows), byrow = TRUE)
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
.modelData <- function(values, phenotype) {
  data <- as.data.frame(values)
  names(data) <- paste0("x", seq_len(ncol(data)))
  data$phenotype <- factor(phenotype, levels = .PHENOTYPES)
  data
}

# Ported from old amRml's buildRecipe(): drop constant features, then normalize.
.buildRecipe <- function(data) {
  recipes::recipe(phenotype ~ ., data = data) |>
    recipes::step_zv(recipes::all_predictors()) |>
    recipes::step_normalize(recipes::all_predictors())
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
.fitModel <- function(data, n_fold) {
  workflow <- workflows::workflow() |>
    workflows::add_recipe(.buildRecipe(data)) |>
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
    )
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
    Sign = ifelse(unname(coefs) > 0, "POS", "NEG")
  )
}

# One task: shuffle if asked, split, tune, fit, and score on the test genomes.
.fitTask <- function(task, path) {
  started <- Sys.time()
  matrix <- .readMatrixWide(path)
  rows <- matrix$rows
  if (task$labels == "shuffled") {
    rows <- .shuffleLabels(rows)
  }

  split <- .splitTask(rows, task$holdout)
  train <- .modelData(matrix$values[split$train, , drop = FALSE], rows$phenotype[split$train])
  test <- .modelData(matrix$values[split$test, , drop = FALSE], rows$phenotype[split$test])
  fitted <- .fitModel(train, task$n_fold)

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

  list(
    performance = performance,
    importance = tibble::tibble(task_id = task$task_id, importance),
    predictions = predictions,
    undefined = if (length(undefined)) toString(undefined) else NA_character_
  )
}

#' Fit the listed models
#'
#' Fits every task from [expandTasks()], one at a time, with old amRml's method:
#' glmnet logistic regression over 10 penalties and 6 mixtures, tuned by
#' stratified `n_fold` cross-validation and chosen by MCC, then scored on the
#' task's test genomes.
#'
#' @details
#' Each task runs under its own seed. A task that errors is recorded as
#' `fit_failed`, and one with undefined metrics as `metrics_undefined`; the rest
#' carry on. A shuffled task permutes the training matrix's labels before the
#' split.
#'
#' @param tasks An `amr_tasks` from [expandTasks()].
#' @param overwrite Replace results already in the output folder.
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
#' is not from [expandTasks()], or results already exist and `overwrite` is
#' `FALSE`.
#'
#' @examples
#' \dontrun{
#' fit <- fitModels(expandTasks(built))
#' fit$performance
#' }
#' @export
fitModels <- function(tasks, overwrite = FALSE) {
  call <- rlang::current_env()
  .checkArgClass(tasks, "tasks", "amr_tasks", "expandTasks()", call = call)
  out_dir <- tasks$built$out_dir
  results <- file.path(out_dir, c("tasks.parquet", "results"))
  .checkArgFlag(overwrite, "overwrite", call = call)
  if (any(file.exists(results)) && !overwrite) {
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

  for (i in seq_len(nrow(listed))) {
    task <- listed[i, ]
    result <- tryCatch(
      withr::with_seed(
        task$seed,
        .fitTask(task, file.path(out_dir, files[[i]])),
        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
        .rng_sample_kind = "Rejection"
      ),
      error = function(e) e
    )

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
    importance = dplyr::bind_rows(lapply(fitted, `[[`, "importance")),
    predictions = dplyr::bind_rows(lapply(fitted, `[[`, "predictions")),
    built = tasks$built
  )

  # If no task was fitted, the results tables still get every column.
  if (!any(listed$status == "fitted")) {
    fit$performance <- tibble::tibble(
      task_id = character(), num_obs = integer(), res_prop = numeric(),
      n_test = integer(), res_prop_test = numeric(), n_feat = integer(),
      model = character(), n_feats_returned = integer(), n_fold = integer(),
      fit_penalty = numeric(), fit_mixture = numeric(), mcc = numeric(), nmcc = numeric(),
      spec = numeric(), sens = numeric(), log2_apop = numeric(), f1 = numeric(),
      bal_acc = numeric(), auprc = numeric(), run_time_sec = numeric(),
      seed = integer(), date = character()
    )
    fit$importance <- tibble::tibble(
      task_id = character(), Variable = character(), Importance = numeric(),
      Sign = character()
    )
    fit$predictions <- tibble::tibble(
      task_id = character(), genome_id = character(), truth = character(),
      predicted = character(), prob_resistant = numeric()
    )
  }

  # Earlier results are replaced only once every task has run.
  unlink(results, recursive = TRUE)
  dir.create(file.path(out_dir, "results"))
  arrow::write_parquet(fit$tasks, file.path(out_dir, "tasks.parquet"))
  for (table in c("performance", "importance", "predictions")) {
    arrow::write_parquet(fit[[table]], file.path(out_dir, "results", paste0(table, ".parquet")))
  }

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
