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

  train <- split(seq_len(nrow(rows)), rows$phenotype) |>
    purrr::map(~ sort(.x[sample.int(length(.x),
                                    floor((1 - holdout) * length(.x)))])) |>
    unlist(use.names = FALSE)
  list(train = train, test = setdiff(seq_len(nrow(rows)), train))
}

# Each training row's fold, assigned per phenotype as rsample's stratified folds do, but
# without pooling rare phenotypes (`pool = 0.1`), which can leave a fold without them.
.stratifiedFolds <- function(phenotype, n_fold) {
  fold <- integer(length(phenotype))
  for (level in intersect(levels(phenotype), phenotype)) {
    i <- which(phenotype == level)
    fold[i] <- sample(rep_len(seq_len(n_fold), length(i)))
  }
  fold
}

# The training data as a data frame: the label, then features renamed x1, x2, ..., so any
# feature ID works in a formula.
.modelData <- function(values, phenotype, sparse = FALSE, labels = .PHENOTYPES) {
  if (sparse) {
    colnames(values) <- paste0("x", seq_len(ncol(values)))
    data <- sparsevctrs::coerce_to_sparse_tibble(values)
  } else {
    data <- as.data.frame(as.matrix(values))
  }
  names(data) <- paste0("x", seq_len(ncol(data)))
  data$phenotype <- factor(phenotype, levels = labels)
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
.modelSpec <- function(multi = FALSE) {
  spec <- if (multi) {
    parsnip::multinom_reg(penalty = tune::tune(), mixture = tune::tune())
  } else {
    parsnip::logistic_reg(penalty = tune::tune(), mixture = tune::tune())
  }
  parsnip::set_engine(spec, "glmnet")
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

# A compact, stable schema for tune notes and conditions outside tuning.
.emptyModelNotes <- function() {
  tibble::tibble(stage = character(), resample = character(), location = character(),
                 type = character(), note = character(), n = integer())
}

.conditionModelNote <- function(cnd, stage) {
  tibble::tibble(stage = stage, resample = NA_character_, location = NA_character_,
                 type = if (inherits(cnd, "warning")) "warning" else "message",
                 note = conditionMessage(cnd), n = 1L)
}

# Ported from old amRml's tuneGrid(), selectBestModel() and fitBestModel(): tune over
# stratified folds, choose by MCC, and fit on all training rows.
.fitModel <- function(data, n_fold, sparse = FALSE, multi = FALSE, cv_workers = 1L) {
  workflow <- workflows::workflow() |>
    workflows::add_recipe(.buildRecipe(data, sparse = sparse)) |>
    workflows::add_model(.modelSpec(multi))

  # After the recipe, whose step IDs use random numbers, so the folds match rsample's.
  fold <- .stratifiedFolds(data$phenotype, n_fold)

  splits <- sort(unique(fold)) |>
    purrr::map(~ rsample::make_splits(
      list(
        analysis = which(fold != .x),
        assessment = which(fold == .x)
      ),
      data = data
    ))
  resamples <- rsample::manual_rset(splits, ids = paste0("Fold", sort(unique(fold))))

  # tune emits CLI alerts and its own progress bar even with verbose = FALSE.
  # Scope these options to this model call; parent progress remains enabled.
  withr::local_options(list(cli.progress_handlers_only = character()))
  tuned <- .tuneModelFolds(workflow, resamples, .tuningGrid(), cv_workers)
  checked <- .checkTuningPaths(tuned)
  tuned <- checked$tuned
  collected <- tune::collect_notes(tuned)
  notes <- checked$notes
  if (nrow(collected)) {
    ids <- grep("^id", names(collected), value = TRUE)
    collected_notes <- tibble::tibble(
      stage = "tune", resample = if (length(ids)) {
        do.call(paste, c(as.list(collected[ids]), sep = "/"))
      } else rep(NA_character_, nrow(collected)),
      location = as.character(collected$location),
      type = as.character(collected$type), note = as.character(collected$note)
    ) |>
      dplyr::count(stage, resample, location, type, note, name = "n")
    notes <- dplyr::bind_rows(notes, collected_notes)
  }
  tryCatch({
    metrics <- suppressWarnings(tune::collect_metrics(tuned))
    if (!any(is.finite(metrics$mean[metrics$.metric == "mcc"]))) {
      stop("No tuning candidate has a defined MCC at a fitted penalty.")
    }
    best <- tune::select_best(tuned, metric = "mcc")
    fit <- withCallingHandlers(
      parsnip::fit(tune::finalize_workflow(workflow, best), data = data),
      warning = function(w) {
        notes <<- dplyr::bind_rows(notes, .conditionModelNote(w, "final_fit"))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        notes <<- dplyr::bind_rows(notes, .conditionModelNote(m, "final_fit"))
        invokeRestart("muffleMessage")
      }
    )
    engine <- parsnip::extract_fit_engine(fit)
    if (!.penaltyCovered(engine$lambda, best$penalty)) {
      path <- .extendPenaltyPath(engine$lambda, best$penalty)
      spec <- tune::finalize_model(.modelSpec(multi), best) |>
        parsnip::set_engine("glmnet", path_values = path)
      final_workflow <- tune::finalize_workflow(workflow, best) |>
        workflows::update_model(spec)
      notes <- dplyr::bind_rows(notes, .conditionModelNote(simpleMessage(
        "Extended final glmnet path to fit the selected penalty exactly."), "final_fit"))
      fit <- withCallingHandlers(parsnip::fit(final_workflow, data = data),
        warning = function(w) {
          notes <<- dplyr::bind_rows(notes, .conditionModelNote(w, "final_fit"))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          notes <<- dplyr::bind_rows(notes, .conditionModelNote(m, "final_fit"))
          invokeRestart("muffleMessage")
        })
      if (!.penaltyCovered(parsnip::extract_fit_engine(fit)$lambda, best$penalty)) {
        stop("The selected penalty was not reached by the final glmnet fit.")
      }
    }
    list(fit = fit, penalty = best$penalty, mixture = best$mixture, notes = notes)
  }, error = function(e) {
    e$amr_notes <- notes
    stop(e)
  })
}

# Ported from old amRml's .calculate*() helpers, keeping AUPRC. Resistant is the event.
# Each metric is rounded to 2, and nMCC and log2(AUPRC / prior) use the rounded values.
.taskMetrics <- function(scored, multi = FALSE) {
  estimate <- function(metric, ...) {
    round(as.numeric(metric(scored, truth = "phenotype", ...)$.estimate), 2)
  }
  mcc <- estimate(yardstick::mcc, estimate = ".pred_class")
  if (multi) {
    return(tibble::tibble(
      mcc = mcc, nmcc = round((mcc + 1) / 2, 2), spec = NA_real_, sens = NA_real_,
      log2_apop = NA_real_, f1 = NA_real_, bal_acc = NA_real_, auprc = NA_real_
    ))
  }
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

# Fun new helper that logs model behavior metrics, retaining confusion matrix values
# and helping us diagnose weird models that collapse by only predicting one class.
.taskDiagnostics <- function(scored, metrics, multi = FALSE) {
  
  collapsed_model <- dplyr::n_distinct(scored$.pred_class) == 1
  
  collapse_class <- if (collapsed_model) {
    as.character(scored$.pred_class[[1]])
  } else {
    NA_character_
  }
  
  prediction_behavior <- if (!collapsed_model) {
    "normal"
  } else if (multi) {
    "constant_class"
  } else if (collapse_class == "Resistant") {
    "all_resistant"
  } else {
    "all_susceptible"
  }
  
  tp <- sum(
    scored$phenotype == "Resistant" &
      scored$.pred_class == "Resistant"
  )
  
  tn <- sum(
    scored$phenotype == "Susceptible" &
      scored$.pred_class == "Susceptible"
  )
  
  fp <- sum(
    scored$phenotype == "Susceptible" &
      scored$.pred_class == "Resistant"
  )
  
  fn <- sum(
    scored$phenotype == "Resistant" &
      scored$.pred_class == "Susceptible"
  )
  
  n_pred_resistant <- sum(
    scored$.pred_class == "Resistant"
  )
  
  n_pred_susceptible <- sum(
    scored$.pred_class == "Susceptible"
  )
  
  tibble::tibble(
    metric_undefined = any(is.na(if (multi) metrics[c("mcc", "nmcc")] else metrics)),
    collapsed_model = collapsed_model,
    collapse_class = collapse_class,
    prediction_behavior = prediction_behavior,
    tp = if (multi) NA_integer_ else tp,
    tn = if (multi) NA_integer_ else tn,
    fp = if (multi) NA_integer_ else fp,
    fn = if (multi) NA_integer_ else fn,
    n_pred_resistant = if (multi) NA_integer_ else n_pred_resistant,
    n_pred_susceptible = if (multi) NA_integer_ else n_pred_susceptible,
    pred_resistant_fraction =
      if (multi) NA_real_ else round(mean(scored$.pred_class == "Resistant"), 4),
    actual_resistant_fraction =
      if (multi) NA_real_ else round(mean(scored$phenotype == "Resistant"), 4)
  )
}

# Ported from old amRml's .viGlmnet(): the non-zero coefficients at the chosen penalty, on
# normalized features. glmnet models the second level, so "POS" points to Susceptible.
.taskImportance <- function(fit, penalty, features, multi = FALSE) {
  table <- function(coefs, label) {
    coefs <- coefs[names(coefs) != "(Intercept)" & coefs != 0]
    coefs <- coefs[order(-abs(coefs))]
    tibble::tibble(
      label = rep(label, length(coefs)),
      Variable = features[as.integer(sub("^x", "", names(coefs)))],
      Importance = abs(unname(coefs)),
      Sign = as.character(ifelse(unname(coefs) > 0, "POS", "NEG"))
    )
  }

  coefs <- stats::coef(parsnip::extract_fit_engine(fit), s = penalty)
  if (!multi) {
    return(table(coefs[, 1], NA_character_))
  }
  dplyr::bind_rows(lapply(names(coefs), function(label) {
    per_label <- coefs[[label]][, 1]
    table(per_label, label)
  }))
}

# Cross-test tasks can share a fit when the training inputs exactly match
.fitTaskGroups <- function(listed, matrices, index) {
  positions <- match(listed$matrix_id, matrices$matrix_id)
  signatures <- vector("list", nrow(listed))
  signature_set <- rep(NA_character_, nrow(listed))

  for (i in which(matrices$mode_id[positions] %in%
                  c("cross_drug", "cross_year", "cross_country"))) {
    p <- positions[[i]]
    rows <- index$rows[[index$matrix_scopes[[p]]]]
    if (!any(rows$role == "test")) next

    train <- rows$role == "train"
    signatures[[i]] <- list(
      file = matrices$file[[p]],
      columns = index$columns[[p]],
      encoding = matrices$encoding[[p]],
      multi = .isMulticlass(matrices$unit[[p]], matrices$target[[p]]),
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
      found <- representatives[
        purrr::map_lgl(
          representatives,
          ~ identical(signatures[[i]], signatures[[.x]])
        )
      ]
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
.fitTask <- function(task, path, species, rows, columns, encoding, cache,
                     fitted = NULL, return_fitted = FALSE, multi = FALSE, cv_workers = 1L) {
  started <- Sys.time()
  if (!multi && !all(rows$phenotype %in% .PHENOTYPES)) {
    stop("A Resistant/Susceptible target has other phenotypes.")
  }
  labels <- if (multi) sort(unique(rows$phenotype), method = "radix") else .PHENOTYPES
  if (task$labels == "shuffled") {
    rows <- .shuffleLabels(rows)
  }
  split <- .splitTask(rows, task$holdout)
  if (is.null(fitted)) {
    matrix <- .readMatrixWide(path, rows, columns, encoding, cache)
    train_values <- matrix$values[split$train, , drop = FALSE]
    sparsity <- 1 - Matrix::nnzero(train_values) /
      (as.double(nrow(train_values)) * ncol(train_values))
    use_sparse <- sparsity >= .SPARSE_THRESHOLD
    train <- .modelData(train_values, rows$phenotype[split$train], sparse = use_sparse, labels = labels)
    test_values <- matrix$values[split$test, , drop = FALSE]
  } else {
    # Identical training inputs were checked by .fitTaskGroups(). Reuse the
    # representation and coefficients; only load/convert the new test rows.
    use_sparse <- fitted$sparse
    matrix <- .readMatrixWide(path, rows[split$test, , drop = FALSE], columns, encoding, cache)
    test_values <- matrix$values
  }
  test <- .modelData(test_values, rows$phenotype[split$test], sparse = use_sparse, labels = labels)
  if (is.null(fitted)) {
    fitted <- .fitModel(train, task$n_fold, sparse = use_sparse, multi = multi,
                        cv_workers = cv_workers)
    fitted$sparse <- use_sparse
    fitted$importance <- .taskImportance(fitted$fit, fitted$penalty, matrix$features, multi = multi)
  }

  scored <- parsnip::augment(fitted$fit, test)
  metrics <- .taskMetrics(scored, multi = multi)
  diagnostics <- .taskDiagnostics(scored, metrics, multi = multi)
  importance <- fitted$importance
  input <- rows$role == "train"

  performance <- tibble::tibble(
    task_id = task$task_id,
    species = species,
    num_obs = sum(input),
    res_prop = if (multi) NA_real_ else round(mean(rows$phenotype[input] == "Resistant"), 2),
    n_test = length(split$test),
    res_prop_test = if (multi) NA_real_ else round(mean(test$phenotype == "Resistant"), 2),
    n_feat = length(matrix$features),
    n_labels = length(labels),
    model = "LR",
    n_feats_returned = dplyr::n_distinct(importance$Variable),
    n_fold = task$n_fold,
    fit_penalty = fitted$penalty,
    fit_mixture = fitted$mixture,
    metrics,
    diagnostics,
    run_time_sec = round(
      as.numeric(difftime(Sys.time(), started, units = "secs")),
      2
    ),
    seed = task$seed,
    date = as.character(Sys.Date())
  )
  
  predictions <- tibble::tibble(
    task_id = task$task_id,
    species = species,
    genome_id = rows$genome_id[split$test],
    truth = as.character(scored$phenotype),
    predicted = as.character(scored$.pred_class),
    prob_resistant = if (multi) NA_real_ else scored$.pred_Resistant,
    prob_predicted = apply(as.matrix(scored[setdiff(grep("^[.]pred_", names(scored), value = TRUE),
                                                    ".pred_class")]), 1, max)
  )

  # Metrics that can't be computed, such as MCC when every test genome gets one prediction.
  computed <- if (multi) metrics[c("mcc", "nmcc")] else metrics
  undefined <- names(computed)[
    purrr::map_lgl(computed, ~ any(is.na(.x)))
  ]
  if (diagnostics$collapsed_model) {
    undefined <- c(
      undefined,
      paste0(
        "every test genome was predicted ",
        diagnostics$collapse_class
      )
    )
  }

  result <- list(
    performance = performance,
    notes = dplyr::bind_cols(
      tibble::tibble(task_id = rep(task$task_id, nrow(fitted$notes %||% .emptyModelNotes())),
                     species = rep(species, nrow(fitted$notes %||% .emptyModelNotes()))),
      fitted$notes %||% .emptyModelNotes()
    ),
    importance = tibble::tibble(
      task_id = rep(task$task_id, nrow(importance)),
      species = rep(species, nrow(importance)),
      importance
    ),
    predictions = predictions,
    undefined = if (length(undefined)) {
      toString(unique(undefined))
    } else {
      NA_character_
    },
    # The model's predictors are x1, x2, ...; `features` names them, in that order.
    model = list(fit = fitted$fit, features = matrix$features)
  )
  if (return_fitted) result$fitted <- fitted
  result
}

## .runRecord block functionality has been folded into run-manifest.R.

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
  units <- table(factor(
    x$built$matrices$unit[match(t$matrix_id, x$built$matrices$matrix_id)],
    c("drug", "drug_class", "mdr")
  ))
  cat("  by target: ", units[["drug"]], " drug, ", units[["drug_class"]],
      " drug class, ", units[["mdr"]], " MDR\n", sep = "")
  cat("  ", nrow(t), " tasks: ", counts(t), "\n", sep = "")
  cat(
    .modeLines(
      names(by_mode),
      purrr::map_chr(by_mode, counts)
    ),
    sep = "\n"
  )
  invisible(x)
}
