# Running feature rescoring (R/feature-rescoring.R) on a fit, from its native results and the ORB.

#' Rescore features
#'
#' Runs feature rescoring ([runFeatureDyadDiscovery()]) on a fit: its vanilla
#' and MDR models' native performance and label-specific coefficients, and the ORB's
#' virtual `dyad_feature` relation. Feature rescoring needs real and shuffled
#' labels and several seeds, e.g.
#' from `runModelingWorkflow(orb_dir, model_types = c("vanilla", "mdr"), seeds = c(31, 79, 51),
#' labels = c("real", "shuffled"))`.
#'
#' @param fit An `amr_fit` from [fitModels()] or [runModelingWorkflow()].
#' @param min_seeds Minimum matched fitted real/shuffled seeds per vanilla or MDR
#'   target, feature type, and encoding. Default `3L`.
#' @param ... Feature rescoring's thresholds, passed to [runFeatureDyadDiscovery()].
#'
#' @details
#' Binary MDR is scored like a single-drug target. Combination coefficients are
#' normalized, ranked and summarized independently for each \code{label}; model
#' quality is assessed once per fitted model. Binary \code{POS} points toward
#' susceptibility; multinomial \code{POS} points toward its coefficient label.
#' Cross-target overlap compares feature membership without requiring equal signs.
#' The \code{features} and \code{dyads} lists include \code{membership},
#' \code{outcomes}, \code{cross_model_sharedness}, \code{mdr_shared}, and
#' pairwise \code{overlap} tables. Counts distinguish outcome labels from parent
#' model targets and collapse binary/counts copies of the same molecular feature.
#' Overlap covers selected feature sets; absence is not proof of biological absence.
#' \code{models$label_diagnostics} reports one-versus-rest counts and MCC for each
#' combination label and fitted task, separately from overall model quality.
#'
#' @return The `feature_dyad_discovery` list from [runFeatureDyadDiscovery()].
#'
#' @section Errors:
#' Invalid input and missing ORB mappings raise `amrml_invalid_argument`.
#' DuckDB read errors retain their DBI/driver error classes.
#'
#' @examples
#' \dontrun{
#' fit <- runModelingWorkflow("path/to/ORB",
#'   model_types = c("vanilla", "mdr"), seeds = c(31, 79, 51), labels = c("real", "shuffled")
#' )
#' rescored <- rescoreFeatures(fit)
#' rescored$features$top
#' }
#' @export
rescoreFeatures <- function(fit, ..., min_seeds = 3L) {
  call <- rlang::current_env()
  .checkArgClass(fit, "fit", "amr_fit", "fitModels()", call = call)

  .checkArgCount(min_seeds, "min_seeds", call = call)
  scopes <- .rescoreScopes(fit, min_seeds)
  if (!nrow(scopes)) {
    .amrAbort("invalid_argument", paste0(
      "Feature rescoring needs at least ", min_seeds,
      " fitted vanilla or MDR real/shuffled seed pairs per target, feature type, and encoding."
    ), call = call)
  }

  # Native results preserve MDR coefficient labels, compatibility export does not
  selected <- .rescoringInputs(fit, scopes)
  selected_dir <- tempfile("amrml-rescoring-")
  dir.create(selected_dir)
  on.exit(unlink(selected_dir, recursive = TRUE), add = TRUE)
  merged <- file.path(selected_dir, c("performance.parquet", "importance.parquet"))
  .writeParquet(selected$performance, merged[[1L]])
  .writeParquet(selected$importance, merged[[2L]])

  # Locate the ORB
  orb_files <- fit$built$eligibility$profile$orb$files$path
  orb_dirs <- unique(dirname(orb_files))
  orb_databases <- unique(unlist(lapply(orb_dirs, function(dir) {
    list.files(dir, pattern = "_parquet[.]duckdb$", full.names = TRUE)
  }), use.names = FALSE))

  if (length(orb_databases) != 1L) {
    .amrAbort(
      "invalid_argument",
      "Expected exactly one Parquet-backed ORB DuckDB with a virtual `dyad_feature` view.",
      call = call
    )
  }

  # Virtual Parquet paths are relative to the ORB directory; restore caller's 
  # working directory after database disconnect
  previous_wd <- getwd()
  on.exit(setwd(previous_wd), add = TRUE)
  setwd(dirname(orb_databases[[1L]]))
  con <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = basename(orb_databases[[1L]]),
    read_only = TRUE
  )
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (!"dyad_feature" %in% DBI::dbListTables(con) ||
      !all(c("dyad", "feature") %in% DBI::dbListFields(con, "dyad_feature"))) {
    .amrAbort(
      "invalid_argument",
      "ORB DuckDB must expose the virtual `dyad_feature(dyad, feature)` relation.",
      call = call
    )
  }

  result <- runFeatureDyadDiscovery(
    all_top_features_parquet = merged[[2]],
    all_performance_parquet = merged[[1]],
    dyad_con = con,
    ...
  )
  # Label-level diagnostics accompany model-level MCC
  result$models$label_diagnostics <- .mdrLabelDiagnostics(fit)
  result
}

# Count successfully fitted seed pairs within each vanilla or MDR representation
.rescoreScopes <- function(fit, min_seeds = 3L) {
  metadata <- .exportMetadata(fit)
  metadata$matrix_id <- fit$tasks$matrix_id
  metadata <- metadata[metadata$mode_id %in% c("vanilla", "mdr") &
                         fit$tasks$status == "fitted", ]
  metadata |>
    dplyr::group_by(matrix_id, seed) |>
    dplyr::summarise(paired = all(c(FALSE, TRUE) %in% shuffled), .groups = "drop") |>
    dplyr::filter(paired) |>
    dplyr::count(matrix_id, name = "n_paired_seeds") |>
    dplyr::filter(n_paired_seeds >= min_seeds)
}

# Empty outcomes are distinct from failures!
.noRescoringEvidence <- function(message) {
  rlang::abort(message, class = "amrml_no_rescoring_evidence")
}

.rescoreWorkflow <- function(fit, manifest_path, enabled = TRUE, reason = NULL) {
  scopes <- .rescoreScopes(fit)
  status <- "skipped"
  result <- NULL
  outputs <- character()
  if (enabled && !nrow(scopes)) {
    reason <- "MAGIC needs at least 3 fitted vanilla or MDR real/shuffled seed pairs per target, feature type, and encoding."
  }
  if (enabled && nrow(scopes)) {
    .manifestArtifact(manifest_path, "feature_rescoring", status = "building")
    .manifestStage(manifest_path, "rescore_features", status = "running")
    result <- tryCatch(rescoreFeatures(fit), amrml_no_rescoring_evidence = function(e) {
      reason <<- conditionMessage(e)
      NULL
    })
    status <- if (is.null(result)) "empty" else "success"
    if (!is.null(result)) {
      directory <- file.path(.runPaths(fit$built$out_dir)$results, "rescoring")
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
      rds <- file.path(directory, "discovery.rds")
      saveRDS(result, rds)
      outputs <- rds
      for (group in c("models", "features", "dyads")) {
        for (name in names(result[[group]])) {
          path <- file.path(directory, paste0(group, "_", name, ".parquet"))
          .writeParquet(result[[group]][[name]], path)
          outputs <- c(outputs, path)
        }
      }
      reason <- "MAGIC feature/dyad rescoring completed."
    }
  }
  .manifestStage(manifest_path, "rescore_features", status = status,
                 parameters = list(min_paired_seeds = 3L), outputs = outputs,
                 metrics = list(n_eligible_scopes = nrow(scopes)), message = reason)
  .manifestArtifact(manifest_path, "feature_rescoring",
                    status = if (status == "success") "ready" else status,
                    details = list(directory = if (status == "success") "results/rescoring" else NULL,
                                   path_mode = "relative_to_manifest", reason = reason,
                                   min_paired_seeds = 3L,
                                   scoring = "Returned-feature ranking; candidate-aware scoring is not applied."))
  .manifestEvent(manifest_path, reason)
  message(reason)
  fit$rescoring <- result
  fit$rescoring_status <- list(status = status, reason = reason)
  fit
}

.rescoringInputs <- function(fit, scopes) {
  metadata <- .exportMetadata(fit)
  metadata$matrix_id <- fit$tasks$matrix_id
  metadata <- metadata[metadata$matrix_id %in% scopes$matrix_id &
                         fit$tasks$status == "fitted", ]
  pairs <- metadata |>
    dplyr::group_by(matrix_id, seed) |>
    dplyr::filter(all(c(FALSE, TRUE) %in% shuffled)) |>
    dplyr::ungroup()
  layout <- .EXPORT_LAYOUTS$vanilla
  performance <- fit$performance[names(fit$performance) != "seed"]
  importance <- fit$importance
  if (!"label" %in% names(importance)) importance$label <- NA_character_
  multiclass_tasks <- pairs$task_id[.isMulticlass(pairs$unit, pairs$target)]
  if (anyNA(importance$label[importance$task_id %in% multiclass_tasks])) {
    .amrAbort("invalid_argument", "Multiclass coefficients must retain their outcome labels.")
  }
  # n_feats_returned = distinct predictors
  counts <- importance |>
    dplyr::group_by(task_id) |>
    dplyr::summarise(n_distinct_features = dplyr::n_distinct(Variable), .groups = "drop")
  at <- match(performance$task_id, counts$task_id)
  performance$n_feats_returned <- dplyr::coalesce(counts$n_distinct_features[at], 0L)
  list(
    performance = .exportRows(performance, pairs, layout$perf),
    importance = .exportRows(importance, pairs, c(layout$top, "label"))
  )
}

# One-versus-rest diagnostics for multiclass labels. These are held-out
# prediction summaries, not independent fits or a claim of causal importance.
.mdrLabelDiagnostics <- function(fit) {
  metadata <- .exportMetadata(fit)
  metadata$matrix_id <- fit$tasks$matrix_id
  metadata <- metadata[.isMulticlass(metadata$unit, metadata$target) &
                         fit$tasks$status == "fitted", ]
  predictions <- dplyr::inner_join(fit$predictions, metadata, by = "task_id",
                                   relationship = "many-to-one")
  empty <- tibble::tibble(task_id = character(), species = character(),
    unit = character(), target = character(), seed = integer(), shuffled = logical(),
    feature_type = character(), encoding = character(), label = character(),
    n_truth = integer(), n_predicted = integer(), tp = integer(), tn = integer(),
    fp = integer(), fn = integer(), sensitivity = numeric(), precision = numeric(),
    label_mcc = numeric())
  if (!nrow(predictions)) return(empty)
  inventory <- dplyr::bind_rows(
    dplyr::select(predictions, matrix_id, label = truth),
    dplyr::select(predictions, matrix_id, label = predicted)
  ) |>
    dplyr::distinct()
  pieces <- lapply(split(predictions, predictions$task_id), function(rows) {
    labels <- inventory$label[inventory$matrix_id == rows$matrix_id[[1L]]]
    dplyr::bind_rows(lapply(labels, function(label_value) {
      truth_positive <- rows$truth == label_value
      predicted_positive <- rows$predicted == label_value
      tp <- sum(truth_positive & predicted_positive)
      tn <- sum(!truth_positive & !predicted_positive)
      fp <- sum(!truth_positive & predicted_positive)
      fn <- sum(truth_positive & !predicted_positive)
      denominator <- sqrt(as.double(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
      tibble::tibble(task_id = rows$task_id[[1L]], species = rows$species[[1L]],
        unit = rows$unit[[1L]], target = rows$target[[1L]], seed = rows$seed[[1L]],
        shuffled = rows$shuffled[[1L]], feature_type = rows$feature_type[[1L]],
        encoding = rows$encoding[[1L]], label = label_value,
        n_truth = tp + fn, n_predicted = tp + fp, tp = tp, tn = tn, fp = fp, fn = fn,
        sensitivity = if (tp + fn > 0L) tp / (tp + fn) else NA_real_,
        precision = if (tp + fp > 0L) tp / (tp + fp) else NA_real_,
        label_mcc = if (denominator > 0) (as.double(tp) * tn - as.double(fp) * fn) /
          denominator else NA_real_)
    }))
  })
  dplyr::bind_rows(c(list(empty), pieces))
}
