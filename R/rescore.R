# Running feature rescoring (R/feature-rescoring.R) on a fit, from its merged results and the ORB.

#' Rescore features
#'
#' Runs feature rescoring ([runFeatureDyadDiscovery()]) on a fit: its vanilla
#' models' merged results, from [exportMergedResults()], and the ORB's
#' virtual `dyad_feature` relation. Feature rescoring needs real and shuffled
#' labels and several seeds, e.g.
#' from `runModelingWorkflow(orb_dir, model_types = "vanilla", seeds = c(31, 79, 51),
#' labels = c("real", "shuffled"))`.
#'
#' @param fit An `amr_fit` from [fitModels()] or [runModelingWorkflow()].
#' @param min_seeds Minimum matched fitted real/shuffled seeds per vanilla
#'   drug, feature type, and encoding. Default `3L`.
#' @param ... Feature rescoring's thresholds, passed to [runFeatureDyadDiscovery()].
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
#'   model_types = "vanilla", seeds = c(31, 79, 51), labels = c("real", "shuffled")
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
      " fitted vanilla real/shuffled seed pairs per drug, feature type, and encoding."
    ), call = call)
  }

  # Complete any merged-result export before changing into the ORB directory
  merged <- file.path(
    .runPaths(fit$built$out_dir)$merged,
    c("all_perf.parquet", "all_top_features.parquet")
  )
  if (!all(file.exists(merged))) {
    exportMergedResults(fit, overwrite = TRUE)
  }

  # Use only representations with sufficient paired & fitted seeds
  scope_keys <- c("species", "drug_label", "drug_or_class", "feature_type", "feature_subtype")
  metadata <- .exportMetadata(fit)
  metadata$matrix_id <- fit$tasks$matrix_id
  eligible <- metadata[metadata$matrix_id %in% scopes$matrix_id, ]
  eligible <- eligible[!duplicated(eligible[c("unit", "target", "feature_type", "encoding")]), ]
  allowed <- tibble::tibble(
    drug_label = eligible$unit, drug_or_class = eligible$target,
    feature_type = eligible$feature_type, feature_subtype = eligible$encoding
  )
  selected <- lapply(merged, function(path) {
    dplyr::semi_join(arrow::read_parquet(normalizePath(path)), allowed,
                     by = scope_keys[-1L])
  })
  selected_dir <- tempfile("amrml-rescoring-")
  dir.create(selected_dir)
  on.exit(unlink(selected_dir, recursive = TRUE), add = TRUE)
  merged <- file.path(selected_dir, basename(merged))
  for (i in seq_along(merged)) .writeParquet(selected[[i]], merged[[i]])

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

  runFeatureDyadDiscovery(
    all_top_features_parquet = merged[[2]],
    all_performance_parquet = merged[[1]],
    dyad_con = con,
    ...
  )
}

# Count successfully fitted seed pairs within each vanilla representation
.rescoreScopes <- function(fit, min_seeds = 3L) {
  metadata <- .exportMetadata(fit)
  metadata$matrix_id <- fit$tasks$matrix_id
  metadata <- metadata[metadata$mode_id == "vanilla" & metadata$unit != "mdr" &
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
    reason <- "MAGIC needs at least 3 fitted vanilla real/shuffled seed pairs per drug, feature type, and encoding."
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
