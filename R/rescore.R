# Running feature rescoring (R/feature-rescoring.R) on a fit, from its merged results and the ORB.

# Whether a fit has fitted vanilla models, other than MDR, with both real and shuffled labels,
# which feature rescoring needs.
.hasRescorableFits <- function(fit) {
  modes <- .taskModes(fit$tasks, fit$built$matrices)
  unit <- fit$built$matrices$unit[match(fit$tasks$matrix_id, fit$built$matrices$matrix_id)]
  vanilla <- fit$tasks[modes %in% "vanilla" & unit != "mdr" & fit$tasks$status == "fitted", ]
  all(c("real", "shuffled") %in% vanilla$labels)
}

# A completed run folder, reopened as an `amr_fit` without refitting. A folder is complete
# once the fit's results and completion receipt are written beside `workflow.rds`.
.readFitRun <- function(run_dir, call) {
  if (!rlang::is_string(run_dir) || !dir.exists(run_dir)) {
    .amrAbort("invalid_argument", "`run_dir` must be an existing amRml run folder.", call = call)
  }
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  paths <- .runPaths(run_dir)
  needed <- c(paths$workflow, paths$completion, paths$tasks, paths$performance,
              paths$importance, paths$predictions)
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    .amrAbort(
      "invalid_argument",
      c(
        "Not a completed amRml run folder.",
        x = run_dir,
        x = paste0("Missing: ", paste(substring(missing, nchar(run_dir) + 2L), collapse = ", ")),
        i = "Resume state is removed after a successful run unless `save_state = TRUE`."
      ),
      call = call
    )
  }

  built <- readRDS(paths$workflow)$built
  # The folder may have moved since the run.
  built$out_dir <- run_dir
  read <- function(path) arrow::read_parquet(path)
  structure(
    list(
      tasks = read(paths$tasks),
      performance = read(paths$performance),
      importance = read(paths$importance),
      predictions = read(paths$predictions),
      built = built
    ),
    class = "amr_fit"
  )
}

# The newest run under `<ORB>/amRml/` that is complete and can be rescored.
.findFitRun <- function(orb_directory, call) {
  runs <- list.dirs(file.path(orb_directory, "amRml"), recursive = FALSE)
  # Run folders are named by date and time, so name order is run order.
  for (run in sort(runs, decreasing = TRUE, method = "radix")) {
    fit <- tryCatch(.readFitRun(run, call), error = function(e) NULL)
    if (!is.null(fit) && .hasRescorableFits(fit)) return(fit)
  }
  .amrAbort(
    "invalid_argument",
    c(
      "No completed run in this ORB can be rescored.",
      x = file.path(orb_directory, "amRml"),
      i = paste0(
        "Run e.g. `runModelingWorkflow(orb_dir, model_types = \"vanilla\", ",
        "seeds = c(31, 79, 51), labels = c(\"real\", \"shuffled\"))`."
      )
    ),
    call = call
  )
}

# The result's 16 tables, kept to four. `consistent`, `all_seeds`, `unique`, `shared` and the
# `high_*` tables are subsets of `features$top` and `dyads$top`, so they become flag columns
# on one table each instead of repeated rows.
.rescoredTables <- function(result) {
  min_shared <- result$parameters$min_models_shared
  in_subset <- function(x, subset, keys) {
    do.call(paste, c(unname(as.list(x[keys])), sep = "\r")) %in%
      do.call(paste, c(unname(as.list(subset[keys])), sep = "\r"))
  }

  fit_keys <- c(
    "species", "drug_label", "drug_or_class", "seed", "feature_type", "feature_subtype",
    "model", "fit_penalty", "fit_mixture"
  )
  fits <- dplyr::bind_rows(
    result$models$good_performance,
    result$models$good_shuffled_separation,
    result$models$qualified
  )
  fits <- dplyr::distinct(fits, dplyr::across(dplyr::all_of(fit_keys)), .keep_all = TRUE)
  fits$is_good_performance <- in_subset(fits, result$models$good_performance, fit_keys)
  fits$is_separated_from_shuffled <- in_subset(fits, result$models$good_shuffled_separation, fit_keys)
  fits$is_qualified <- in_subset(fits, result$models$qualified, fit_keys)

  feature_keys <- c("species", "drug_label", "drug_or_class", "feature_type", "feature_subtype", "variable")
  features <- dplyr::left_join(
    result$features$top,
    result$features$sharedness[c(
      "species", "drug_label", "feature_type", "variable", "n_drug_or_class", "drug_or_class_csv"
    )],
    by = c("species", "drug_label", "feature_type", "variable")
  )
  features$is_consistent <- in_subset(features, result$features$consistent, feature_keys)
  features$is_high_contribution <- in_subset(features, result$features$high_contribution, feature_keys)
  features$is_all_seeds <- in_subset(features, result$features$all_seeds, feature_keys)
  features$is_shared <- features$n_drug_or_class >= min_shared
  features$is_unique <- features$n_drug_or_class == 1L

  dyad_keys <- c("species", "drug_label", "drug_or_class", "dyad")
  dyads <- dplyr::left_join(
    result$dyads$top,
    result$dyads$sharedness[c("species", "drug_label", "dyad", "n_drug_or_class", "drug_or_class_csv")],
    by = c("species", "drug_label", "dyad")
  )
  dyads$is_high_score <- in_subset(dyads, result$dyads$high_score, dyad_keys)
  dyads$is_shared <- dyads$n_drug_or_class >= min_shared
  dyads$is_unique <- dyads$n_drug_or_class == 1L

  list(
    optimal_models = fits,
    model_summary = result$models$summary,
    qualified_features = features,
    qualified_dyads = dyads
  )
}

# Zstandard compresses these tables well; Snappy is the fallback where an Arrow build lacks it.
.parquetCompression <- function() {
  if (arrow::codec_is_available("zstd")) "zstd" else "snappy"
}

# The merged performance and top-feature tables with every per-seed calculation that does not
# depend on a threshold added as columns. Thresholds (the core share, the MCC cut-off, the rank
# quantile and the rest) are then filters over these columns. Original columns are unchanged.
.completeRescoring <- function(perf_parquet, top_parquet) {
  perf <- arrow::read_parquet(perf_parquet)
  top <- arrow::read_parquet(top_parquet)

  fit_keys <- c("species", "drug_label", "drug_or_class", "feature_type", "feature_subtype", "seed")

  # Real-label MCC against its shuffled counterpart, and the shuffled background it is judged by.
  # Defined only on real-label rows; the shuffled rows keep NA.
  shuffled_mcc <- perf[perf$shuffled, c(fit_keys, "mcc")]
  names(shuffled_mcc)[names(shuffled_mcc) == "mcc"] <- "shuffled_mcc"
  background <- assessModelReliability(perf_parquet)$by_seed[c(
    fit_keys, "mean_shuffled_mcc", "sd_shuffled_mcc", "shuffled_upper", "above_background"
  )]
  perf$.row <- seq_len(nrow(perf))
  perf <- dplyr::left_join(perf, dplyr::mutate(shuffled_mcc, shuffled = FALSE), by = c(fit_keys, "shuffled"))
  perf <- dplyr::left_join(perf, dplyr::mutate(background, shuffled = FALSE), by = c(fit_keys, "shuffled"))
  perf$mcc_diff <- perf$mcc - perf$shuffled_mcc
  perf$feat_return_ratio <- perf$n_feats_returned / perf$n_feat
  perf <- perf[order(perf$.row), ]
  perf$.row <- NULL

  # Contribution, rank, rank score and cumulative share within each model's seed. The core group
  # is not stored: it follows from `cum_contrib - contribution < core_contribution_threshold`.
  scored <- scoreFeaturesWithinSeed(
    top_parquet, filter_model = FALSE, all_performance_parquet = perf_parquet
  )
  scored <- scored[c(
    fit_keys, "variable", "contribution", "rank", "n_features", "rank_score", "cum_contrib"
  )]
  # Not `variable`, which a case-insensitive reader such as DuckDB would confuse with `Variable`
  names(scored)[names(scored) == "variable"] <- "variable_clean"
  scored$shuffled <- FALSE
  top$variable_clean <- .cleanFeatureNames(top$feature_type, top$Variable)
  join_keys <- c(fit_keys, "shuffled", "variable_clean")
  if (anyDuplicated(scored[join_keys])) {
    rlang::abort("Scored features are not unique within each seed, so they cannot be matched back.")
  }
  rows <- nrow(top)
  top <- dplyr::left_join(top, scored, by = join_keys)
  if (nrow(top) != rows) rlang::abort("Matching scores back to the features changed the row count.")

  list(all_perf = tibble::as_tibble(perf), all_top_features = tibble::as_tibble(top))
}

# Write the complete tables to `<run_dir>/magic/`, unless the copies there are newer than the
# merged results they came from. Returns the folder.
.writeCompleteRescoring <- function(merged, run_dir) {
  out_dir <- file.path(run_dir, "magic")
  files <- file.path(out_dir, basename(merged))
  if (all(file.exists(files)) && min(file.mtime(files)) >= max(file.mtime(merged))) {
    return(out_dir)
  }

  tables <- .completeRescoring(merged[[1L]], merged[[2L]])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Write beside the target and rename, so an interrupted write never leaves a half file.
  for (name in names(tables)) {
    target <- file.path(out_dir, paste0(name, ".parquet"))
    partial <- paste0(target, ".partial")
    arrow::write_parquet(tables[[name]], partial, compression = .parquetCompression())
    file.rename(partial, target)
  }
  rlang::inform(paste0("Per-seed calculations written to: ", out_dir))
  out_dir
}

# Write a result to `<run_dir>/magic/<date-time>/`, a new folder for every call so an earlier
# result is never replaced. Returns that folder.
.writeRescored <- function(result, run_dir) {
  tables <- .rescoredTables(result)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  out_dir <- file.path(run_dir, "magic", stamp)
  n <- 1L
  while (dir.exists(out_dir)) {
    n <- n + 1L
    out_dir <- file.path(run_dir, "magic", paste0(stamp, "-", n))
  }
  dir.create(out_dir, recursive = TRUE)

  compression <- .parquetCompression()
  for (name in names(tables)) {
    arrow::write_parquet(
      tables[[name]], file.path(out_dir, paste0(name, ".parquet")), compression = compression
    )
  }
  jsonlite::write_json(
    list(
      run = basename(run_dir),
      written = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      amRml_version = tryCatch(
        as.character(utils::packageVersion("amRml")), error = function(e) NA_character_
      ),
      compression = compression,
      rows = lapply(tables, nrow),
      parameters = result$parameters
    ),
    file.path(out_dir, "parameters.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"
  )
  normalizePath(out_dir, winslash = "/", mustWork = TRUE)
}

#' Rescore features
#'
#' Runs feature rescoring ([runFeatureDyadDiscovery()]) on a fit: its vanilla
#' models' merged results, from [exportMergedResults()], and the ORB's
#' virtual `dyad_feature` relation. Feature rescoring needs real and shuffled
#' labels and several seeds, e.g.
#' from `runModelingWorkflow(orb_dir, model_types = "vanilla", seeds = c(31, 79, 51),
#' labels = c("real", "shuffled"))`.
#'
#' @details
#' No fit object is needed. With none of `fit`, `run_dir` and `orb_dir`, an ORB is
#' chosen from those amRdata registered, as in [readORB()] (this needs an
#' interactive session). The newest completed run in the ORB's `amRml/` folder with
#' fitted vanilla models on both real and shuffled labels is then reopened from disk
#' and rescored. The run used is reported.
#'
#' @param fit An `amr_fit` from [fitModels()] or [runModelingWorkflow()]. Optional.
#' @param orb_dir An amRdata ORB directory, whose newest rescorable run is used.
#'   Used when `fit` and `run_dir` are not supplied; if it is `NULL` too, an ORB is chosen
#'   from the registry.
#' @param run_dir A completed run folder (e.g. `<ORB>/amRml/<date-time>/`) to
#'   reopen. Supply at most one of `fit`, `run_dir` and `orb_dir`.
#' @param min_seeds Minimum matched fitted real/shuffled seeds per vanilla
#'   drug, feature type, and encoding. Default `3L`.
#' @param ... Feature rescoring's thresholds, passed to [runFeatureDyadDiscovery()].
#'
#' @return The `feature_dyad_discovery` list from [runFeatureDyadDiscovery()]. Its
#'   `"output_dir"` attribute is the folder the results were written to.
#'
#' @section Output files:
#' Every call writes its results; there is no option to skip it. Everything goes in
#' `<run>/magic/`, as Parquet files compressed with Zstandard.
#'
#' The per-seed calculations that need no threshold are saved once for the run, as
#' `all_perf.parquet` and `all_top_features.parquet`. They are the merged results with columns
#' added, and every original column is unchanged. They are rewritten only when the merged
#' results are newer. Thresholds are then filters over these columns.
#' \itemize{
#'   \item `all_perf.parquet` adds `shuffled_mcc`, `mcc_diff`, `feat_return_ratio`, and the
#'   shuffled background from [assessModelReliability()]: `mean_shuffled_mcc`,
#'   `sd_shuffled_mcc`, `shuffled_upper` and `above_background`. These are set on real-label
#'   rows, and are `NA` where a seed has no shuffled counterpart.
#'   \item `all_top_features.parquet` adds `variable_clean` (the name used by the ORB's
#'   dyad relation), `contribution`, `rank`, `n_features`, `rank_score` and `cum_contrib`, set on
#'   real-label rows. A feature is in the core group when
#'   `cum_contrib - contribution < core_contribution_threshold`.
#' }
#'
#' Each call also writes its own result to a new folder, `<run>/magic/<date-time>/`, so a result
#' is never replaced by a later call with other thresholds. It holds:
#' \itemize{
#'   \item `optimal_models.parquet`: every fit that passes at least one model criterion, with
#'   `is_good_performance`, `is_separated_from_shuffled` and `is_qualified` flags.
#'   \item `model_summary.parquet`: one row per species and drug or drug class.
#'   \item `qualified_features.parquet`: the selected features, with `is_consistent`,
#'   `is_high_contribution`, `is_all_seeds`, `is_shared` and `is_unique` flags, and the number of
#'   drugs or classes each feature is selected for.
#'   \item `qualified_dyads.parquet`: the dyads reached by those features, with `is_high_score`, `is_shared`
#'   and `is_unique` flags.
#'   \item `parameters.json`: the thresholds used, the run, the package version and row counts.
#' }
#' The result's subset tables (for example `features$consistent`) are flags on these tables, not
#' separate files, to keep the folder small. The dyads table is the largest and grows with the
#' number of selected features.
#'
#' @section Errors:
#' Invalid input, an unusable or incomplete run folder, no rescorable run, and
#' missing ORB mappings raise `amrml_invalid_argument`.
#' DuckDB read errors retain their DBI/driver error classes.
#'
#' @examples
#' \dontrun{
#' # From a fit in memory
#' fit <- runModelingWorkflow("path/to/ORB",
#'   model_types = "vanilla", seeds = c(31, 79, 51), labels = c("real", "shuffled")
#' )
#' rescored <- rescoreFeatures(fit)
#' rescored$features$top
#'
#' # Without one: the newest rescorable run of an ORB, or of a registered ORB
#' rescoreFeatures(orb_dir = "path/to/ORB")
#' rescoreFeatures()
#' }
#' @export
rescoreFeatures <- function(fit = NULL, orb_dir = NULL, run_dir = NULL, ..., min_seeds = 3L) {
  call <- rlang::current_env()
  if (sum(!is.null(fit), !is.null(orb_dir), !is.null(run_dir)) > 1L) {
    .amrAbort(
      "invalid_argument",
      "Supply at most one of `fit`, `run_dir` and `orb_dir`.",
      call = call
    )
  }

  orb_directory <- NULL
  if (!is.null(fit)) {
    .checkArgClass(fit, "fit", "amr_fit", "fitModels()", call = call)
  } else if (!is.null(run_dir)) {
    fit <- .readFitRun(run_dir, call)
  } else {
    # No fit given: find the ORB (from the registry if need be), then its newest usable run.
    if (is.null(orb_dir)) orb_dir <- .selectORB(call, argument = "orb_dir")
    orb_directory <- readORB(orb_dir)$directory
    fit <- .findFitRun(orb_directory, call)
    rlang::inform(paste0("Rescoring run: ", fit$built$out_dir))
  }

  if (!.hasRescorableFits(fit)) {
    .amrAbort(
      "invalid_argument",
      c(
        "Feature rescoring needs fitted vanilla models with both real and shuffled labels.",
        i = paste0(
          "Run e.g. `runModelingWorkflow(orb_dir, model_types = \"vanilla\", ",
          "seeds = c(31, 79, 51), labels = c(\"real\", \"shuffled\"))`."
        )
      ),
      call = call
    )
  }
  .checkArgCount(min_seeds, "min_seeds", call = call)
  scopes <- .rescoreScopes(fit, min_seeds)
  if (!nrow(scopes)) {
    .amrAbort("invalid_argument", paste0(
      "Feature rescoring needs at least ", min_seeds,
      " fitted vanilla real/shuffled seed pairs per drug, feature type, and encoding."
    ), call = call)
  }

  # Resolve the run folder before changing into the ORB directory below
  run_folder <- normalizePath(fit$built$out_dir, winslash = "/", mustWork = TRUE)

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
  orb_dirs <- orb_directory %||%
    unique(dirname(fit$built$eligibility$profile$orb$files$path))
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

  out_dir <- .writeRescored(result, run_folder)
  attr(result, "output_dir") <- out_dir
  rlang::inform(paste0("Rescoring results written to: ", out_dir))
  result
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
