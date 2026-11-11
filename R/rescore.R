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
  needed <- c(
    "workflow.rds", "tasks.parquet",
    file.path("results", c(
      "fit-completion.rds", "performance.parquet", "importance.parquet", "predictions.parquet"
    ))
  )
  missing <- needed[!file.exists(file.path(run_dir, needed))]
  if (length(missing)) {
    .amrAbort(
      "invalid_argument",
      c(
        "Not a completed amRml run folder.",
        x = run_dir,
        x = paste0("Missing: ", paste(missing, collapse = ", "))
      ),
      call = call
    )
  }

  built <- readRDS(file.path(run_dir, "workflow.rds"))$built
  # The folder may have moved since the run.
  built$out_dir <- run_dir
  read <- function(...) arrow::read_parquet(file.path(run_dir, ...))
  structure(
    list(
      tasks = read("tasks.parquet"),
      performance = read("results", "performance.parquet"),
      importance = read("results", "importance.parquet"),
      predictions = read("results", "predictions.parquet"),
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
        "Run e.g. `runModelingWorkflow(orb_dir, modes = \"vanilla\", ",
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

  # Zstandard compresses these tables well; Snappy is the fallback where an Arrow build lacks it.
  compression <- if (arrow::codec_is_available("zstd")) "zstd" else "snappy"
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
#' from `runModelingWorkflow(orb_dir, modes = "vanilla", seeds = c(31, 79, 51),
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
#' @param ... Feature rescoring's thresholds, passed to [runFeatureDyadDiscovery()].
#'
#' @return The `feature_dyad_discovery` list from [runFeatureDyadDiscovery()]. Its
#'   `"output_dir"` attribute is the folder the results were written to.
#'
#' @section Output files:
#' Every call writes its results; there is no option to skip it. Each call writes a new folder, `<run>/magic/<date-time>/`, so a result is never replaced by a
#' later call with other thresholds. It holds Parquet files, compressed with Zstandard:
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
#'   modes = "vanilla", seeds = c(31, 79, 51), labels = c("real", "shuffled")
#' )
#' rescored <- rescoreFeatures(fit)
#' rescored$features$top
#'
#' # Without one: the newest rescorable run of an ORB, or of a registered ORB
#' rescoreFeatures(orb_dir = "path/to/ORB")
#' rescoreFeatures()
#' }
#' @export
rescoreFeatures <- function(fit = NULL, orb_dir = NULL, run_dir = NULL, ...) {
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
    if (is.null(orb_dir)) orb_dir <- .selectORB(call)
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
          "Run e.g. `runModelingWorkflow(orb_dir, modes = \"vanilla\", ",
          "seeds = c(31, 79, 51), labels = c(\"real\", \"shuffled\"))`."
        )
      ),
      call = call
    )
  }

  # Resolve the run folder before changing into the ORB directory below
  run_folder <- normalizePath(fit$built$out_dir, winslash = "/", mustWork = TRUE)

  # Complete any merged-result export before changing into the ORB directory
  merged <- file.path(
    fit$built$out_dir, "merged", c("all_perf.parquet", "all_top_features.parquet")
  )
  if (!all(file.exists(merged))) {
    exportMergedResults(fit, overwrite = TRUE)
  }

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
