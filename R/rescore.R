# Running feature rescoring (R/feature-rescoring.R) on a fit, from its merged results and the ORB.

#' Rescore features
#'
#' Runs feature rescoring ([runFeatureDyadDiscovery()]) on a fit: its vanilla
#' models' merged results, from [exportMergedResults()], and the ORB's
#' virtual `dyad_feature` relation. Feature rescoring needs real and shuffled
#' labels and several seeds, e.g.
#' from `runModelingWorkflow(orb_dir, modes = "vanilla", seeds = c(31, 79, 51),
#' labels = c("real", "shuffled"))`.
#'
#' @param fit An `amr_fit` from [fitModels()] or [runModelingWorkflow()].
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
#'   modes = "vanilla", seeds = c(31, 79, 51), labels = c("real", "shuffled")
#' )
#' crap <- rescoreFeatures(fit)
#' crap$features$top
#' }
#' @export
rescoreFeatures <- function(fit, ...) {
  call <- rlang::current_env()
  .checkArgClass(fit, "fit", "amr_fit", "fitModels()", call = call)

  modes <- .taskModes(fit$tasks, fit$built$matrices)
  vanilla <- fit$tasks[modes %in% "vanilla" & fit$tasks$status == "fitted", ]
  if (!all(c("real", "shuffled") %in% vanilla$labels)) {
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

  # Complete any merged-result export before changing into the ORB directory
  merged <- file.path(
    fit$built$out_dir, "merged", c("all_perf.parquet", "all_top_features.parquet")
  )
  if (!all(file.exists(merged))) {
    exportMergedResults(fit, overwrite = TRUE)
  }

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
