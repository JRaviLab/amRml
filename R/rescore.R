# Running feature rescoring (R/feature-rescoring.R) on a fit, from its merged results and the ORB.

#' Rescore features
#'
#' Runs feature rescoring ([runFeatureDyadDiscovery()]) on a fit: its vanilla
#' models' merged results, from [exportMergedResults()], and the ORB's
#' protein-dyad map (`dyad_feature.parquet`). Feature rescoring needs real and
#' shuffled labels and several seeds, e.g. from `runModelingWorkflow(orb_dir,
#' modes = "vanilla", seeds = c(31, 79, 51), labels = c("real", "shuffled"))`.
#'
#' @param fit An `amr_fit` from [fitModels()] or [runModelingWorkflow()].
#' @param ... Feature rescoring's thresholds, passed to [runFeatureDyadDiscovery()].
#'
#' @return The `feature_dyad_discovery` list from [runFeatureDyadDiscovery()].
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus `amrml_invalid_argument`: `fit` is
#' not from [fitModels()], it has no fitted vanilla models with both real and
#' shuffled labels, or the ORB declares no `dyad_feature.parquet`.
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

  orb_files <- fit$built$eligibility$profile$orb$files$path
  dyads <- orb_files[basename(orb_files) == "dyad_feature.parquet"]
  if (!length(dyads)) {
    .amrAbort(
      "invalid_argument",
      "The ORB declares no `dyad_feature.parquet`, which maps features to protein dyads.",
      call = call
    )
  }

  merged <- file.path(
    fit$built$out_dir, "merged", c("all_perf.parquet", "all_top_features.parquet")
  )
  if (!all(file.exists(merged))) {
    exportMergedResults(fit, overwrite = TRUE)
  }

  runFeatureDyadDiscovery(
    all_top_features_parquet = merged[[2]],
    all_performance_parquet = merged[[1]],
    dyad_feature_parquet = dyads[[1]],
    ...
  )
}
