# Temporary: a fit's results in old amRml's merged format, until amRml's column names match.

# Old amRml's performance columns, without and with a separate test set.
.OLD_PERF <- c(
  "num_obs", "res_prop", "n_feat", "model", "n_feats_returned", "n_fold", "fit_penalty",
  "fit_mixture", "mcc", "nmcc", "spec", "sens", "log2_apop", "f1", "bal_acc", "run_time_sec",
  "seed", "date"
)
.OLD_PERF_TEST <- c(
  "num_obs_original_ml_input_tibble", "num_obs_test_data", "res_prop_original_ml_input_tibble",
  "n_feat", "res_prop_test_data", "model", "n_feats_returned", "n_fold", "fit_penalty",
  "fit_mixture", "mcc", "nmcc", "spec", "sens", "log2_apop", "f1", "bal_acc", "run_time_sec",
  "seed", "date"
)
.OLD_TOP <- c("Variable", "Importance", "Sign")

# Each file layout's columns, in old amRml's order; `top = NULL` where old amRml merged no top
# features. `shuffled` and an integer `seed` are ours wherever old amRml had none.
.EXPORT_LAYOUTS <- list(
  vanilla = list(
    perf = c(
      "shuffled", "species", "drug_label", "drug_or_class", "feature_type", "feature_subtype",
      .OLD_PERF
    ),
    top = c(
      "shuffled", "species", "drug_label", "drug_or_class", "feature_type", "feature_subtype",
      "seed", .OLD_TOP
    )
  ),
  stratified = list(
    perf = c(
      .OLD_PERF, "shuffled", "species", "drug_label", "strat_label", "drug_or_class",
      "strat_value", "feature_type", "feature_subtype"
    ),
    top = c(
      .OLD_TOP, "shuffled", "species", "drug_label", "strat_label", "drug_or_class",
      "strat_value", "feature_type", "feature_subtype", "seed"
    )
  ),
  cross_drug = list(
    perf = c(
      .OLD_PERF_TEST, "shuffled", "species", "drug", "test_drug", "feature_type",
      "feature_subtype"
    ),
    top = NULL
  ),
  cross = list(
    perf = c(
      .OLD_PERF_TEST, "shuffled", "species", "drug_label", "drug_or_class", "strat_value",
      "strat_value_test", "feature_type", "feature_subtype"
    ),
    top = NULL
  ),
  loo = list(
    perf = c(
      .OLD_PERF_TEST, "shuffled", "species", "drug_label", "drug_or_class",
      "leaveout_strat_value", "feature_type", "feature_subtype"
    ),
    top = c(
      .OLD_TOP, "shuffled", "species", "drug_label", "drug_or_class", "leaveout_strat_value",
      "feature_type", "feature_subtype", "seed"
    )
  )
)

# Our column behind each old column whose name differs; the rest share their name.
.EXPORT_SOURCES <- c(
  drug_label = "unit", drug_or_class = "target", feature_subtype = "encoding",
  strat_value = "train_group", strat_value_test = "test_group",
  leaveout_strat_value = "test_group", drug = "train_group", test_drug = "test_group",
  num_obs_original_ml_input_tibble = "num_obs", num_obs_test_data = "n_test",
  res_prop_original_ml_input_tibble = "res_prop", res_prop_test_data = "res_prop_test"
)


# Each task's metadata from the tasks and matrix list; `species` comes with the results.
.exportMetadata <- function(fit) {
  m <- fit$built$matrices[match(fit$tasks$matrix_id, fit$built$matrices$matrix_id), ]
  grouping <- .MODES$grouping[match(m$mode_id, .MODES$mode_id)]

  tibble::tibble(
    task_id = fit$tasks$task_id,
    mode_id = m$mode_id,
    seed = fit$tasks$seed,
    shuffled = fit$tasks$labels == "shuffled",
    strat_label = ifelse(grouping %in% c("year", "country"), grouping, NA_character_),
    unit = m$unit,
    target = m$target,
    train_group = m$train_group,
    test_group = m$test_group,
    feature_type = m$feature_type,
    encoding = m$encoding
  )
}

# A mode's layout, chosen by its registry flags.
.exportLayout <- function(mode) {
  shape <- if (mode$grouping == "none") {
    "vanilla"
  } else if (mode$cross_test) {
    if (mode$grouping == "drug") "cross_drug" else "cross"
  } else if (mode$LOO) {
    "loo"
  } else {
    "stratified"
  }
  .EXPORT_LAYOUTS[[shape]]
}

# A mode's file names, from its registry flags as old amRml named them: all_, <grouping>_,
# cross_<grouping>_ or LOO_<grouping>_.
.exportFiles <- function(mode) {
  prefix <- if (mode$grouping == "none") {
    "all"
  } else if (mode$cross_test) {
    paste0("cross_", mode$grouping)
  } else if (mode$LOO) {
    paste0("LOO_", mode$grouping)
  } else {
    mode$grouping
  }
  c(perf = paste0(prefix, "_perf.parquet"), top = paste0(prefix, "_top_features.parquet"))
}

# One file's rows: results joined to their tasks' metadata, under the layout's column names.
.exportRows <- function(results, metadata, columns) {
  rows <- dplyr::inner_join(metadata, results, by = "task_id")
  sources <- ifelse(columns %in% names(.EXPORT_SOURCES), .EXPORT_SOURCES[columns], columns)
  stats::setNames(rows[sources], columns)
}

#' Export results in old amRml's merged format
#'
#' Writes a fit's performance and top features as the merged tables old amRml's
#' `mergeMLresults()` produced, one file per mode, with its file names and
#' columns: `all_perf.parquet` and `all_top_features.parquet` for vanilla,
#' `year_perf.parquet`, `cross_drug_perf.parquet`, `LOO_year_perf.parquet` and
#' so on. Metadata comes from the tasks rather than file names, and every file
#' has `shuffled` and an integer `seed`. This export is temporary: it bridges to
#' old amRml's format until amRml's own column names match it.
#'
#' @details
#' Left out: old amRml's `filename`, `train_prop`, `val_prop` and
#' `prop_vi_top_feats` bounds, and this package's added columns (`task_id`,
#' `auprc`), which stay in `results/`.
#'
#' @param fit An `amr_fit` from [fitModels()].
#' @param overwrite Replace files already in the output folder's `merged/`.
#'
#' @return The written paths, invisibly; the files are in `<out_dir>/merged/`.
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus `amrml_invalid_argument`: `fit`
#' is not from [fitModels()], or `merged/` already has files and `overwrite` is
#' `FALSE`.
#'
#' @examples
#' \dontrun{
#' exportMergedResults(fitModels(tasks))
#' }
#' @export
exportMergedResults <- function(fit, overwrite = FALSE) {
  call <- rlang::current_env()
  .checkArgClass(fit, "fit", "amr_fit", "fitModels()", call = call)
  .checkArgFlag(overwrite, "overwrite", call = call)
  merged <- file.path(fit$built$out_dir, "merged")
  if (length(list.files(merged)) && !overwrite) {
    .amrAbort(
      "invalid_argument",
      "The output folder already has merged results; set `overwrite = TRUE` to replace them.",
      observed = list(out_dir = fit$built$out_dir),
      call = call
    )
  }

  unlink(merged, recursive = TRUE)
  dir.create(merged)
  metadata <- .exportMetadata(fit)
  # The seed comes from the tasks, so top features get it too.
  performance <- fit$performance[names(fit$performance) != "seed"]
  written <- character()
  for (mode_id in intersect(.MODES$mode_id, metadata$mode_id)) {
    mode <- .MODES[.MODES$mode_id == mode_id, ]
    layout <- .exportLayout(mode)
    files <- file.path(merged, .exportFiles(mode))
    mode_metadata <- metadata[metadata$mode_id == mode_id, ]

    arrow::write_parquet(.exportRows(performance, mode_metadata, layout$perf), files[[1]])
    written <- c(written, files[[1]])
    if (!is.null(layout$top)) {
      arrow::write_parquet(.exportRows(fit$importance, mode_metadata, layout$top), files[[2]])
      written <- c(written, files[[2]])
    }
  }

  invisible(written)
}
