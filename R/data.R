#' Demo ML input tibble
#'
#' Stratified subset (30 Resistant + 30 Susceptible) of the AMP-genes-binary
#' matrix from the bundled `Sfl_parquet.duckdb`, restricted to 80 feature
#' columns.
#'
#' @format A tibble with 60 rows and 82 columns: `genome_id`,
#'   `genome_drug.resistant_phenotype`, and 80 binary feature columns.
#' @source `inst/scripts/make_demo_data.R`.
"demo_ml_tibble"

#' Demo LR fit
#'
#' A tuned logistic-regression workflow fitted on [demo_ml_tibble].
#'
#' @format A fitted `workflow` object (output of [fitBestModel()]).
#' @source `inst/scripts/make_demo_data.R`.
"demo_fit"
