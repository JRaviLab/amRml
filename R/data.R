# FLAG for Abhirupa: data/demo_fit.rda and data/demo_ml_tibble.rda were
# deleted in b8f45b1 ("Add example data and remove old, unused data"), but
# the roxygen docs below (and @examples in plot_ml.R, run_ml_pipeline.R,
# ife_ml.R, core_ml.R that call data(demo_ml_tibble)) still reference them.
# Confirmed via `exists("demo_ml_tibble")` after devtools::load_all() -> FALSE,
# the object genuinely isn't in the package anymore, so devtools::document()
# fatally errors trying to build .Rd files for data that doesn't exist ('demo_ml_tibble'
# is not an exported object from 'namespace:amRml'). Commented out both blocks
# below as the minimal unblock so NAMESPACE regeneration can proceed for the
# rest of the package. The real question -- restore the .rda files via
# inst/scripts/make_demo_data.R, or retire these two blocks + the four
# @examples blocks elsewhere that still call data(demo_ml_tibble) -- is still
# open and needs her call; nothing below was deleted, just disabled.

# #' Demo ML input tibble
# #'
# #' Stratified subset (30 Resistant + 30 Susceptible) of the AMP-genes-binary
# #' matrix from the bundled `Sfl_parquet.duckdb`, restricted to 80 feature
# #' columns.
# #'
# #' @format A tibble with 60 rows and 82 columns: `genome_id`,
# #'   `genome_drug.resistant_phenotype`, and 80 binary feature columns.
# #' @source `inst/scripts/make_demo_data.R`.
# #' @examples
# #' data(demo_ml_tibble)
# #' dim(demo_ml_tibble)
# "demo_ml_tibble"
#
# #' Demo LR fit
# #'
# #' A tuned logistic-regression workflow fitted on `demo_ml_tibble`.
# #'
# #' @format A fitted `workflow` object (output of [fitBestModel()]).
# #' @source `inst/scripts/make_demo_data.R`.
# #' @examples
# #' data(demo_fit)
# #' class(demo_fit)
# "demo_fit"
