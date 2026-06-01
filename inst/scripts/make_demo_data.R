## Build the data() fixtures used by amRml example blocks:
##   demo_ml_tibble : 60-row, 80-feature ml_input_tibble (Sfl, AMP, genes, binary)
##   demo_fit       : a tuned LR workflow fitted on demo_ml_tibble
##
## Run from the package root:
##   Rscript inst/scripts/make_demo_data.R

suppressPackageStartupMessages({
  library(amRml)
  library(dplyr)
  library(rsample)
})

set.seed(42)
#Loading actual sfl data
fixture <- system.file("extdata", "Sfl_parquet.duckdb", package = "amRml")
out_dir <- file.path(tempdir(), "amRml_demo_build")
if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
dir.create(out_dir, recursive = TRUE)

#Generating the inputs that will be used as dataset to load for tests
generateMLInputs(
  parquet_duckdb_path = fixture,
  out_path            = out_dir,
  n_fold              = 5,
  split               = c(1, 0),
  min_n               = 25,
  verbosity           = "minimal"
)

target <- "genome_drug.resistant_phenotype"
full   <- loadMLInputTibble(file.path(
  out_dir, "matrix", "Sfl_drug_AMP_genes_binary_sparse.parquet"
))

#Subsetting to save on size and to run tests quickly
demo_ml_tibble <- full |>
  dplyr::group_by(.data[[target]]) |>
  dplyr::slice_sample(n = 30) |>
  dplyr::ungroup() |>
  dplyr::select(dplyr::all_of(c(
    "genome_id", target,
    head(setdiff(colnames(full), c("genome_id", target)), 80)
  )))

data_split <- splitMLInputTibble(demo_ml_tibble, split = c(1, 0), seed = 1)
train_data <- rsample::training(data_split)

wflow    <- buildWflow(buildLRModel(), buildRecipe(train_data, use_pca = FALSE))
grid     <- buildTuningGrid("LR", penalty_vec = 10^c(-3, -1), mix_vec = c(0, 0.5, 1))
tune_res <- tuneGrid(wflow, data_split, grid, n_fold = 2)
demo_fit <- fitBestModel(selectBestModel(tune_res, wflow, "mcc"), train_data)

usethis::use_data(demo_ml_tibble, demo_fit, overwrite = TRUE, compress = "xz")
