# Regression tests for control-flow bugs in
# createMLinputList()/runMLmodels().

test_that("runMLmodels() exits cleanly instead of crashing on no files", {
  tmp <- withr::local_tempdir()
  # No matrix parquet files are created, so createMLinputList() returns an
  # empty tibble. runMLmodels() must message and return early instead of
  # trying to use the empty result (which previously crashed with
  # "a character vector argument expected").
  result <- NULL
  expect_message(
    result <- runMLmodels(
      path = tmp, stratify_by = NULL, LOO = FALSE, cross_test = FALSE
    ),
    "No files found"
  )
  expect_null(result)
})

# `stratify_by = NULL` means "no stratification" when LOO is FALSE, and "leave
# one drug out" when LOO is TRUE. There is no stratify_by = "drug", so the drug
# case has to be spelled NULL. Every caller has to know that rule, and the ones
# that got it wrong are what the tests below cover.

test_that("createMLResultDir() maps every LOO/stratify_by pair to a real dir", {
  tmp <- withr::local_tempdir()
  matrix_dir <- function(loo, strat) {
    basename(createMLResultDir(
      tmp,
      stratify_by = strat, LOO = loo, cross_test = FALSE, MDR = FALSE
    )$matrix_path)
  }

  expect_equal(matrix_dir(FALSE, NULL), "matrix")
  expect_equal(matrix_dir(FALSE, "year"), "matrix_year")
  expect_equal(matrix_dir(FALSE, "country"), "matrix_country")
  expect_equal(matrix_dir(TRUE, NULL), "LOO_matrix_drug")
  expect_equal(matrix_dir(TRUE, "year"), "LOO_matrix_year")
  expect_equal(matrix_dir(TRUE, "country"), "LOO_matrix_country")
})

test_that("createMLResultDir() keeps LOO result dirs unsuffixed", {
  # Only matrix_path takes "_drug". Drug LOO results are documented as
  # LOO_ML_performance/ and LOO_ML_top_features/, unlike their year/country
  # siblings, so they must not become LOO_ML_drug_*.
  tmp <- withr::local_tempdir()
  paths <- createMLResultDir(
    tmp,
    stratify_by = NULL, LOO = TRUE, cross_test = FALSE, MDR = FALSE
  )
  expect_equal(basename(paths$ML_performance), "LOO_ML_performance")
  expect_equal(basename(paths$ML_top_features), "LOO_ML_top_features")
})

test_that("createMLinputList() accepts leave-one-drug-out and finds its files", {
  # This used to hard-error: the LOO check demanded stratify_by be "year" or
  # "country", which rejected the drug case outright.
  tmp <- withr::local_tempdir()
  paths <- createMLResultDir(
    tmp,
    stratify_by = NULL, LOO = TRUE, cross_test = FALSE, MDR = FALSE
  )
  for (drug in c("AMP", "CIP")) {
    file.create(file.path(
      paths$matrix_path,
      sprintf("Sfl_drug_leaveout_%s_gene_binary_sparse.parquet", drug)
    ))
  }

  out <- createMLinputList(
    tmp,
    stratify_by = NULL, LOO = TRUE, cross_test = FALSE, MDR = FALSE
  )

  expect_equal(nrow(out), 2)
  # One model per left-out drug, named from the file so the two do not collide.
  expect_setequal(
    out$output_prefix,
    c("Sfl_drug_leaveout_AMP_gene_binary", "Sfl_drug_leaveout_CIP_gene_binary")
  )
  expect_true(all(grepl("LOO_matrix_drug", out$matrix_path)))
})

test_that("createMLinputList() still rejects a bad stratify_by under LOO", {
  # NULL is valid (leave-one-drug-out), but a typo must still be caught.
  tmp <- withr::local_tempdir()
  expect_error(
    createMLinputList(tmp, LOO = TRUE, cross_test = FALSE, stratify_by = "bananas"),
    "stratify_by"
  )
})

test_that("leave-one-drug-out cross testing fails loudly, not silently", {
  # This mode needs a test set generateMLInputs() does not produce yet, so it
  # must error rather than quietly return zero rows.
  tmp <- withr::local_tempdir()
  paths <- createMLResultDir(
    tmp,
    stratify_by = NULL, LOO = TRUE, cross_test = TRUE, MDR = FALSE
  )
  file.create(file.path(
    paths$matrix_path,
    "Sfl_drug_leaveout_AMP_gene_binary_sparse.parquet"
  ))

  expect_error(
    createMLinputList(
      tmp,
      stratify_by = NULL, LOO = TRUE, cross_test = TRUE, MDR = FALSE
    ),
    "not supported yet"
  )
})
