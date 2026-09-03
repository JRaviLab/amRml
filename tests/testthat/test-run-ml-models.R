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

# LOO comes in three kinds, and each names the stratum it holds out:
#   stratify_by = "year"     -> LOO_matrix_year/
#   stratify_by = "country"  -> LOO_matrix_country/
#   stratify_by = "drug"     -> LOO_matrix_drug/
# NULL is not a LOO mode: it means "no stratification", which is only valid
# when LOO is FALSE.

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
  expect_equal(matrix_dir(TRUE, "drug"), "LOO_matrix_drug")
  expect_equal(matrix_dir(TRUE, "year"), "LOO_matrix_year")
  expect_equal(matrix_dir(TRUE, "country"), "LOO_matrix_country")
})

test_that("createMLResultDir() suffixes drug LOO result dirs like its siblings", {
  # Every stratum suffixes both the matrix path and the result dirs, so drug
  # follows year/country rather than being a special case.
  tmp <- withr::local_tempdir()
  paths <- createMLResultDir(
    tmp,
    stratify_by = "drug", LOO = TRUE, cross_test = FALSE, MDR = FALSE
  )
  expect_equal(basename(paths$ML_performance), "LOO_ML_drug_performance")
  expect_equal(basename(paths$ML_top_features), "LOO_ML_drug_top_features")
})

test_that("createMLinputList() rejects LOO without a stratum", {
  # NULL means "no stratification", which cannot be a leave-one-out mode:
  # every LOO matrix directory carries a suffix.
  tmp <- withr::local_tempdir()
  expect_error(
    createMLinputList(tmp, LOO = TRUE, cross_test = FALSE, stratify_by = NULL),
    "must be 'year', 'country', or 'drug'"
  )
})

test_that("createMLinputList() accepts leave-one-drug-out and finds its files", {
  # This used to hard-error: the LOO check demanded stratify_by be "year" or
  # "country", so there was no way to ask for the drug case at all.
  tmp <- withr::local_tempdir()
  paths <- createMLResultDir(
    tmp,
    stratify_by = "drug", LOO = TRUE, cross_test = FALSE, MDR = FALSE
  )
  for (drug in c("AMP", "CIP")) {
    file.create(file.path(
      paths$matrix_path,
      sprintf("Sfl_drug_leaveout_%s_gene_binary_sparse.parquet", drug)
    ))
  }

  out <- createMLinputList(
    tmp,
    stratify_by = "drug", LOO = TRUE, cross_test = FALSE, MDR = FALSE
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
  # A typo must still be caught, not silently treated as a stratum.
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
    stratify_by = "drug", LOO = TRUE, cross_test = TRUE, MDR = FALSE
  )
  file.create(file.path(
    paths$matrix_path,
    "Sfl_drug_leaveout_AMP_gene_binary_sparse.parquet"
  ))

  expect_error(
    createMLinputList(
      tmp,
      stratify_by = "drug", LOO = TRUE, cross_test = TRUE, MDR = FALSE
    ),
    "not supported yet"
  )
})
