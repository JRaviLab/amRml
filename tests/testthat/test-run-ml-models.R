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

test_that("createMLinputList() rejects LOO w/o year/country (no cross-test)", {
  tmp <- withr::local_tempdir()
  expect_error(
    createMLinputList(tmp, LOO = TRUE, cross_test = FALSE, stratify_by = NULL),
    "stratify_by must be"
  )
})

test_that("createMLinputList() allows LOO+cross-test with stratify_by = NULL", {
  # LOO + cross_test with stratify_by = NULL is a distinct mode
  # (leave-one-drug-out cross-testing, "Case A" in the cross_test && LOO
  # branch) and is intentionally exempt from the year/country requirement
  # above.
  tmp <- withr::local_tempdir()
  paths <- createMLResultDir(
    tmp,
    stratify_by = NULL, LOO = TRUE, cross_test = TRUE, MDR = FALSE
  )

  # A single matrix filename that satisfies both the general filename parser
  # and the LOO-specific parser used inside Case A. File content is
  # irrelevant - createMLinputList() only inspects filenames.
  file.create(file.path(
    paths$matrix_path,
    "Csp_drug_leaveout_leaveout_genes_binary_sparse.parquet"
  ))

  out <- createMLinputList(
    tmp,
    stratify_by = NULL, LOO = TRUE, cross_test = TRUE, MDR = FALSE
  )

  # Before the return() fix, Case A built the right result but never
  # returned it, so execution fell through into the stratify_by != NULL
  # branch, which self-joins the single file against itself and always
  # filters it out (ref_file != test_file), silently returning 0 rows.
  expect_equal(nrow(out), 1)
  expect_true(grepl("_drug_leaveout_", out$output_prefix))
})
