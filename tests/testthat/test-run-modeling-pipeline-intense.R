# Unit test for runModelingPipelineIntense()'s control flow. Mocks out the
# heavy steps (matrix generation, model training) so this stays fast and only
# checks that the wrapper actually calls them in order, rather than exercising
# real matrix generation / model fitting (already covered elsewhere).

test_that("runModelingPipelineIntense() generates matrices before training", {
  tmp_root <- withr::local_tempdir()
  parquet_dir <- file.path(tmp_root, "Test_species")
  dir.create(parquet_dir)

  generate_calls <- list()
  local_mocked_bindings(
    generateMLInputs = function(...) {
      generate_calls[[length(generate_calls) + 1]] <<- list(...)
      invisible(NULL)
    },
    runMLmodels = function(...) invisible(NULL),
    runMDRmodels = function(...) invisible(NULL)
  )

  runModelingPipelineIntense(
    parquet_dir = parquet_dir,
    n_seeds = 1,
    verbose = FALSE
  )

  # The matrix-generation step must actually run (regression test: this step
  # was previously commented out, so every later step silently had nothing
  # to work with).
  expect_length(generate_calls, 1)
  expect_equal(generate_calls[[1]]$parquet_dir, normalizePath(parquet_dir))
  expect_equal(generate_calls[[1]]$out_path, normalizePath(tmp_root))
})
