# Regression tests for the DuckDB connection cleanup bug in
# .parquet2Matrix()/.parquet2MDRMatrix(): the outer on.exit() handler checked
# a variable named `con` instead of `con0`, so if setup failed before the
# per-group loop ever created `con`, the real error was masked by a second,
# confusing "object 'con' not found" error thrown while evaluating the
# on.exit handler itself - on top of never closing con0.

make_feature_only_parquet_dir <- function(.local_envir = parent.frame()) {
  # A directory with a valid feature-count parquet (so .getFeatureTypes()
  # succeeds and con0 gets created) but no metadata.parquet, so
  # .register_parquet_views(con0, ...) fails naturally right after con0 is
  # created - exactly the failure window the bug affects.
  # Cleanup is tied to the caller's frame (.local_envir), not this helper's
  # own frame, so the directory survives until the calling test finishes.
  dir <- withr::local_tempdir(.local_envir = .local_envir)
  arrow::write_parquet(
    tibble::tibble(
      genome_id = character(),
      gene = character(),
      value = double()
    ),
    file.path(dir, "gene_count.parquet")
  )
  dir
}

test_that(".parquet2Matrix() surfaces the real setup error, not a masked one", {
  parquet_dir <- make_feature_only_parquet_dir()
  out_dir <- withr::local_tempdir()

  expect_error(
    amRml:::.parquet2Matrix(
      parquet_dir = parquet_dir,
      path = out_dir,
      n_fold = 5,
      split = c(0.8, 0.2)
    ),
    "metadata parquet not found"
  )
})

test_that(".parquet2MDRMatrix() surfaces the real setup error, not masked", {
  parquet_dir <- make_feature_only_parquet_dir()
  out_dir <- withr::local_tempdir()

  expect_error(
    amRml:::.parquet2MDRMatrix(
      parquet_dir = parquet_dir,
      path = out_dir,
      min_n = 25
    ),
    "metadata parquet not found"
  )
})
