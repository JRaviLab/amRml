# Unit tests for ML-prep helpers in prep_ml.R.

test_that(".getTargetVarName returns the correct target symbol", {
  bin <- make_binary_ml_tibble()
  mc <- make_multiclass_ml_tibble()
  expect_identical(
    as.character(.getTargetVarName(bin)),
    "genome_drug.resistant_phenotype"
  )
  expect_identical(
    as.character(.getTargetVarName(mc)),
    "resistant_classes"
  )
})

test_that(".getTargetVarName fails on invalid tibbles", {
  expect_error(.getTargetVarName(data.frame(x = 1)), "tibble")
  bad <- tibble::tibble(
    genome_id = "g1",
    genome_drug.resistant_phenotype = "R",
    resistant_classes = "A"
  )
  expect_error(.getTargetVarName(bad), "not both")
})

test_that("getNumFeat counts feature columns correctly", {
  bin <- make_binary_ml_tibble(n_genomes = 6, n_feats = 4)
  expect_equal(getNumFeat(bin), 4)

  mc <- make_multiclass_ml_tibble()
  # resistant_classes is a target, not a feature, despite its name prefix.
  expect_equal(getNumFeat(mc), 2)
})

test_that("shuffleLabels preserves labels and genome ordering", {
  bin <- make_binary_ml_tibble(n_genomes = 20, n_feats = 3, seed = 7)
  shuffled <- shuffleLabels(bin, seed = 42)

  expect_identical(shuffled$genome_id, bin$genome_id)
  # Features untouched
  expect_identical(shuffled$feat_1, bin$feat_1)
  # Same label multiset
  expect_identical(
    sort(as.character(shuffled$genome_drug.resistant_phenotype)),
    sort(as.character(bin$genome_drug.resistant_phenotype))
  )
  # Deterministic given the same seed
  expect_identical(shuffleLabels(bin, seed = 42), shuffled)
})

test_that("shuffleLabels rejects non-numeric seed", {
  bin <- make_binary_ml_tibble()
  expect_error(shuffleLabels(bin, seed = "x"), "numeric")
})

test_that("calculateMinSamples computes CV formula", {
  # Pure CV, balanced classes, 5 folds -> expect ceiling(5 / 0.5), i.e. 10.
  expect_identical(
    calculateMinSamples(n_fold = 5, split = c(1, 0), res_prop = 0.5),
    10
  )
  # smallest_n_obs_rs is a multiplier on top of the base requirement.
  expect_identical(
    calculateMinSamples(
      n_fold = 5, split = c(1, 0), res_prop = 0.5, smallest_n_obs_rs = 3
    ),
    30
  )
})

test_that("calculateMinSamples computes train/val/test formula", {
  # Split (0.6, 0.2) gives lowest partition 0.2; res_prop 0.25 gives
  # lowest class prop 0.25; expect ceiling(1 / (0.25 * 0.2)) = 20.
  expect_identical(
    calculateMinSamples(n_fold = NULL, split = c(0.6, 0.2), res_prop = 0.25),
    20
  )
})

test_that("calculateMinSamples returns Inf for degenerate classes", {
  expect_true(is.infinite(
    calculateMinSamples(n_fold = 5, split = c(1, 0), res_prop = 0)
  ))
  expect_true(is.infinite(
    calculateMinSamples(n_fold = 5, split = c(1, 0), res_prop = 1)
  ))
})

test_that("loadMLInputTibble pivots long parquet into wide tibble", {
  skip_if_not_installed("arrow")

  long <- tibble::tibble(
    genome_id = c("g1", "g1", "g2", "g2"),
    feature_id = c("f1", "f2", "f1", "f2"),
    value = c(1L, 0L, 0L, 1L),
    genome_drug.resistant_phenotype = c(
      "Resistant", "Resistant", "Susceptible", "Susceptible"
    )
  )
  parquet_path <- tempfile(fileext = ".parquet")
  on.exit(unlink(parquet_path))
  arrow::write_parquet(long, parquet_path)

  wide <- loadMLInputTibble(parquet_path)

  expect_true(tibble::is_tibble(wide))
  expect_setequal(
    colnames(wide),
    c("genome_id", "genome_drug.resistant_phenotype", "f1", "f2")
  )
  expect_identical(nrow(wide), 2L)
  expect_s3_class(wide$genome_drug.resistant_phenotype, "factor")
})

test_that("loadMLInputTibble errors on bad path / missing columns", {
  expect_error(loadMLInputTibble("/no/such/file.parquet"), "No file")

  skip_if_not_installed("arrow")
  bad <- tibble::tibble(
    genome_id = "g1",
    feature_id = "f1",
    value = 1L
    # no target column
  )
  bad_path <- tempfile(fileext = ".parquet")
  on.exit(unlink(bad_path))
  arrow::write_parquet(bad, bad_path)
  expect_error(loadMLInputTibble(bad_path), "target")
})
