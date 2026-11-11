testthat::test_that("feature rescoring retains paired and orphan dyads", {
  testthat::skip_if_not_installed("duckdb")

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(con, "dyad_feature", data.frame(
    dyad = c("P1|g1", "P2|orphan", "orphan|g2"),
    feature = c("protein:P1", "protein:P2", "gene:g2")
  ))

  top <- tibble::tibble(
    species = "Tst",
    drug_label = "drug",
    drug_or_class = "AMP",
    feature_type = c("protein", "protein", "gene"),
    feature_subtype = "binary",
    variable = c("P1", "P2", "g2"),
    sign = "NEG",
    median_rank_score = c(0.9, 0.8, 0.7),
    median_contribution = c(0.5, 0.4, 0.3),
    subtype_csv = "binary"
  )

  result <- summariseDyads(top, con)

  testthat::expect_setequal(
    result$dyad,
    c("P1|g1", "P2|orphan", "orphan|g2")
  )
  testthat::expect_equal(result$n_variables, rep(1L, 3L))
  testthat::expect_equal(result$frequency, rep(1L, 3L))
  testthat::expect_equal(
    result$dyad_median_rank_score,
    c(0.9, 0.8, 0.7)
  )
})


test_that("shuffled comparisons pair independently tuned fits and enforce direction", {
  path <- withr::local_tempfile(fileext = ".parquet")
  perf <- tibble::tibble(
    species = "Tst", drug_label = "drug", drug_or_class = "AMP",
    feature_type = "gene", feature_subtype = "binary", model = "LR",
    seed = rep(1:3, each = 2), shuffled = rep(c(FALSE, TRUE), 3),
    fit_penalty = rep(c(0.01, 0.1), 3), fit_mixture = rep(c(1, 0), 3),
    mcc = c(0.4, 0.8, -0.6, 0.1, 0.7, 0.1)
  )
  arrow::write_parquet(perf, path)
  positive <- filterOptimalModel(path)
  expect_equal(positive$seed, 3L)
  expect_equal(positive$shuffled_MCC, 0.1)
  expect_equal(positive$fit_penalty, 0.01)
  negative <- filterOptimalModel(path, direction = "negative", MCC_threshold = 0.5)
  expect_equal(negative$seed, 2L)
  expect_equal(negative$MCC_diff, -0.7)
  expect_equal(nrow(filterOptimalModel(path, direction = "negative", MCC_threshold = 0.7)), 0L)
  perf$mcc[6] <- NA_real_
  arrow::write_parquet(perf, path)
  expect_equal(nrow(filterOptimalModel(path)), 0L)
})

test_that("inverse cross-test reliability requires repeated negative evidence within scope", {
  path <- withr::local_tempfile(fileext = ".parquet")
  perf <- tidyr::expand_grid(seed = 1:3, strat_value_test = c("inverse", "positive", "mixed"),
                             shuffled = c(FALSE, TRUE)) |>
    dplyr::mutate(species = "Tst", drug_label = "drug", drug_or_class = "AMP",
                  strat_value = "train", feature_type = "gene", feature_subtype = "binary",
                  mcc = dplyr::case_when(
                    shuffled ~ 0.1 * (seed - 2),
                    strat_value_test == "inverse" ~ -0.6,
                    strat_value_test == "positive" ~ 0.05,
                    seed == 1 ~ -0.6,
                    TRUE ~ 0.2
                  ))
  arrow::write_parquet(perf, path)
  inverse <- assessModelReliability(path, direction = "negative")$by_representation
  expect_equal(inverse$strat_value_test[inverse$supported], "inverse")
  expect_equal(inverse$n_supported_seeds[inverse$strat_value_test == "mixed"], 1L)
  expect_true(all(inverse$direction == "negative"))
  positive <- assessModelReliability(path)$by_representation
  expect_equal(positive$strat_value_test[positive$supported], "mixed")
})
