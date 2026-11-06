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
