# Unit tests for parse_ml_filename() in merge_ml_results.R.
#
# parse_ml_filename() handles only unstratified result filenames. Year/country
# stratified filenames are aggregated by buildPerfPqYearCountry() instead.

test_that("parse_ml_filename parses an unstratified drug filename", {
  out <- parse_ml_filename("Csp_drug_AMX_genes_binary_42_top_features.parquet")

  expect_false(out$shuffled)
  expect_equal(out$species, "Csp")
  expect_equal(out$drug_label, "drug")
  expect_equal(out$drug_or_class, "AMX")
  expect_equal(out$feature_type, "genes")
  expect_equal(out$feature_subtype, "binary")
  expect_equal(out$seed, 42L)
})

test_that("parse_ml_filename parses an unstratified drug_class filename", {
  out <- parse_ml_filename(
    "Csp_drug_class_AMINOGLYCOSIDES_genes_binary_42_performance.parquet"
  )

  expect_equal(out$drug_label, "drug_class")
  expect_equal(out$drug_or_class, "AMINOGLYCOSIDES")
  expect_equal(out$seed, 42L)
})

test_that("parse_ml_filename detects a shuffled run", {
  out <- parse_ml_filename(
    "shuffled_Csp_drug_AMX_genes_binary_42_top_features.parquet"
  )

  expect_true(out$shuffled)
  expect_equal(out$drug_or_class, "AMX")
})

test_that("parse_ml_filename errors on a non-'drug' token after species", {
  expect_error(
    parse_ml_filename("Csp_widget_AMX_genes_binary_42_performance.parquet"),
    "expected 'drug' token"
  )
})

test_that("parse_ml_filename rejects a year-stratified filename", {
  # The strat label sits between "drug" and the drug value in stratified
  # filenames, e.g. "<species>_drug_year_<drug>_<year_range>_...". These are
  # handled by buildPerfPqYearCountry(), not this parser.
  expect_error(
    parse_ml_filename(
      "Csp_drug_year_AMX_2010-2015_genes_binary_year_42_performance.parquet"
    ),
    "does not support the 'year' stratified filename"
  )
})

test_that("parse_ml_filename rejects a country-stratified filename", {
  expect_error(
    parse_ml_filename(
      "Csp_drug_class_country_AMINOGLYCOSIDES_USA_genes_binary_country_7_top_features.parquet"
    ),
    "does not support the 'country' stratified filename"
  )
})

# ---------------------------------------------------------------------------
# buildPerfPqYearCountry() / buildTopFeatsPqYearCountry() filename parsing.
# Stratified filenames carry the strat label once, in the prefix:
#   <species>_<drug|drug_class>_<year|country>_<value>_<range>_<feat>_<sub>_<seed>_...
# ---------------------------------------------------------------------------

test_that("buildPerfPqYearCountry parses year-stratified drug filenames", {
  dir <- withr::local_tempdir()
  arrow::write_parquet(
    tibble::tibble(metric = "roc_auc", estimate = 0.9),
    file.path(dir, "Csp_drug_year_AMX_2010-2015_genes_binary_42_performance.parquet")
  )
  arrow::write_parquet(
    tibble::tibble(metric = "roc_auc", estimate = 0.8),
    file.path(dir, "Csp_drug_year_CIP_2016-2020_genes_binary_42_performance.parquet")
  )

  out <- buildPerfPqYearCountry(dir)

  expect_setequal(out$drug_or_class, c("AMX", "CIP"))
  expect_setequal(out$strat_value, c("2010-2015", "2016-2020"))
  expect_true(all(out$species == "Csp"))
  expect_true(all(out$drug_label == "drug"))
  expect_true(all(out$strat_label == "year"))
  expect_true(all(out$feature_type == "genes"))
  expect_true(all(out$feature_subtype == "binary"))
  expect_false("seed_from_name" %in% names(out))
  expect_true(file.exists(file.path(dir, "year_perf.parquet")))
})

test_that("buildPerfPqYearCountry parses country-stratified drug_class filenames", {
  dir <- withr::local_tempdir()
  arrow::write_parquet(
    tibble::tibble(metric = "roc_auc", estimate = 0.7),
    file.path(
      dir,
      "Csp_drug_class_country_AMINOGLYCOSIDES_USA_genes_binary_7_performance.parquet"
    )
  )

  out <- buildPerfPqYearCountry(dir)

  expect_equal(out$drug_label, "drug_class")
  expect_equal(out$strat_label, "country")
  expect_equal(out$drug_or_class, "AMINOGLYCOSIDES")
  expect_equal(out$strat_value, "USA")
  expect_true(file.exists(file.path(dir, "country_perf.parquet")))
})

test_that("buildTopFeatsPqYearCountry parses stratified filenames", {
  dir <- withr::local_tempdir()
  arrow::write_parquet(
    tibble::tibble(Feature = "gyrA", Importance = "0.42"),
    file.path(dir, "Csp_drug_year_AMX_2010-2015_genes_binary_42_top_features.parquet")
  )

  buildTopFeatsPqYearCountry(dir)

  res <- arrow::read_parquet(file.path(dir, "year_top_features.parquet"))
  expect_equal(res$species, "Csp")
  expect_equal(res$drug_label, "drug")
  expect_equal(res$strat_label, "year")
  expect_equal(res$drug_or_class, "AMX")
  expect_equal(res$strat_value, "2010-2015")
  expect_type(res$Importance, "double")
})
