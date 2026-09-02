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
