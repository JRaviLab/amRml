# Unit tests for parse_ml_filename() in merge_ml_results.R.

test_that("parse_ml_filename parses an unstratified drug filename", {
  out <- parse_ml_filename("Csp_drug_AMX_genes_binary_42_top_features.tsv")

  expect_false(out$shuffled)
  expect_equal(out$species, "Csp")
  expect_equal(out$drug_label, "drug")
  expect_equal(out$drug_or_class, "AMX")
  expect_true(is.na(out$strat_label))
  expect_true(is.na(out$strat_value))
  expect_equal(out$feature_type, "genes")
  expect_equal(out$feature_subtype, "binary")
  expect_equal(out$seed, 42L)
})

test_that("parse_ml_filename parses an unstratified drug_class filename", {
  out <- parse_ml_filename("Csp_drug_class_AMINOGLYCOSIDES_genes_binary_42_performance.tsv")

  expect_equal(out$drug_label, "drug_class")
  expect_equal(out$drug_or_class, "AMINOGLYCOSIDES")
  expect_true(is.na(out$strat_label))
  expect_equal(out$seed, 42L)
})

test_that("parse_ml_filename detects a shuffled run", {
  out <- parse_ml_filename("shuffled_Csp_drug_AMX_genes_binary_42_top_features.tsv")

  expect_true(out$shuffled)
  expect_equal(out$drug_or_class, "AMX")
})

test_that("parse_ml_filename parses a year-stratified drug filename", {
  # The strat label sits between "drug" and the drug value in the actual
  # filenames written by the matrix-generation code, e.g.
  # "<species>_drug_year_<drug>_<year_range>_...".
  out <- parse_ml_filename(
    "Csp_drug_year_AMX_2010-2015_genes_binary_year_42_performance.tsv"
  )

  expect_equal(out$species, "Csp")
  expect_equal(out$drug_label, "drug")
  expect_equal(out$drug_or_class, "AMX")
  expect_equal(out$strat_label, "year")
  expect_equal(out$strat_value, "2010-2015")
  expect_equal(out$feature_type, "genes")
  expect_equal(out$feature_subtype, "binary")
  expect_equal(out$seed, 42L)
})

test_that("parse_ml_filename parses a country-stratified drug_class filename", {
  out <- parse_ml_filename(
    "Csp_drug_class_country_AMINOGLYCOSIDES_USA_genes_binary_country_7_top_features.tsv"
  )

  expect_equal(out$drug_label, "drug_class")
  expect_equal(out$drug_or_class, "AMINOGLYCOSIDES")
  expect_equal(out$strat_label, "country")
  expect_equal(out$strat_value, "USA")
  expect_equal(out$seed, 7L)
})
