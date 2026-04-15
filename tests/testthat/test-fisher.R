# Unit tests for Fisher / Benjamini-Hochberg helpers in run_Fisher_tests.R.

# Build a tiny wide binary matrix where feat_good perfectly predicts the
# phenotype and feat_bad is independent of it.
make_fisher_fixture <- function() {
  tibble::tibble(
    genome_id = paste0("g", 1:8),
    genome_drug.resistant_phenotype = c(
      "Resistant", "Resistant", "Resistant", "Resistant",
      "Susceptible", "Susceptible", "Susceptible", "Susceptible"
    ),
    feat_good = c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L),
    feat_bad = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L)
  )
}

test_that("encodePhenotype encodes Susceptible=0 / Resistant=1", {
  df <- make_fisher_fixture()
  out <- encodePhenotype(df)

  expect_named(out, c("df", "target"))
  expect_type(out$target, "integer")
  expect_identical(out$target, df$genome_drug.resistant_phenotype |>
    (\(x) ifelse(x == "Resistant", 1L, 0L))())
  expect_identical(
    out$df$genome_drug.resistant_phenotype,
    out$target
  )
})

test_that("encodePhenotype filters unknown labels", {
  df <- make_fisher_fixture()
  df$genome_drug.resistant_phenotype[1] <- "Intermediate"
  out <- encodePhenotype(df)
  expect_identical(nrow(out$df), 7L)
  expect_false(any(is.na(out$target)))
})

test_that("encodePhenotype errors if target column is missing", {
  bad <- tibble::tibble(genome_id = "g1", feat_good = 1L)
  expect_error(encodePhenotype(bad))
})

test_that("runFisherTests orders genes by ascending p-value", {
  df <- make_fisher_fixture()
  enc <- encodePhenotype(df)
  out <- runFisherTests(enc$df, enc$target)

  expect_setequal(colnames(out), c("gene", "p_value", "alternative"))
  # feat_good perfectly separates classes, so its p-value must be
  # strictly smaller than feat_bad's.
  expect_identical(out$gene[1], "feat_good")
  expect_lt(
    out$p_value[out$gene == "feat_good"],
    out$p_value[out$gene == "feat_bad"]
  )
  expect_true(all(out$alternative == "two.sided"))
})

test_that("applyBenjaminiHochberg adds BH columns and flags significance", {
  df <- make_fisher_fixture()
  enc <- encodePhenotype(df)
  fisher <- runFisherTests(enc$df, enc$target)
  bh <- applyBenjaminiHochberg(fisher, Q = 0.05)

  expect_true(all(
    c("adj_p_value", "sig_after_bh", "Q") %in% colnames(bh)
  ))
  expect_true(all(bh$Q == 0.05))
  expect_type(bh$sig_after_bh, "logical")
  # adj_p_value column should be positioned right after p_value.
  cn <- colnames(bh)
  expect_identical(cn[which(cn == "p_value") + 1], "adj_p_value")
  expect_identical(cn[which(cn == "adj_p_value") + 1], "sig_after_bh")
})

test_that("applyBenjaminiHochberg requires p_value column", {
  expect_error(applyBenjaminiHochberg(tibble::tibble(x = 1)), "p_value")
})

test_that("computeFeatureFreq adds class-conditional frequencies in [0, 1]", {
  df <- make_fisher_fixture()
  enc <- encodePhenotype(df)
  fisher <- runFisherTests(enc$df, enc$target) |>
    applyBenjaminiHochberg(Q = 0.05)
  out <- computeFeatureFreq(enc$df, fisher, enc$target)

  expect_true(all(
    c("freq_susceptible_gene_pres", "freq_resistant_gene_pres") %in%
      colnames(out)
  ))
  expect_true(all(out$freq_susceptible_gene_pres >= 0 &
    out$freq_susceptible_gene_pres <= 1))
  expect_true(all(out$freq_resistant_gene_pres >= 0 &
    out$freq_resistant_gene_pres <= 1))

  # feat_good is present in 100% of Resistant and 0% of Susceptible genomes.
  good <- out[out$gene == "feat_good", ]
  expect_identical(good$freq_resistant_gene_pres, 1)
  expect_identical(good$freq_susceptible_gene_pres, 0)
})
