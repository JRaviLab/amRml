# Unit tests for internal argument validators in arg_check_ml.R.

test_that(".checkArgPath rejects non-character and missing files", {
  expect_error(amRml:::.checkArgPath(42), "character")
  expect_error(amRml:::.checkArgPath("/does/not/exist.parquet"), "No file")
})

test_that(".checkArgPath accepts an existing path", {
  tmp <- tempfile()
  file.create(tmp)
  on.exit(unlink(tmp))
  expect_silent(amRml:::.checkArgPath(tmp))
})

test_that(".checkArgTibble enforces tibble type and non-empty dimensions", {
  expect_error(amRml:::.checkArgTibble(data.frame(x = 1)), "tibble")
  expect_error(
    amRml:::.checkArgTibble(tibble::tibble()),
    "0 rows"
  )
})

test_that(".checkArgTibble(ml=TRUE) enforces exactly one target column", {
  both <- tibble::tibble(
    genome_id = "g1",
    genome_drug.resistant_phenotype = "Resistant",
    resistant_classes = "A"
  )
  neither <- tibble::tibble(genome_id = "g1", feat = 1L)
  ok <- tibble::tibble(
    genome_id = "g1",
    genome_drug.resistant_phenotype = "Resistant"
  )

  expect_error(amRml:::.checkArgTibble(both, ml = TRUE), "not both")
  expect_error(amRml:::.checkArgTibble(neither, ml = TRUE), "not both")
  expect_silent(amRml:::.checkArgTibble(ok, ml = TRUE))
})

test_that(".checkArgSeed requires numeric", {
  expect_error(amRml:::.checkArgSeed("x"), "numeric")
  expect_silent(amRml:::.checkArgSeed(42))
})

test_that(".checkArgNFold requires whole number > 1", {
  expect_error(amRml:::.checkArgNFold("x"), "numeric")
  expect_error(amRml:::.checkArgNFold(1), "whole")
  expect_error(amRml:::.checkArgNFold(2.5), "whole")
  expect_silent(amRml:::.checkArgNFold(5))
})

test_that(".checkArgSplit allows pure-CV and train/val/test, rejects others", {
  expect_silent(amRml:::.checkArgSplit(c(1, 0)))
  expect_silent(amRml:::.checkArgSplit(c(0.6, 0.2)))

  expect_error(amRml:::.checkArgSplit(c(0.8)), "length 2")
  expect_error(amRml:::.checkArgSplit(c(NA_real_, 0.2)), "NA")
  expect_error(amRml:::.checkArgSplit(c(-0.1, 0.5)), "\\[0, 1\\]")
  expect_error(amRml:::.checkArgSplit(c(0.8, 0.3)), "Unsupported|<= 1")
  # train>0, val=0, test>0 is neither pure-CV nor train/val/test → rejected
  expect_error(amRml:::.checkArgSplit(c(0.7, 0)), "Unsupported")
})

test_that(".checkArgSplitAndMode enforces mode/n_fold consistency", {
  expect_identical(
    amRml:::.checkArgSplitAndMode(c(1, 0), n_fold = 5),
    "cv"
  )
  expect_identical(
    amRml:::.checkArgSplitAndMode(c(0.6, 0.2), n_fold = NULL),
    "splits"
  )
  expect_error(
    amRml:::.checkArgSplitAndMode(c(1, 0), n_fold = NULL),
    "n_fold"
  )
  expect_error(
    amRml:::.checkArgSplitAndMode(c(0.6, 0.2), n_fold = 5),
    "n_fold = NULL"
  )
})

test_that(".checkArgResProp enforces [0, 1] numeric", {
  expect_error(amRml:::.checkArgResProp("x"), "numeric")
  expect_error(amRml:::.checkArgResProp(-0.01), "\\[0, 1\\]")
  expect_error(amRml:::.checkArgResProp(1.01), "\\[0, 1\\]")
  expect_silent(amRml:::.checkArgResProp(0.3))
})

test_that(".checkArgSmallestNObsRS requires positive whole number", {
  expect_error(amRml:::.checkArgSmallestNObsRS("x"), "numeric")
  expect_error(amRml:::.checkArgSmallestNObsRS(0), "whole")
  expect_error(amRml:::.checkArgSmallestNObsRS(2.5), "whole")
  expect_silent(amRml:::.checkArgSmallestNObsRS(3))
})

test_that(".checkArgModel accepts only LR", {
  expect_silent(amRml:::.checkArgModel("LR"))
})

test_that(".checkArgPenaltyVec and .checkArgMixVec enforce ranges", {
  expect_error(amRml:::.checkArgPenaltyVec("x"), "numeric")
  expect_error(amRml:::.checkArgPenaltyVec(c(0.1, -1)), "less than zero")
  expect_silent(amRml:::.checkArgPenaltyVec(c(1e-3, 1, 10)))

  expect_error(amRml:::.checkArgMixVec("x"), "numeric")
  expect_error(amRml:::.checkArgMixVec(c(0, 1.1)), "\\[0, 1\\]")
  expect_silent(amRml:::.checkArgMixVec(c(0, 0.5, 1)))
})

test_that(".checkArgNFeat requires numeric >= 1 and not NA", {
  expect_error(amRml:::.checkArgNFeat(NA), "must be specified")
  expect_error(amRml:::.checkArgNFeat("x"), "numeric")
  expect_error(amRml:::.checkArgNFeat(0), "at least 1")
  expect_silent(amRml:::.checkArgNFeat(5))
})

test_that(".checkArgSelectBestMetric accepts only known metrics", {
  for (m in c("f_meas", "pr_auc", "mcc", "bal_accuracy")) {
    expect_silent(amRml:::.checkArgSelectBestMetric(m))
  }
  expect_error(amRml:::.checkArgSelectBestMetric("roc_auc"), "one of")
})

test_that(".checkArgPropVITopFeats enforces length-2 ordered [0,1]", {
  expect_error(amRml:::.checkArgPropVITopFeats(0.5), "length")
  expect_error(amRml:::.checkArgPropVITopFeats(c(-0.1, 0.5)), "\\[0, 1\\]")
  expect_error(amRml:::.checkArgPropVITopFeats(c(0.3, 0.3)), "less than")
  expect_error(amRml:::.checkArgPropVITopFeats(c(0.5, 0.2)), "less than")
  expect_silent(amRml:::.checkArgPropVITopFeats(c(0.1, 0.2)))
})

test_that(".checkArgPercentRemovalVec enforces ascending values in [0, 100]", {
  expect_error(amRml:::.checkArgPercentRemovalVec(c("a", "b")), "numeric")
  expect_error(amRml:::.checkArgPercentRemovalVec(c(-1, 10)), "\\[0, 100\\]")
  expect_error(amRml:::.checkArgPercentRemovalVec(c(10, 200)), "\\[0, 100\\]")
  expect_error(amRml:::.checkArgPercentRemovalVec(c(10, 10, 20)), "ascending")
  expect_error(amRml:::.checkArgPercentRemovalVec(c(30, 20)), "ascending")
  expect_silent(amRml:::.checkArgPercentRemovalVec(c(10, 20, 30)))
})

test_that("logical-only validators reject non-logicals", {
  logical_checkers <- list(
    amRml:::.checkArgUsePCA,
    amRml:::.checkArgMultiClass,
    amRml:::.checkArgShuffleLabels,
    amRml:::.checkArgReturnTuneRes,
    amRml:::.checkArgReturnFit,
    amRml:::.checkArgReturnPred,
    amRml:::.checkArgVerbose,
    amRml:::.checkArgByVI,
    amRml:::.checkArgByNum,
    amRml:::.checkArgReturnFeats
  )
  for (chk in logical_checkers) {
    expect_error(chk("TRUE"), "logical")
    expect_silent(chk(TRUE))
    expect_silent(chk(FALSE))
  }
})

test_that("plot-axis validators enforce allowed strings", {
  expect_error(amRml:::.checkArgXDefaultEval(1), "character")
  expect_error(amRml:::.checkArgXDefaultEval("nope"), "train_prop")
  expect_silent(amRml:::.checkArgXDefaultEval("train_prop"))

  expect_error(amRml:::.checkArgYDefaultEval(1), "character")
  expect_error(amRml:::.checkArgYDefaultEval("nope"), "avg_f1_score")
  expect_silent(amRml:::.checkArgYDefaultEval("avg_nmcc"))

  expect_error(amRml:::.checkArgXYLabs(1, "y"), "character")
  expect_silent(amRml:::.checkArgXYLabs("x", "y"))
})
