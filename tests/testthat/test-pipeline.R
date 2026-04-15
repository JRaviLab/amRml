# End-to-end tests for runMLPipeline(). Uses a small tuning grid to keep runtime low while still exercising the full
# split -> tune -> fit -> predict -> extract path.
#
# Tests run in pure-CV mode (split = c(1, 0)) because the classical
# train/val/test path in runMLPipeline currently zeroes n_fold before calling
# tuneGrid(), which still requires it.

skip_pipeline_if_missing <- function() {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("tune")
  skip_if_not_installed("workflows")
  skip_if_not_installed("yardstick")
  skip_if_not_installed("vip")
}

minimal_grid_args <- list(
  penalty_vec = c(0.01, 0.1),
  mix_vec = c(0, 1)
)

run_pipeline <- function(fx, ...) {
  suppressMessages(suppressWarnings(
    runMLPipeline(
      fx,
      model = "LR",
      split = c(1, 0),
      n_fold = 2,
      n_top_feats = 2,
      prop_vi_top_feats = NA,
      penalty_vec = minimal_grid_args$penalty_vec,
      mix_vec = minimal_grid_args$mix_vec,
      verbose = FALSE,
      ...
    )
  ))
}

test_that("runMLPipeline returns expected structure in CV mode", {
  skip_pipeline_if_missing()
  fx <- make_pipeline_fixture(n_per_class = 25)

  res <- run_pipeline(fx, seed = 42)

  expect_named(res, c("performance_tibble", "top_feat_tibble"),
    ignore.order = TRUE
  )
  expect_s3_class(res$performance_tibble, "tbl_df")
  expect_identical(nrow(res$performance_tibble), 1L)
  expect_true(all(
    c(
      "num_obs", "res_prop", "model", "train_prop", "val_prop",
      "n_fold", "nmcc", "run_time_sec"
    ) %in%
      colnames(res$performance_tibble)
  ))
  expect_identical(res$performance_tibble$model, "LR")
  expect_identical(res$performance_tibble$train_prop, 1)
  expect_identical(res$performance_tibble$val_prop, 0)
  expect_identical(res$performance_tibble$n_fold, 2)
  expect_gte(res$performance_tibble$nmcc, 0)
  expect_lte(res$performance_tibble$nmcc, 1)

  expect_s3_class(res$top_feat_tibble, "tbl_df")
  # With n_top_feats = 2 the pipeline pads up to 2 features even if vip
  # assigned importance to fewer.
  expect_identical(nrow(res$top_feat_tibble), 2L)
})

test_that("runMLPipeline recovers the informative signal", {
  skip_pipeline_if_missing()
  fx <- make_pipeline_fixture(n_per_class = 30, seed = 99)
  res <- run_pipeline(fx, seed = 99)

  expect_gt(res$performance_tibble$nmcc, 0.6)
})

test_that("runMLPipeline returns fit/pred/tune_res when requested", {
  skip_pipeline_if_missing()
  fx <- make_pipeline_fixture(n_per_class = 25)
  res <- run_pipeline(
    fx,
    seed = 1,
    return_fit = TRUE,
    return_pred = TRUE,
    return_tune_res = TRUE
  )

  expect_true(all(c("fit", "pred", "tune_res") %in% names(res)))
  expect_s3_class(res$fit, "workflow")
  expect_true(".pred_class" %in% colnames(res$pred))
  expect_s3_class(res$tune_res, "tune_results")
})

test_that("runMLPipeline shuffle_labels produces a valid baseline result", {
  skip_pipeline_if_missing()
  fx <- make_pipeline_fixture(n_per_class = 25)
  baseline <- run_pipeline(fx, seed = 1, shuffle_labels = TRUE)

  expect_identical(nrow(baseline$performance_tibble), 1L)
  expect_false(is.na(baseline$performance_tibble$nmcc))
})

test_that("runMLPipeline rejects non-LR models for multi-class data", {
  skip_pipeline_if_missing()
  mc <- tibble::tibble(
    genome_id = paste0("g", 1:12),
    resistant_classes = factor(rep(c("A", "B", "C"), each = 4)),
    feat_1 = rep(c(1L, 0L), 6),
    feat_2 = rep(c(0L, 1L), 6)
  )
  expect_error(
    suppressMessages(runMLPipeline(
      mc,
      model = "RF", split = c(1, 0), n_fold = 2, verbose = FALSE
    )),
    "[Mm]ulti-class"
  )
})

test_that("removeTopFeats drops requested feature columns", {
  fx <- make_pipeline_fixture(n_per_class = 10)
  top <- tibble::tibble(
    Variable = c("feat_inf_1", "feat_noise_1"),
    Importance = c(0.9, 0.1),
    Sign = c("POS", "NEG")
  )
  stripped <- removeTopFeats(fx, top)
  expect_false("feat_inf_1" %in% colnames(stripped))
  expect_false("feat_noise_1" %in% colnames(stripped))
  expect_true(all(
    c("genome_id", "genome_drug.resistant_phenotype") %in% colnames(stripped)
  ))
})

test_that("removeTopFeats errors when no features are provided", {
  fx <- make_pipeline_fixture(n_per_class = 5)
  empty <- tibble::tibble(
    Variable = character(),
    Importance = numeric(),
    Sign = character()
  )
  expect_error(removeTopFeats(fx, empty), "No more feats")
})
