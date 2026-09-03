# Unit tests for functions in core_ml.R. These cover
# splitting, recipe/model/workflow/grid construction, and predict/eval helpers
# without running the full grid-tuning pipeline.

skip_if_missing_deps <- function() {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("rsample")
  skip_if_not_installed("recipes")
  skip_if_not_installed("parsnip")
}

test_that("splitMLInputTibble returns correct rsplit class per mode", {
  skip_if_missing_deps()
  fx <- make_pipeline_fixture()

  # Classical train/val/test split
  sp <- splitMLInputTibble(fx, split = c(0.6, 0.2), seed = 42)
  expect_s3_class(sp, "initial_validation_split")

  # Pure CV -> initial_split holdout
  sp_cv <- splitMLInputTibble(fx, split = c(1, 0), seed = 42)
  expect_s3_class(sp_cv, "initial_split")
})

test_that("splitMLInputTibble validates arguments", {
  skip_if_missing_deps()
  fx <- make_pipeline_fixture()
  expect_error(splitMLInputTibble(data.frame(x = 1)), "tibble")
  expect_error(splitMLInputTibble(fx, split = c(0.8, 0.3)), "Unsupported|<= 1")
  expect_error(splitMLInputTibble(fx, seed = "x"), "numeric")
})

test_that("buildRecipe returns a recipe with zv + normalize steps", {
  skip_if_missing_deps()
  fx <- make_pipeline_fixture()
  rec <- buildRecipe(fx)
  expect_s3_class(rec, "recipe")
  step_types <- vapply(rec$steps, function(s) class(s)[1], character(1))
  expect_true("step_zv" %in% step_types)
  expect_true("step_normalize" %in% step_types)
  expect_false("step_pca" %in% step_types)

  rec_pca <- buildRecipe(fx, use_pca = TRUE, pca_threshold = 0.9)
  step_types_pca <- vapply(rec_pca$steps, function(s) class(s)[1], character(1))
  expect_true("step_pca" %in% step_types_pca)
})

test_that("buildLRModel switches between binary and multi-class", {
  skip_if_missing_deps()
  lr <- buildLRModel(multi_class = FALSE)
  expect_s3_class(lr, "logistic_reg")

  mn <- buildLRModel(multi_class = TRUE)
  expect_s3_class(mn, "multinom_reg")

  expect_error(buildLRModel(multi_class = "yes"), "logical")
})

test_that("buildWflow assembles a workflow with model + recipe", {
  skip_if_missing_deps()
  fx <- make_pipeline_fixture()
  wflow <- buildWflow(buildLRModel(), buildRecipe(fx))
  expect_s3_class(wflow, "workflow")

  # Non-parsnip inputs trigger an error in .checkArgParsnipMod (though the
  # internal check doesn't guard against single-class inputs cleanly).
  expect_error(buildWflow("not a model", buildRecipe(fx)))
  expect_error(buildWflow(buildLRModel(), "not a recipe"), "recipe")
})

test_that("buildTuningGrid produces penalty x mixture cross-product", {
  grid <- buildTuningGrid(
    model = "LR",
    penalty_vec = c(0.01, 0.1, 1),
    mix_vec = c(0, 0.5, 1)
  )
  expect_s3_class(grid, "tbl_df")
  expect_identical(nrow(grid), 9L)
  expect_setequal(colnames(grid), c("penalty", "mixture"))
  expect_setequal(unique(grid$penalty), c(0.01, 0.1, 1))
  expect_setequal(unique(grid$mixture), c(0, 0.5, 1))

  expect_error(buildTuningGrid(model = "SVM"), "one of")
})

test_that("predictML augments test data with .pred_class", {
  skip_if_missing_deps()
  fx <- make_pipeline_fixture()
  set.seed(7)

  # Fit a LR so we can exercise predict/confmat helpers.
  mod <- parsnip::logistic_reg(penalty = 0.01, mixture = 0) |>
    parsnip::set_engine("glmnet")
  rec <- buildRecipe(fx)
  wflow <- buildWflow(mod, rec)
  fit <- wflow |> parsnip::fit(data = fx)

  preds <- predictML(fit, fx)
  expect_true(".pred_class" %in% colnames(preds))
  expect_s3_class(preds$.pred_class, "factor")

  cm <- getConfusionMatrix(preds)
  expect_s3_class(cm, "conf_mat")

  # log2(AUPRC/prior) warns on balanced classes; unrelated to what we test.
  mets <- suppressWarnings(calculateEvalMets(preds))
  # Returns a length-6 vector: f1, auprc, bal_acc, mcc, nmcc, log2_apop.
  expect_length(mets, 6)
  expect_true(all(!is.na(mets[1:4])))
})

test_that("extractTopFeats rejects conflicting / missing selectors", {
  skip_if_missing_deps()
  fx <- make_pipeline_fixture()
  mod <- parsnip::logistic_reg(penalty = 0.01, mixture = 0) |>
    parsnip::set_engine("glmnet")
  fit <- buildWflow(mod, buildRecipe(fx)) |> parsnip::fit(data = fx)

  expect_error(
    extractTopFeats(fit, prop_vi_top_feats = NA, n_top_feats = NA),
    "specify either"
  )
  # With both numeric, n_top_feats takes precedence (prop gets NA'd internally).
  top <- extractTopFeats(fit, prop_vi_top_feats = c(0.1, 0.2), n_top_feats = 2)
  expect_s3_class(top, "tbl_df")
  expect_setequal(colnames(top), c("Variable", "Importance", "Sign"))
  expect_lte(nrow(top), 2)
})

# A multi-class (MDR-style) fit, where glmnet returns one coefficient matrix per
# class. This is the case that must not reach `.viGlmnet()`.
make_multiclass_fit <- function() {
  set.seed(42)
  n <- 30
  mc <- tibble::tibble(
    genome_id = paste0("g", seq_len(n)),
    resistant_classes = factor(rep(c("A", "B", "C"), each = n / 3))
  )
  for (j in 1:8) mc[[paste0("feat_", j)]] <- as.integer(rbinom(n, 1, 0.5))

  mod <- parsnip::multinom_reg(penalty = 0.01, mixture = 0.5) |>
    parsnip::set_engine("glmnet")
  buildWflow(mod, buildRecipe(mc)) |> parsnip::fit(data = mc)
}

test_that("extractTopFeats handles multi-class fits with per-class columns", {
  skip_if_missing_deps()
  fit <- make_multiclass_fit()

  expect_warning(
    top <- extractTopFeats(fit, n_top_feats = 5),
    "multi-class"
  )

  # One column per class, not the binary Variable/Importance/Sign triple.
  expect_setequal(colnames(top), c("Variable", "A", "B", "C"))
  expect_gt(nrow(top), 0)
})

test_that(".viGlmnet reports coefficients from the tuned model", {
  skip_if_missing_deps()
  fx <- make_pipeline_fixture()
  # LASSO, so some coefficients are driven to exactly zero.
  mod <- parsnip::logistic_reg(penalty = 0.01, mixture = 1) |>
    parsnip::set_engine("glmnet")
  fit <- buildWflow(mod, buildRecipe(fx)) |> parsnip::fit(data = fx)

  vi <- .viGlmnet(fit)
  expect_setequal(colnames(vi), c("Variable", "Importance", "Sign"))
  expect_true(all(vi$Sign %in% c("POS", "NEG")))
  # Downstream slicing relies on the decreasing sort.
  expect_equal(vi$Importance, sort(vi$Importance, decreasing = TRUE))

  # Importance is |coefficient| at the tuned penalty, with zeros dropped.
  gf <- parsnip::extract_fit_engine(fit)
  tuned <- as.numeric(.getFitHps(fit)["penalty"])
  expected <- stats::coef(gf, s = tuned)[, 1, drop = TRUE]
  expected <- expected[names(expected) != "(Intercept)"]
  expected <- expected[expected != 0]
  expect_equal(vi$Importance, sort(abs(unname(expected)), decreasing = TRUE))
  expect_lt(nrow(vi), length(coef(gf, s = tuned)) - 1)
})


test_that(".viGlmnet drops zeroed-out features", {
  # Mock the three external calls so a zero coefficient is guaranteed, rather
  # than depending on a real fit happening to produce one.
  local_mocked_bindings(
    extract_fit_engine = function(fit) list(lambda = c(0.5, 0.1, 0.01)),
    .package = "parsnip"
  )
  local_mocked_bindings(
    .getFitHps = function(fit) c(penalty = 0.01, mixture = 1)
  )
  local_mocked_bindings(
    coef = function(object, s, ...) {
      matrix(
        c(1, 2, -3, 0),
        dimnames = list(
          c("(Intercept)", "pos_feat", "neg_feat", "zero_feat"), NULL
        )
      )
    },
    .package = "stats"
  )

  vi <- .viGlmnet(fit = "placeholder")

  expect_setequal(vi$Variable, c("pos_feat", "neg_feat"))
  expect_equal(unname(vi$Sign[vi$Variable == "pos_feat"]), "POS")
  expect_equal(unname(vi$Sign[vi$Variable == "neg_feat"]), "NEG")
  expect_false("zero_feat" %in% vi$Variable)
})

test_that("extractTopFeats includes the last feature at prop_vi_top_feats = c(0, 1)", {
  skip_if_missing_deps()
  fx <- make_pipeline_fixture()
  mod <- parsnip::logistic_reg(penalty = 0.01, mixture = 0) |>
    parsnip::set_engine("glmnet")
  fit <- buildWflow(mod, buildRecipe(fx)) |> parsnip::fit(data = fx)

  # The default prop_vi_top_feats = c(0, 1) is documented to return all
  # features; the strict "<" upper-bound comparison used to drop the least
  # important one, since its cumulative importance exactly equals the total.
  all_feats <- .viGlmnet(fit)
  top <- extractTopFeats(fit, prop_vi_top_feats = c(0, 1), n_top_feats = NA)

  expect_equal(nrow(top), nrow(all_feats))
  expect_setequal(top$Variable, all_feats$Variable)
})
make_log2apop_fixture <- function(n_resistant, n_susceptible, seed = 1) {
  set.seed(seed)
  phenotype <- factor(
    c(rep("Resistant", n_resistant), rep("Susceptible", n_susceptible)),
    levels = c("Resistant", "Susceptible")
  )
  pred_resistant <- ifelse(
    phenotype == "Resistant",
    stats::runif(length(phenotype), 0.6, 0.95),
    stats::runif(length(phenotype), 0.05, 0.4)
  )
  tibble::tibble(
    genome_drug.resistant_phenotype = phenotype,
    .pred_class = factor(
      ifelse(pred_resistant > 0.5, "Resistant", "Susceptible"),
      levels = c("Resistant", "Susceptible")
    ),
    .pred_Resistant = pred_resistant
  )
}

test_that(".calculateLog2APOP warns (not just messages) on roughly balanced classes", {
  fx <- make_log2apop_fixture(n_resistant = 5, n_susceptible = 5)
  expect_warning(amRml:::.calculateLog2APOP(fx), "roughly balanced")
})

test_that(".calculateLog2APOP warns (not just messages) on imbalanced classes", {
  fx <- make_log2apop_fixture(n_resistant = 9, n_susceptible = 1)
  expect_warning(amRml:::.calculateLog2APOP(fx), "imbalanced")

})
