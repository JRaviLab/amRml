# Bakes `data` with one scaling step, from a fixed seed; returns the data and the RNG state after.
fxBakeScaled <- function(data, add_step) {
  set.seed(1)
  recipe <- recipes::prep(add_step(recipes::recipe(phenotype ~ ., data = data)), training = data)
  list(data = recipes::bake(recipe, new_data = data), seed = .Random.seed)
}

test_that(".stepScale() bakes exactly as step_scale() and step_normalize()", {
  set.seed(42)
  values <- Matrix::rsparsematrix(40, 12, density = 0.3)
  values[, 1] <- 0
  phenotype <- rep(.PHENOTYPES, 20)

  sparse <- .modelData(values, phenotype, sparse = TRUE)
  ours <- fxBakeScaled(sparse, function(r) .stepScale(r, recipes::all_predictors()))
  # recipes warns about the constant column; ours sets its sd to 1 the same way.
  expect_warning(
    theirs <- fxBakeScaled(sparse, function(r) recipes::step_scale(r, recipes::all_predictors())),
    "zero variance"
  )
  expect_identical(ours, theirs)

  dense <- .modelData(values, phenotype)
  ours <- fxBakeScaled(dense, function(r) .stepScale(r, recipes::all_predictors(), center = TRUE))
  expect_warning(
    theirs <- fxBakeScaled(dense, function(r) recipes::step_normalize(r, recipes::all_predictors())),
    "zero variance"
  )
  expect_identical(ours, theirs)
})
