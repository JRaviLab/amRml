cvTestData <- function(multi, sparse) {
  withr::with_seed(71, {
    labels <- if (multi) c("a", "b", "c") else .PHENOTYPES
    y <- factor(rep(labels, each = 30), levels = labels)
    p <- if (sparse) 30L else 8L
    x <- matrix(if (sparse) rbinom(length(y) * p, 1, .06) else
                  rnorm(length(y) * p), nrow = length(y))
    x[, 1] <- 3 * as.integer(y)
    x[, 2] <- 2 * (y == labels[1])
    .modelData(Matrix::Matrix(x, sparse = TRUE), y, sparse, labels)
  })
}

test_that("serial and socket CV conserve binary and multinomial fits", {
  for (multi in c(FALSE, TRUE)) for (sparse in c(FALSE, TRUE)) {
    d <- cvTestData(multi, sparse)
    serial <- withr::with_seed(123, .fitModel(d, 3, sparse, multi, cv_workers = 1))
    sockets <- withr::with_seed(123, .fitModel(d, 3, sparse, multi, cv_workers = 2))
    expect_equal(serial$penalty, sockets$penalty)
    expect_equal(serial$mixture, sockets$mixture)
    expect_equal(predict(serial$fit, d), predict(sockets$fit, d))
    expect_equal(stats::coef(parsnip::extract_fit_engine(serial$fit), s = serial$penalty),
                 stats::coef(parsnip::extract_fit_engine(sockets$fit), s = sockets$penalty))
  }
})

test_that("penalty diagnostics preserve the full legacy tuning results and Ridge", {
  d <- cvTestData(FALSE, FALSE)
  w <- workflows::workflow() |>
    workflows::add_recipe(.buildRecipe(d)) |>
    workflows::add_model(.modelSpec())
  folds <- rsample::vfold_cv(d, v = 2, strata = phenotype)
  tuned <- .tuneModelFolds(w, folds, .tuningGrid())
  before <- tune::collect_metrics(tuned)
  notes <- .tuningPathNotes(tuned)
  expect_true(any(grepl("outside the fitted", notes$note)))
  expect_identical(tune::collect_metrics(tuned), before)
  expect_equal(nrow(before), 60L)
  expect_equal(sum(before$mixture == 0 & is.finite(before$mean)), 10L)
  expect_true(all(is.finite(before$mean)))
})

test_that("CV parallelism is off by default and rejects non-boolean flags", {
  expect_identical(formals(fitModels)$parallel_cv, FALSE)
  expect_identical(formals(runModelingWorkflow)$parallel_cv, FALSE)
  tasks <- fxFitTasks()
  expect_error(fitModels(tasks, parallel_cv = 2), class = "amrml_invalid_argument")
})

test_that("parallel CV uses one outer task and a shared CPU budget for every family", {
  tasks <- fxFitTasks(labels = "real")
  tasks$built$matrices$mode_id <- "vanilla"
  tasks$built$matrices$unit <- c("drug", "mdr")
  tasks$built$matrices$target <- c("DRA", "resistant_classes")
  parent <- Sys.getpid()
  local_mocked_bindings(.runTaskGroup = function(job, ...) {
    expect_identical(Sys.getpid(), parent)
    expect_equal(job$cv_workers, 2)
    .saveFitCheckpoint(rep(list(list(error = "fixture")), nrow(job$tasks)), job$checkpoint)
    TRUE
  })
  fitModels(tasks, parallel_cv = TRUE, BPPARAM = BiocParallel::SnowParam(2))
})

test_that("prepared arrays retain the CV toggle and collect equivalent results", {
  tasks <- fxFitTasks(labels = "real")
  serial <- fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam())
  plan <- fxRunFit(tasks, overwrite = TRUE, execution = "prepare", n_chunks = 1,
                   parallel_cv = TRUE, BPPARAM = BiocParallel::SnowParam(2))
  expect_true(plan$parallel_cv)
  withr::local_envvar(SLURM_CPUS_PER_TASK = "2")
  suppressWarnings(runTaskChunk(tasks$built$out_dir, 1, threads = 2))
  collected <- collectTaskArray(tasks$built$out_dir)
  expect_equal(fxComparableFit(serial), fxComparableFit(collected))
})

test_that("below-path penalties retain legacy endpoint coefficients in the final fit", {
  local_mocked_bindings(.tuningGrid = function() {
    tibble::tibble(penalty = 1e-12, mixture = 0)
  })
  d <- cvTestData(FALSE, FALSE)
  result <- withr::with_seed(123, .fitModel(d, 2))
  legacy <- workflows::workflow() |>
    workflows::add_recipe(.buildRecipe(d)) |>
    workflows::add_model(tune::finalize_model(.modelSpec(),
      tibble::tibble(penalty = 1e-12, mixture = 0))) |>
    parsnip::fit(data = d)
  expect_equal(result$penalty, 1e-12)
  expect_equal(result$mixture, 0)
  expect_false(.penaltyCovered(parsnip::extract_fit_engine(result$fit)$lambda,
                              result$penalty))
  expect_equal(stats::coef(parsnip::extract_fit_engine(result$fit), s = result$penalty),
               stats::coef(parsnip::extract_fit_engine(legacy), s = result$penalty))
  expect_equal(predict(result$fit, d), predict(legacy, d))
  expect_true(any(result$notes$stage == "tune"))
  expect_true(any(result$notes$stage == "final_fit"))
  expect_true(all(grepl("compatibility", result$notes$note)))
})
