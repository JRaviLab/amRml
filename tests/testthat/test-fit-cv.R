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
    expect_true(.penaltyCovered(parsnip::extract_fit_engine(sockets$fit)$lambda,
                                sockets$penalty))
  }
})

test_that("unfitted penalties cannot qualify through endpoint substitution", {
  d <- cvTestData(FALSE, FALSE)
  w <- workflows::workflow() |>
    workflows::add_recipe(.buildRecipe(d)) |>
    workflows::add_model(.modelSpec())
  folds <- rsample::vfold_cv(d, v = 2, strata = phenotype)
  grid <- tibble::tibble(penalty = c(1e-12, .1), mixture = 1)
  tuned <- .tuneModelFolds(w, folds, grid)
  checked <- .checkTuningPaths(tuned)
  expect_true(any(grepl("outside the fitted", checked$notes$note)))
  for (metrics in checked$tuned$.metrics) {
    expect_true(all(is.na(metrics$.estimate[metrics$penalty == 1e-12])))
  }
  expect_equal(tune::select_best(checked$tuned, metric = "mcc")$penalty, .1)
  expect_false(.penaltyCovered(c(10, 1), .1))
  path <- .extendPenaltyPath(c(10, 1), .1)
  expect_true(.penaltyCovered(path, .1))
  expect_true(all(c(10, 1, .1) %in% path))
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

test_that("an entirely unsupported penalty grid fails instead of selecting a clamped fit", {
  local_mocked_bindings(.tuningGrid = function() {
    tibble::tibble(penalty = 1e-12, mixture = 1)
  })
  expect_error(.fitModel(cvTestData(FALSE, FALSE), 2),
               "No tuning candidate has a defined MCC at a fitted penalty")
})
