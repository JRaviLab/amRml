test_that("saved models survive resume and reproduce recorded predictions", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  install <- .restoreFitResults
  local_mocked_bindings(.restoreFitResults = function(...) stop("injected installation failure"))
  expect_error(fitModels(tasks, BPPARAM = BiocParallel::SerialParam()), "injected installation")
  expect_error(fitModels(tasks, resume = TRUE, save_fits = FALSE,
                         BPPARAM = BiocParallel::SerialParam()),
               class = "amrml_invalid_argument")
  local_mocked_bindings(.restoreFitResults = install,
                       .runTaskGroup = function(...) stop("must not refit"))
  fit <- suppressWarnings(fitModels(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam()))
  expect_true(file.exists(file.path(tasks$built$out_dir, "run.json")))
  expect_true(all(fit$performance$species == "Tst"))
  expect_true(all(fit$importance$species == "Tst"))
  index <- readRDS(file.path(tasks$built$out_dir, "matrix_index.rds"))
  cache <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(fit$tasks))) {
    id <- fit$tasks$task_id[[i]]
    saved <- readRDS(file.path(tasks$built$out_dir, "results", "models", paste0(id, ".rds")))
    matrix <- .readMatrixWide(file.path(tasks$built$out_dir, "features", "0001.rds"),
                              index$rows[[i]], index$columns[[i]], "counts", cache)
    test <- matrix$rows$role == "test"
    data <- .modelData(matrix$values[test, , drop = FALSE], matrix$rows$phenotype[test])
    scored <- parsnip::augment(saved$fit, data)
    expect_equal(saved$features, matrix$features)
    expect_equal(scored$.pred_Resistant,
                 fit$predictions$prob_resistant[fit$predictions$task_id == id])
  }
})

test_that("save_fits FALSE leaves no saved models and can reopen results", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  fit <- suppressWarnings(fitModels(tasks, save_fits = FALSE,
                                    BPPARAM = BiocParallel::SerialParam()))
  expect_false(dir.exists(file.path(tasks$built$out_dir, "results", "models")))
  local_mocked_bindings(.runTaskGroup = function(...) stop("must not refit"))
  expect_equal(fxComparableFit(fit),
               fxComparableFit(fitModels(tasks, save_fits = FALSE, resume = TRUE,
                                         BPPARAM = BiocParallel::SerialParam())))
})

test_that("zero returned coefficients retain the export schema", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid,
                       .taskImportance = function(...) tibble::tibble(
                         Variable = character(), Importance = numeric(), Sign = character()))
  tasks <- fxFitTasks(labels = "real")
  tasks$built$matrices$mode_id <- "vanilla"
  fit <- suppressWarnings(fitModels(tasks, save_fits = FALSE,
                                    BPPARAM = BiocParallel::SerialParam()))
  expect_true(all(fit$tasks$status == "fitted"))
  expect_true(all(fit$performance$n_feats_returned == 0L))
  expect_equal(nrow(fit$importance), 0L)
  expect_type(fit$importance$species, "character")
  expect_no_error(exportMergedResults(fit))
  top <- arrow::read_parquet(file.path(tasks$built$out_dir, "merged", "all_top_features.parquet"))
  expect_equal(nrow(top), 0L)
  expect_true(all(c("species", "Variable", "Sign") %in% names(top)))
})

test_that("workflow resume finishes a failed export without refitting", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  metadata <- fxRows(10, 10)
  orb <- fxOrb(metadata = metadata, parquets = list(gene_count = fxSignalGenes(10, 10)))
  root <- withr::local_tempdir()
  withr::local_dir(root)
  args <- list(orb_dir = orb, out_dir = "run", modes = "vanilla", seeds = 31L,
               min_genomes = 4, n_fold = 2, save_fits = FALSE,
               BPPARAM = BiocParallel::SerialParam())
  exporter <- exportMergedResults
  local_mocked_bindings(exportMergedResults = function(...) stop("injected export failure"))
  expect_error(suppressWarnings(do.call(runModelingWorkflow, args)), "injected export")
  local_mocked_bindings(exportMergedResults = exporter,
                       .fitModel = function(...) stop("completed models must not refit"))
  fit <- do.call(runModelingWorkflow, c(args, list(resume = TRUE)))
  expect_true(all(fit$tasks$status == "fitted"))
  expect_equal(fit$built$out_dir, normalizePath("run"))
  expect_true(all(file.exists(file.path("run", "merged",
                     c("all_perf.parquet", "all_top_features.parquet")))))
})

test_that("saved sparse workflows reproduce held-out predictions", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(sparse = TRUE, labels = "real")
  fit <- suppressWarnings(fitModels(tasks, BPPARAM = BiocParallel::SerialParam()))
  index <- readRDS(file.path(tasks$built$out_dir, "matrix_index.rds"))
  for (i in seq_len(nrow(fit$tasks))) {
    id <- fit$tasks$task_id[[i]]
    saved <- readRDS(file.path(tasks$built$out_dir, "results", "models", paste0(id, ".rds")))
    matrix <- .readMatrixWide(file.path(tasks$built$out_dir, "features", "0001.rds"),
                              index$rows[[i]], index$columns[[i]], "counts",
                              new.env(parent = emptyenv()))
    test <- matrix$rows$role == "test"
    data <- .modelData(matrix$values[test, , drop = FALSE], matrix$rows$phenotype[test],
                       sparse = TRUE)
    scored <- parsnip::augment(saved$fit, data)
    expect_equal(scored$.pred_Resistant,
                 fit$predictions$prob_resistant[fit$predictions$task_id == id])
  }
})
