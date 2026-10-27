test_that("dense and sparse fits agree across serial, sockets and reuse settings", {
  skip_on_cran()
  # Socket workers load amRml from the library, so with devtools::load_all() they can't
  # find internal functions, or run an older installed version. Compare sockets only
  # when amRml itself is installed, as in R CMD check.
  use_sockets <- !exists(".__DEVTOOLS__", envir = asNamespace("amRml"), inherits = FALSE)
  for (sparse in c(FALSE, TRUE)) {
    tasks <- fxFitTasks(sparse = sparse)
    serial <- fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam())
    independent <- fxRunFit(tasks, overwrite = TRUE, reuse_fits = FALSE,
                           BPPARAM = BiocParallel::SerialParam())
    expect_true(all(serial$tasks$status == "fitted"))
    expect_equal(fxComparableFit(serial), fxComparableFit(independent))
    if (use_sockets) {
      sockets <- fxRunFit(tasks, overwrite = TRUE,
                         BPPARAM = BiocParallel::SnowParam(workers = 2L, tasks = 2L))
      expect_equal(fxComparableFit(serial), fxComparableFit(sockets))
    }
  }
})

test_that("interrupted installation resumes without fitting completed groups again", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks()
  install <- .restoreFitResults
  local_mocked_bindings(.restoreFitResults = function(...) stop("injected installation failure"))
  expect_error(fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam()), "injected installation")
  expect_true(dir.exists(file.path(tasks$built$out_dir, "fit-checkpoints")))
  local_mocked_bindings(.restoreFitResults = install,
                       .runTaskGroup = function(...) stop("completed groups must not refit"))
  resumed <- fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam())
  expect_true(all(resumed$tasks$status == "fitted"))
  expect_false(dir.exists(file.path(tasks$built$out_dir, "fit-checkpoints")))
  reopened <- fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam())
  expect_equal(fxComparableFit(resumed), fxComparableFit(reopened))
})

test_that("resume rejects changed matrices and completed output corruption", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  fit <- fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam())
  model_path <- file.path(tasks$built$out_dir, "features", "0001.rds")
  original <- readBin(model_path, "raw", n = file.info(model_path)$size)
  store <- readRDS(model_path)
  store$values[1, 1] <- store$values[1, 1] + 1
  saveRDS(store, model_path)
  expect_error(fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam()),
               class = "amrml_invalid_argument")
  writeBin(original, model_path)
  writeLines("corrupted", file.path(tasks$built$out_dir, "results", "performance.parquet"))
  expect_error(fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam()),
               class = "amrml_invalid_argument")
})

test_that("failed tasks do not discard successful peers", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  path <- file.path(tasks$built$out_dir, "matrix_index.rds")
  index <- readRDS(path)
  index$rows[[2]]$genome_id[[36]] <- "missing-genome"
  saveRDS(index, path)
  fit <- fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam())
  expect_equal(fit$tasks$status, c("fitted", "failed"))
  expect_equal(fit$tasks$rule_id[[2]], "fit_failed")
  expect_match(fit$tasks$message[[2]], "missing from feature store")
  expect_equal(nrow(fit$performance), 1L)
})

test_that("all-failed and empty plans still write typed result tables", {
  local_mocked_bindings(.fitModel = function(...) stop("injected model failure"))
  tasks <- fxFitTasks(labels = "real")
  fit <- fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam())
  expect_true(all(fit$tasks$status == "failed"))
  expect_equal(nrow(fit$performance), 0L)
  expect_type(fit$importance$Sign, "character")
  expect_true(all(c("genome_id", "truth") %in% names(fit$predictions)))
  tasks$tasks <- tasks$tasks[0, ]
  empty <- fxRunFit(tasks, overwrite = TRUE, BPPARAM = BiocParallel::SerialParam())
  expect_equal(nrow(empty$tasks), 0L)
  expect_equal(nrow(empty$performance), 0L)
})

test_that("installation rollback restores every previous output", {
  out <- withr::local_tempdir()
  stage <- withr::local_tempdir()
  dir.create(file.path(out, "results"))
  writeLines("previous tasks", file.path(out, "tasks.parquet"))
  writeLines("previous results", file.path(out, "results", "performance.parquet"))
  writeLines("new tasks", file.path(stage, "tasks.parquet"))
  before <- tools::md5sum(file.path(out, c("tasks.parquet", "results/performance.parquet")))
  expect_warning(expect_error(.restoreFitResults(stage, out, c("tasks.parquet", "results")),
                               "Could not install result"))
  expect_equal(tools::md5sum(names(before)), before)
})

test_that("fit flags are validated before creating checkpoint files", {
  tasks <- fxFitTasks()
  for (flag in c("resume", "reuse_fits", "overwrite")) {
    args <- c(list(tasks = tasks, BPPARAM = BiocParallel::SerialParam()),
              stats::setNames(list(NA), flag))
    expect_error(do.call(fitModels, args), class = "amrml_invalid_argument")
  }
  expect_false(dir.exists(file.path(tasks$built$out_dir, "fit-checkpoints")))
})

test_that("partial workflow artifacts block an accidental fresh run", {
  orb <- withr::local_tempdir()
  for (artifact in c("features", "matrix_index.rds", "fit-checkpoints", "workflow.rds")) {
    out <- withr::local_tempdir()
    file.create(file.path(out, artifact))
    expect_error(.checkArgOutDir(out, orb, overwrite = FALSE),
                 class = "amrml_invalid_argument")
  }
})

test_that("groups share a worker-local feature-store cache", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  index <- readRDS(file.path(tasks$built$out_dir, "matrix_index.rds"))
  paths <- file.path(tasks$built$out_dir, c("group1.rds", "group2.rds"))
  cache <- new.env(parent = emptyenv())
  extra <- if ("save_fits" %in% names(formals(.runTaskGroup))) {
    list(species = "Tst", model_dir = NULL, save_fits = FALSE)
  } else list()
  run <- function(i) suppressWarnings(do.call(.runTaskGroup, c(list(
    job = .fitGroupJob(i, tasks$tasks, tasks$built$matrices, index,
                       tasks$built$out_dir, paths[[i]]), cache = cache
  ), extra)))
  expect_true(run(1L))
  expect_gt(length(cache), 0L)
  source <- file.path(tasks$built$out_dir, "features", "0001.rds")
  expect_true(file.rename(source, paste0(source, ".moved")))
  expect_true(run(2L))
  expect_null(readRDS(paths[[2]])[[1]]$error)
})
