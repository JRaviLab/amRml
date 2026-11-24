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

test_that("parallel completion is reported before the remaining worker finishes", {
  skip_on_os("windows")
  tasks <- fxFitTasks(labels = "real")
  manifest <- .manifestStart(tasks$built$out_dir,
                            tasks$built$eligibility$profile$orb, list())
  local_mocked_bindings(.runTaskGroup = function(job, ...) {
    Sys.sleep(if (grepl("000001.rds$", job$checkpoint)) 1.2 else 3.5)
    .saveFitCheckpoint(rep(list(list(error = "test failure")), nrow(job$tasks)),
                       job$checkpoint)
    TRUE
  })
  intermediate <- FALSE
  record <- .manifestEvent
  parent <- Sys.getpid()
  local_mocked_bindings(.manifestEvent = function(path, message, ..., details = list()) {
    expect_equal(Sys.getpid(), parent)
    if (identical(as.integer(details$n_completed), 1L)) {
      intermediate <<- TRUE
      expect_false(file.exists(file.path(.runPaths(tasks$built$out_dir)$checkpoints,
                                         "group-000002.rds")))
    }
    record(path, message, ..., details = details)
  })
  fit <- fitModels(tasks, reuse_fits = FALSE, BPPARAM = BiocParallel::MulticoreParam(2))
  expect_true(intermediate)
  events <- jsonlite::read_json(manifest)$runs[[1]]$events
  expect_equal(tail(events, 1)[[1]]$details$n_completed, 2)
  expect_true(all(fit$tasks$status == "failed"))
})

test_that("interrupted installation resumes without fitting completed groups again", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks()
  install <- .restoreFitResults
  local_mocked_bindings(.restoreFitResults = function(...) stop("injected installation failure"))
  expect_error(fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam()), "injected installation")
  expect_true(dir.exists(.runPaths(tasks$built$out_dir)$checkpoints))
  local_mocked_bindings(.restoreFitResults = install,
                       .runTaskGroup = function(...) stop("completed groups must not refit"))
  resumed <- fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam())
  expect_true(all(resumed$tasks$status == "fitted"))
  expect_false(dir.exists(.runPaths(tasks$built$out_dir)$checkpoints))
  reopened <- fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam())
  expect_equal(fxComparableFit(resumed), fxComparableFit(reopened))
})

test_that("resume rejects changed matrices and completed output corruption", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  fit <- fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam())
  model_path <- file.path(.runPaths(tasks$built$out_dir)$features, "0001.rds")
  original <- readBin(model_path, "raw", n = file.info(model_path)$size)
  store <- readRDS(model_path)
  store$values[1, 1] <- store$values[1, 1] + 1
  saveRDS(store, model_path)
  expect_error(fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam()),
               class = "amrml_invalid_argument")
  writeBin(original, model_path)
  writeLines("corrupted", .runPaths(tasks$built$out_dir)$performance)
  expect_error(fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam()),
               class = "amrml_invalid_argument")
})

test_that("failed tasks do not discard successful peers", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  path <- .runPaths(tasks$built$out_dir)$matrix_index
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
  expect_false(dir.exists(.runPaths(tasks$built$out_dir)$checkpoints))
})

test_that("partial workflow artifacts block an accidental fresh run", {
  orb <- withr::local_tempdir()
  for (artifact in c("metadata", "results", "state")) {
    out <- withr::local_tempdir()
    dir.create(file.path(out, artifact))
    expect_error(.checkArgOutDir(out, orb, overwrite = FALSE),
                 class = "amrml_invalid_argument")
  }
})

test_that("groups share a worker-local feature-store cache", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  index <- readRDS(.runPaths(tasks$built$out_dir)$matrix_index)
  paths <- file.path(tasks$built$out_dir, c("group1.rds", "group2.rds"))
  cache <- new.env(parent = emptyenv())
  extra <- list(species = "Tst", model_dir = NULL, save_fits = FALSE)
  run <- function(i) suppressWarnings(do.call(.runTaskGroup, c(list(
    job = .fitGroupJob(i, tasks$tasks, tasks$built$matrices, index,
                       tasks$built$out_dir, paths[[i]]), cache = cache
  ), extra)))
  expect_true(run(1L))
  expect_gt(length(cache), 0L)
  source <- file.path(.runPaths(tasks$built$out_dir)$features, "0001.rds")
  expect_true(file.rename(source, paste0(source, ".moved")))
  expect_true(run(2L))
  expect_null(readRDS(paths[[2]])[[1]]$error)
})


test_that("year and country reuse avoid repeated preparation without changing results", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  for (mode in c("cross_year", "cross_country")) for (sparse in c(FALSE, TRUE)) {
    tasks <- fxFitTasks(sparse = sparse)
    tasks$built$matrices$mode_id <- mode
    counts <- new.env(parent = emptyenv())
    counts$data <- counts$importance <- 0L
    model_data <- .modelData
    importance <- .taskImportance
    local_mocked_bindings(
      .modelData = function(...) {
        counts$data <- counts$data + 1L
        model_data(...)
      },
      .taskImportance = function(...) {
        counts$importance <- counts$importance + 1L
        importance(...)
      }
    )
    reused <- fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam())
    expect_true(all(reused$tasks$status == "fitted"))
    expect_equal(counts$data, 6L) # Two training preparations, four test sets.
    expect_equal(counts$importance, 2L)
    counts$data <- counts$importance <- 0L
    independent <- fxRunFit(tasks, overwrite = TRUE, reuse_fits = FALSE,
                            BPPARAM = BiocParallel::SerialParam())
    expect_equal(counts$data, 8L)
    expect_equal(counts$importance, 4L)
    expect_equal(fxComparableFit(reused), fxComparableFit(independent))
    # Restore these mocks before the next iteration captures the original functions.
    testthat::local_mocked_bindings(.modelData = model_data, .taskImportance = importance)
  }
})

test_that("cross-test grouping checks exact training inputs and fitting settings", {
  tasks <- fxFitTasks(labels = "real")
  matrices <- tasks$built$matrices
  index <- readRDS(.runPaths(tasks$built$out_dir)$matrix_index)
  for (mode in c("cross_drug", "cross_year", "cross_country")) {
    matrices$mode_id <- mode
    expect_length(.fitTaskGroups(tasks$tasks, matrices, index), 1L)
    for (field in c("seed", "labels", "n_fold", "holdout")) {
      listed <- tasks$tasks
      listed[[field]][2] <- switch(field, seed = 79L, labels = "shuffled",
                                   n_fold = 3L, holdout = 0.3)
      expect_length(.fitTaskGroups(listed, matrices, index), 2L)
    }
    changed <- matrices
    changed$encoding[2] <- "binary"
    expect_length(.fitTaskGroups(tasks$tasks, changed, index), 2L)
    changed <- index
    changed$columns[[2]] <- rev(changed$columns[[2]])
    expect_length(.fitTaskGroups(tasks$tasks, matrices, changed), 2L)
    changed <- index
    changed$rows[[2]]$phenotype[1] <- "Susceptible"
    expect_length(.fitTaskGroups(tasks$tasks, matrices, changed), 2L)
    changed <- index
    train <- which(changed$rows[[2]]$role == "train")
    changed$rows[[2]][train, ] <- changed$rows[[2]][rev(train), ]
    expect_length(.fitTaskGroups(tasks$tasks, matrices, changed), 2L)
  }
})

test_that("resuming older cross-test checkpoints retains their group numbering", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  tasks$built$matrices$mode_id <- "cross_year"
  grouping <- .fitTaskGroups
  install <- .restoreFitResults
  local_mocked_bindings(.fitTaskGroups = function(listed, ...) as.list(seq_len(nrow(listed))),
                       .restoreFitResults = function(...) stop("injected installation failure"))
  expect_error(fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam()), "injected installation")
  local_mocked_bindings(.fitTaskGroups = grouping, .restoreFitResults = install,
                       .fitModel = function(...) stop("completed tasks must not refit"))
  fit <- fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam())
  expect_true(all(fit$tasks$status == "fitted"))
  reopened <- fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam())
  expect_equal(fxComparableFit(fit), fxComparableFit(reopened))
})
