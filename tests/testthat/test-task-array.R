test_that("chunk maps cover complete groups and validate array indices", {
  tasks <- fxFitTasks()
  groups <- splitTaskGroups(tasks, 100)
  expect_length(groups, 2L)
  expect_identical(sort(unlist(groups)), 1:2)
  expect_identical(taskChunk(groups, "2"), groups[[2]])
  expect_length(splitTaskGroups(tasks, 1), 1)
  expect_length(splitTaskGroups(tasks, reuse_fits = FALSE), 4)
  for (bad in list(NA, Inf, 1.5, 0, numeric(), c(1, 2), "1.5", "abc")) {
    expect_error(taskChunk(groups, bad), class = "amrml_invalid_argument")
  }
  for (bad in list(NA, Inf, 1.5, 0)) {
    expect_error(splitTaskGroups(tasks, bad), class = "amrml_invalid_argument")
  }
  tasks$tasks <- tasks$tasks[0, ]
  expect_identical(splitTaskGroups(tasks), list())
  tasks <- fxFitTasks(labels = "real")
  expect_identical(splitTaskGroups(tasks), list(1L))
})

test_that("batched workers preserve shared state, retry and collect equivalent fits", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks()
  serial <- fxRunFit(tasks, BPPARAM = BiocParallel::SerialParam())
  paths <- .runPaths(tasks$built$out_dir)
  manifest <- .manifestStart(tasks$built$out_dir, tasks$built$eligibility$profile$orb, list())
  before <- tools::md5sum(c(paths$tasks, paths$performance, manifest))
  plan <- fxRunFit(tasks, execution = "prepare", n_chunks = 1, overwrite = TRUE,
                   BPPARAM = BiocParallel::SerialParam())
  expect_length(plan$chunks, 1)
  expect_length(plan$chunks[[1]], 2)
  expect_error(fxRunFit(tasks, resume = TRUE, BPPARAM = BiocParallel::SerialParam()),
               "array run")
  expect_error(collectTaskArray(tasks$built$out_dir), "missing task groups")
  suppressWarnings(runTaskChunk(tasks$built$out_dir, 1))
  expect_true(dir.exists(paths$features))
  expect_true(dir.exists(paths$checkpoints))
  expect_identical(tools::md5sum(names(before)), before)
  checkpoints <- file.path(paths$checkpoints, sprintf("group-%06d.rds", 1:2))
  expect_true(all(file.exists(checkpoints)))
  checkpoint_before <- file.info(checkpoints)$mtime
  suppressWarnings(runTaskChunk(tasks$built$out_dir, 1))
  expect_identical(file.info(checkpoints)$mtime, checkpoint_before)
  local_mocked_bindings(.fitModel = function(...) stop("collector must not fit"))
  collected <- withVisible(collectTaskArray(tasks$built$out_dir))
  expect_false(collected$visible)
  fit <- collected$value
  expect_equal(fxComparableFit(fit), fxComparableFit(serial))
  expect_false(dir.exists(paths$checkpoints))
  expect_length(list.files(paths$models), 4)
})

test_that("collection refuses active, missing and corrupt groups", {
  tasks <- fxFitTasks()
  local_mocked_bindings(.runTaskGroup = function(job, ...) {
    .saveFitCheckpoint(rep(list(list(error = "test failure")), nrow(job$tasks)), job$checkpoint)
    TRUE
  })
  fitModels(tasks, execution = "prepare", n_chunks = 2, save_state = TRUE)
  out <- tasks$built$out_dir
  lock <- file.path(.runPaths(out)$checkpoints, "chunk-000001.lock")
  dir.create(lock)
  expect_error(runTaskChunk(out, 1), "already running")
  expect_error(collectTaskArray(out), "still active")
  unlink(lock, recursive = TRUE)
  runTaskChunk(out, 1)
  expect_error(collectTaskArray(out), "missing task groups")
  runTaskChunk(out, 2)
  checkpoint <- file.path(.runPaths(out)$checkpoints, "group-000002.rds")
  saveRDS("corrupt", checkpoint)
  expect_error(collectTaskArray(out), "Invalid fitted task checkpoint")
  expect_true(dir.exists(.runPaths(out)$features))
  expect_false(dir.exists(.runPaths(out)$results))
})

test_that("array settings and CPU budgets are checked before fitting", {
  tasks <- fxFitTasks()
  for (bad in list(numeric(), Inf, NA, 0, 1.2)) {
    expect_error(fitModels(tasks, group_subset = bad, execution = "worker", resume = TRUE),
                 class = "amrml_invalid_argument")
  }
  expect_error(fitModels(tasks, group_subset = 1), "required only")
  expect_error(fitModels(tasks, execution = "worker", group_subset = 1), "resume")
  expect_error(fitModels(tasks, execution = "collect"), "resume")
  withr::local_envvar(SLURM_CPUS_PER_TASK = "2")
  expect_error(runTaskChunk(tasks$built$out_dir, 1, threads = 4), "exceeds")
})

test_that("prepared workflow arrays finish merged exports and cleanup", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid, .registerManifest = function(path) path)
  orb <- fxOrb(metadata = fxRows(20, 20), parquets = list(gene_count = fxSignalGenes(20, 20)))
  out <- file.path(withr::local_tempdir(), "run")
  plan <- suppressWarnings(runModelingWorkflow(
    orb_dir = orb, out_dir = out, model_types = "vanilla", seeds = 31L,
    min_genomes = 4, n_fold = 2, feature_types = "gene", execution = "prepare", n_chunks = 2
  ))
  paths <- .runPaths(out)
  expect_false(file.exists(paths$tasks))
  before <- tools::md5sum(paths$manifest)
  for (id in rev(seq_along(plan$chunks))) suppressWarnings(runTaskChunk(out, id))
  expect_identical(tools::md5sum(paths$manifest), before)
  exporter <- exportMergedResults
  local_mocked_bindings(exportMergedResults = function(...) stop("injected export failure"))
  expect_error(collectTaskArray(out), "injected export failure")
  local_mocked_bindings(exportMergedResults = exporter,
                       .fitModel = function(...) stop("completed groups must not refit"))
  fit <- collectTaskArray(out)
  expect_true(all(fit$tasks$status == "fitted"))
  expect_true(file.exists(file.path(paths$merged, "all_perf.parquet")))
  expect_false(dir.exists(paths$state))
  manifest <- jsonlite::read_json(paths$manifest)
  expect_identical(vapply(manifest$runs, `[[`, character(1), "status"), c("prepared", "failed", "success"))
  expect_identical(manifest$artifacts$amRviz_input$status, "ready")
})


test_that("separate concurrent workers retain isolated checkpoints and failed task identity", {
  skip_on_os("windows")
  tasks <- fxFitTasks()
  listed <- tasks$tasks
  local_mocked_bindings(.runTaskGroup = function(job, ...) {
    .saveFitCheckpoint(lapply(job$tasks$task_id, function(id) list(error = paste("failure", id))),
                       job$checkpoint)
    TRUE
  })
  plan <- fitModels(tasks, execution = "prepare", n_chunks = 2, save_state = TRUE)
  out <- tasks$built$out_dir
  workers <- lapply(seq_along(plan$chunks), function(id) {
    parallel::mcparallel(runTaskChunk(out, id), mc.set.seed = FALSE)
  })
  results <- parallel::mccollect(workers)
  expect_false(any(vapply(results, inherits, logical(1), "try-error")))
  fit <- collectTaskArray(out)
  expect_identical(fit$tasks$task_id, listed$task_id)
  expect_identical(fit$tasks$message, paste("failure", listed$task_id))
  expect_true(all(fit$tasks$status == "failed"))
})

test_that("empty array plans collect without scheduling a worker", {
  tasks <- fxFitTasks()
  tasks$tasks <- tasks$tasks[0, ]
  plan <- fitModels(tasks, execution = "prepare")
  expect_identical(plan$chunks, list())
  fit <- collectTaskArray(tasks$built$out_dir)
  expect_equal(nrow(fit$tasks), 0L)
  expect_false(dir.exists(.runPaths(tasks$built$out_dir)$state))
})


test_that("the default workflow stays local when SLURM variables are present", {
  withr::local_envvar(SLURM_ARRAY_TASK_ID = "not-an-array-id",
                     SLURM_CPUS_PER_TASK = "not-a-cpu-count")
  local_mocked_bindings(
    .tuningGrid = fxSmallGrid, .registerManifest = function(path) path,
    runTaskChunk = function(...) stop("default workflow must stay local"),
    collectTaskArray = function(...) stop("default workflow must stay local")
  )
  orb <- fxOrb(metadata = fxRows(20, 20), parquets = list(gene_count = fxSignalGenes(20, 20)))
  out <- file.path(withr::local_tempdir(), "local-run")
  fit <- suppressWarnings(runModelingWorkflow(
    orb_dir = orb, out_dir = out, model_types = "vanilla", seeds = 31L,
    min_genomes = 4, n_fold = 2, feature_types = "gene", save_state = TRUE
  ))
  expect_s3_class(fit, "amr_fit")
  expect_true(all(fit$tasks$status == "fitted"))
  manifest <- jsonlite::read_json(.runPaths(out)$manifest)
  stages <- manifest$runs[[1]]$stages
  fitting <- Filter(function(x) identical(x$name, "fit_models"), stages)[[1]]
  expect_identical(fitting$parameters$backend, "SerialParam")
  expect_true(file.exists(file.path(.runPaths(out)$merged, "all_perf.parquet")))
  expect_false(file.exists(file.path(.runPaths(out)$checkpoints, "array.rds")))
})


test_that("lock metadata identifies owners without removing existing locks", {
  lock <- file.path(withr::local_tempdir(), "chunk.lock")
  withr::local_envvar(SLURM_ARRAY_JOB_ID = "321", SLURM_ARRAY_TASK_ID = "7")
  .acquireFitLock(lock, "Chunk")
  owner <- jsonlite::read_json(file.path(lock, "owner.json"))
  expect_identical(owner$array_job_id, "321")
  expect_identical(owner$array_task_id, "7")
  expect_equal(owner$pid, Sys.getpid())
  expect_error(.acquireFitLock(lock, "Chunk"), "321")
  expect_true(dir.exists(lock))
})

test_that("allocation failures retain task identity and durable diagnostics", {
  tasks <- fxFitTasks()
  local_mocked_bindings(.fitTask = function(...) stop("cannot allocate vector of size 8 Gb"))
  fit <- suppressWarnings(fitModels(tasks, save_state = FALSE))
  expect_true(all(fit$tasks$status == "failed"))
  expect_true(all(fit$tasks$rule_id == "fit_oom"))
  expect_identical(fit$tasks$task_id, tasks$tasks$task_id)
  directory <- file.path(tasks$built$out_dir, "diagnostics")
  lines <- unlist(lapply(list.files(directory, full.names = TRUE), readLines))
  records <- lapply(lines, jsonlite::fromJSON)
  events <- vapply(records, `[[`, character(1), "event")
  expect_equal(sum(events == "oom_error"), nrow(tasks$tasks))
  starts <- records[events == "group_start"]
  expect_true(all(vapply(starts, function(x) length(x$matrix_ids) > 0L, logical(1))))
  expect_true(all(vapply(starts, function(x) is.numeric(x$r_heap_used_mb), logical(1))))
  expect_false(dir.exists(.runPaths(tasks$built$out_dir)$state))
})


test_that("memory-limit errors are recognized and diagnostic failures do not replace fit errors", {
  expect_true(.isMemoryError("vector memory limit of 16.0 Gb reached, see mem.maxVSize()"))
  expect_true(.isMemoryError("std::bad_alloc"))
  expect_false(.isMemoryError("Model fitting failed: no usable predictors"))
  blocker <- tempfile(tmpdir = withr::local_tempdir())
  writeLines("file, not a directory", blocker)
  withr::local_options(warn = 2)
  expect_message(.fitDiagnostic(file.path(blocker, "logs"), "task_error"),
                 "Could not write fit diagnostics")
})
