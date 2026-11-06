test_that("saved models survive resume and reproduce recorded predictions", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  install <- .restoreFitResults
  local_mocked_bindings(.restoreFitResults = function(...) stop("injected installation failure"))
  expect_error(fitModels(tasks, save_fits = TRUE, save_state = TRUE,
                         BPPARAM = BiocParallel::SerialParam()), "injected installation")
  expect_error(fitModels(tasks, resume = TRUE, save_fits = FALSE,
                         BPPARAM = BiocParallel::SerialParam()),
               class = "amrml_invalid_argument")
  local_mocked_bindings(.restoreFitResults = install,
                       .runTaskGroup = function(...) stop("must not refit"))
  fit <- suppressWarnings(fitModels(
    tasks, resume = TRUE, save_fits = TRUE, save_state = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  ))
  expect_true(file.exists(.runPaths(tasks$built$out_dir)$completion))
  expect_true(all(fit$performance$species == "Tst"))
  expect_true(all(fit$importance$species == "Tst"))
  index <- readRDS(.runPaths(tasks$built$out_dir)$matrix_index)
  cache <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(fit$tasks))) {
    id <- fit$tasks$task_id[[i]]
    saved <- readRDS(file.path(tasks$built$out_dir, "results", "models", paste0(id, ".rds")))
    matrix <- .readMatrixWide(file.path(.runPaths(tasks$built$out_dir)$features, "0001.rds"),
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
  fit <- suppressWarnings(fitModels(tasks, save_fits = FALSE, save_state = TRUE,
                                    BPPARAM = BiocParallel::SerialParam()))
  expect_false(dir.exists(file.path(tasks$built$out_dir, "results", "models")))
  local_mocked_bindings(.runTaskGroup = function(...) stop("must not refit"))
  expect_equal(fxComparableFit(fit),
               fxComparableFit(fitModels(tasks, save_fits = FALSE, save_state = TRUE,
                                         resume = TRUE,
                                         BPPARAM = BiocParallel::SerialParam())))
})

test_that("successful fits remove models and RDS state by default", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  fit <- suppressWarnings(fitModels(tasks, BPPARAM = BiocParallel::SerialParam()))
  paths <- .runPaths(tasks$built$out_dir)
  expect_true(all(fit$tasks$status == "fitted"))
  expect_false(dir.exists(paths$models))
  expect_false(dir.exists(paths$state))
  expect_true(all(file.exists(c(paths$tasks, paths$performance,
                                paths$importance, paths$predictions))))
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
  top <- arrow::read_parquet(file.path(
    .runPaths(tasks$built$out_dir)$merged, "all_top_features.parquet"
  ))
  expect_equal(nrow(top), 0L)
  expect_true(all(c("species", "Variable", "Sign") %in% names(top)))
})

test_that("workflow validates threads and constructs the requested backend", {
  for (threads in list(0, -1, 1.5, NA, "8")) {
    expect_error(runModelingWorkflow(threads = threads), class = "amrml_invalid_argument")
  }
  expect_error(runModelingWorkflow(threads = 8, BPPARAM = BiocParallel::SerialParam()),
               "Supply either", class = "amrml_invalid_argument")
  local_mocked_bindings(.checkArgClass = function(value, name, ...) {
    expect_equal(name, "BPPARAM")
    expect_equal(BiocParallel::bpnworkers(value), 2)
    stop("backend checked")
  })
  expect_error(runModelingWorkflow(threads = 2), "backend checked")
})

test_that("full feature selection warns before writing when struct is available", {
  local_mocked_bindings(
    readORB = function(...) list(directory = tempdir(),
                                 feature_tables = data.frame(feature_type = "struct")),
    .checkArgOutDir = function(...) stop("before writing")
  )
  expect_warning(
    expect_error(runModelingWorkflow(orb_dir = tempdir(), feature_types = "full"),
                 "before writing"),
    "struct models.*very slow"
  )
})

test_that("workflow resume finishes a failed export without refitting", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  metadata <- fxRows(20, 20)
  orb <- fxOrb(metadata = metadata, parquets = list(
    gene_count = fxSignalGenes(20, 20),
    protein_count = dplyr::rename(fxSignalGenes(20, 20), protein = gene)
  ))
  root <- withr::local_tempdir()
  withr::local_dir(root)
  args <- list(orb_dir = orb, out_dir = "run", model_types = "vanilla", seeds = 31L,
               min_genomes = 4, n_fold = 2, save_fits = FALSE, feature_types = "gene",
               BPPARAM = BiocParallel::SerialParam())
  exporter <- exportMergedResults
  local_mocked_bindings(exportMergedResults = function(...) stop("injected export failure"))
  expect_error(suppressWarnings(do.call(runModelingWorkflow, args)), "injected export")
  local_mocked_bindings(exportMergedResults = exporter,
                       .fitModel = function(...) stop("completed models must not refit"),
                       .selectORB = function(...) stop("resume must use saved ORB"),
                       .registerManifest = function(path) path)
  resume_args <- args
  resume_args$orb_dir <- NULL
  expect_error(do.call(runModelingWorkflow, c(resume_args,
    list(resume = TRUE, min_fold_minority = 3))), "saved workflow does not match")
  # Old runs predate the configurable minima and must retain their saved plan.
  snapshot <- .runPaths(normalizePath("run"))$workflow
  legacy <- readRDS(snapshot)
  legacy$settings$min_fold_minority <- NULL
  legacy$settings$min_analysis_minority <- NULL
  legacy$settings$easy_mode <- NULL
  legacy$settings$rescore <- NULL
  # An omitted seed must retain the legacy single-seed plan after the default change.
  resume_args$seeds <- NULL
  .saveFitCheckpoint(legacy, snapshot)
  fit <- do.call(runModelingWorkflow, c(resume_args, list(resume = TRUE)))
  expect_true(all(fit$tasks$status == "fitted"))
  expect_equal(fit$built$out_dir, normalizePath("run"))
  expect_true(all(file.exists(file.path("run", "results", "merged",
                     c("all_perf.parquet", "all_top_features.parquet")))))
  expect_setequal(list.files("run"), c("manifest.json", "metadata", "results"))
  manifest <- jsonlite::read_json(file.path("run", "manifest.json"), simplifyVector = FALSE)
  expect_identical(manifest$manifest_type, "amR_ml")
  matrix_stage <- Filter(
    function(stage) identical(stage$name, "build_matrices"), manifest$runs[[1]]$stages
  )[[1]]
  expect_identical(vapply(matrix_stage$inputs, function(x) basename(x$path), character(1)),
                   "gene_count.parquet")
  expect_identical(vapply(manifest$runs, `[[`, character(1), "status"),
                   c("failed", "success"))
  expect_identical(manifest$artifacts$amRviz_input$status, "ready")
  expect_identical(manifest$artifacts$amRviz_input$producer_run_id,
                   manifest$runs[[2]]$run_id)
  expect_identical(manifest$artifacts$amRviz_input$results$performance,
                   "results/performance.parquet")
  expect_identical(manifest$artifacts$resume_state$status, "removed")
  completion <- tail(manifest$runs[[2]]$events, 1)[[1]]
  expect_equal(completion$details$n_tasks, 2)
  expect_equal(completion$details$n_failed, 0)
  expect_gte(completion$details$elapsed_seconds, 0)
  fit_stage <- Filter(
    function(stage) identical(stage$name, "fit_models"), manifest$runs[[2]]$stages
  )[[1]]
  output_paths <- vapply(fit_stage$outputs, `[[`, character(1), "path")
  expect_true("results/performance.parquet" %in% output_paths)
  expect_true(all(vapply(fit_stage$outputs, `[[`, logical(1), "exists")))
  resumed <- vapply(manifest$runs[[2]]$stages, function(stage) {
    identical(stage$message, "Resumed: reused successful output from a previous run.")
  }, logical(1))
  expect_true(any(resumed))
})

test_that("saved sparse workflows reproduce held-out predictions", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(sparse = TRUE, labels = "real")
  fit <- suppressWarnings(fitModels(
    tasks, save_fits = TRUE, save_state = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  ))
  index <- readRDS(.runPaths(tasks$built$out_dir)$matrix_index)
  for (i in seq_len(nrow(fit$tasks))) {
    id <- fit$tasks$task_id[[i]]
    saved <- readRDS(file.path(tasks$built$out_dir, "results", "models", paste0(id, ".rds")))
    matrix <- .readMatrixWide(file.path(.runPaths(tasks$built$out_dir)$features, "0001.rds"),
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


test_that("tune CLI output is quiet and tuning plus final-fit warnings are retained", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  # Exercise tune's real CLI catalog; tune normally disables it under testthat.
  withr::local_envvar(TESTTHAT = "false")
  withr::local_options(list(cli.progress_handlers_only = "cli", cli.progress_show_after = 0))
  data <- withr::with_seed(17, .modelData(
    matrix(stats::rnorm(26 * 3), nrow = 26), c(rep("Resistant", 6), rep("Susceptible", 20))
  ))
  expect_output(expect_message(expect_warning(
    fitted <- .fitModel(data, n_fold = 3), NA), NA), NA)
  expect_true(any(fitted$notes$stage == "tune"))
  expect_true(any(fitted$notes$stage == "final_fit"))
  expect_true(any(grepl("fewer than 8", fitted$notes$note)))
  expect_true(all(fitted$notes$n >= 1L))
  expect_equal(getOption("cli.progress_handlers_only"), "cli")
})

test_that("task diagnostics survive failure, checkpointing, and completed resume", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  tasks <- fxFitTasks(labels = "real")
  fit_task <- .fitTask
  local_mocked_bindings(.fitTask = function(task, ...) {
    warning("retained task warning")
    if (task$matrix_id == "m2") stop("injected model failure")
    fit_task(task, ...)
  })
  fit <- fitModels(tasks, reuse_fits = FALSE, save_state = TRUE,
                    BPPARAM = BiocParallel::SerialParam())
  expect_equal(sum(fit$tasks$status == "failed"), 1L)
  expect_setequal(fit$notes$task_id[fit$notes$note == "retained task warning"], tasks$tasks$task_id)
  expect_true(file.exists(.runPaths(tasks$built$out_dir)$notes))
  reopened <- fitModels(tasks, reuse_fits = FALSE, save_state = TRUE, resume = TRUE,
                         BPPARAM = BiocParallel::SerialParam())
  expect_equal(reopened$notes, fit$notes)
})


test_that("failed tuning retains its collected notes on the error", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  data <- .modelData(matrix(1, nrow = 20, ncol = 2),
                     rep(c("Resistant", "Susceptible"), 10))
  error <- tryCatch(suppressWarnings(.fitModel(data, n_fold = 2)), error = identity)
  expect_s3_class(error, "error")
  expect_gt(nrow(error$amr_notes), 0L)
  expect_true(any(error$amr_notes$type == "error"))
})
