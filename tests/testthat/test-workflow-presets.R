test_that("easy mode selects only a single-seed vanilla preview before planning", {
  orb <- fxOrb(metadata = fxRows(20, 20), parquets = list(
    gene_count = fxSignalGenes(20, 20),
    protein_count = dplyr::rename(fxSignalGenes(20, 20), protein = gene),
    protein_Pfam_count = dplyr::rename(fxSignalGenes(20, 20), Pfam = gene)
  ))
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  out <- file.path(withr::local_tempdir(), "run")
  expect_message(fit <- runModelingWorkflow(orb, out_dir = out, easy_mode = TRUE,
    min_genomes = 4, n_fold = 2), "Easy mode fits vanilla")
  expect_setequal(unique(fit$tasks$seed), 123L)
  expect_setequal(unique(fit$tasks$labels), "real")
  expect_setequal(unique(fit$built$matrices$mode_id), "vanilla")
  expect_setequal(unique(fit$built$matrices$feature_type), c("gene", "protein", "Pfam"))
  expect_identical(fit$rescoring_status$status, "skipped")
  manifest <- jsonlite::read_json(file.path(out, "manifest.json"))
  expect_identical(manifest$artifacts$feature_rescoring$status, "skipped")
  expect_true(manifest$runs[[1]]$parameters$easy_mode)
  expect_false(dir.exists(file.path(out, "state")))
  for (args in list(list(seeds = 1:3), list(labels = "shuffled"),
                    list(feature_types = "full"), list(model_types = "cross_year"))) {
    expect_error(do.call(runModelingWorkflow, c(list(easy_mode = TRUE), args)),
                 "Easy mode requires", class = "amrml_invalid_argument")
  }
})

fxRescoringOrb <- function(env = parent.frame()) {
  orb <- fxOrb(env = env, metadata = fxRows(20, 20),
               parquets = list(gene_count = fxSignalGenes(20, 20)))
  db <- file.path(orb, "Tst_parquet.duckdb")
  unlink(db) # Replace the fixture's empty placeholder with a valid mapping database.
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db)
  DBI::dbWriteTable(con, "dyad_feature", data.frame(
    dyad = c("orphan|geneA", "orphan|geneB"), feature = c("gene:geneA", "gene:geneB")))
  DBI::dbDisconnect(con, shutdown = TRUE)
  fxRecordSizes(orb)
  orb
}

fxWorkflowArgs <- function(orb, out) {
  list(orb_dir = orb, out_dir = out, model_types = "vanilla", feature_types = "gene",
       min_genomes = 4, n_fold = 2, save_state = TRUE)
}

test_that("normal workflow rescoring exports results and resumes without fitting", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  orb <- fxRescoringOrb()
  args <- fxWorkflowArgs(orb, file.path(withr::local_tempdir(), "run"))
  original_wd <- getwd()
  fit <- do.call(runModelingWorkflow, args)
  expect_equal(getwd(), original_wd)
  expect_setequal(unique(fit$tasks$seed), c(123L, 456L, 789L))
  expect_identical(fit$rescoring_status$status, "success")
  expect_s3_class(fit$rescoring, "feature_dyad_discovery")
  expect_true("geneA" %in% fit$rescoring$features$top$variable)
  paths <- file.path(args$out_dir, "results", "rescoring")
  expect_setequal(list.files(paths), c("optimal_models.parquet", "model_summary.parquet",
    "qualified_features.parquet", "qualified_dyads.parquet", "parameters.json"))
  expect_true(all(file.exists(file.path(args$out_dir, "results",
    c("magic_ready_performance.parquet", "magic_ready_importance.parquet")))))
  expect_false(dir.exists(file.path(args$out_dir, "magic")))
  manifest <- jsonlite::read_json(file.path(args$out_dir, "manifest.json"))
  expect_identical(manifest$artifacts$feature_rescoring$status, "ready")
  local_mocked_bindings(.fitModel = function(...) stop("must not refit"))
  resumed <- do.call(runModelingWorkflow, c(args, list(resume = TRUE)))
  expect_equal(resumed$rescoring, fit$rescoring, ignore_attr = "output_dir")
})

test_that("rescoring failure keeps resumable model results and empty results succeed", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  orb <- fxRescoringOrb()
  args <- fxWorkflowArgs(orb, file.path(withr::local_tempdir(), "run"))
  rescoring <- rescoreFeatures
  local_mocked_bindings(rescoreFeatures = function(...) stop("injected rescoring failure"))
  expect_error(do.call(runModelingWorkflow, args), "injected rescoring failure")
  expect_true(file.exists(file.path(args$out_dir, "state", "workflow.rds")))
  expect_true(file.exists(file.path(args$out_dir, "results", "performance.parquet")))
  local_mocked_bindings(.fitModel = function(...) stop("must not refit"),
    rescoreFeatures = function(...) .noRescoringEvidence("No qualifying features."))
  fit <- do.call(runModelingWorkflow, c(args, list(resume = TRUE)))
  expect_identical(fit$rescoring_status$status, "empty")
  manifest <- jsonlite::read_json(file.path(args$out_dir, "manifest.json"))
  expect_identical(tail(manifest$runs, 1)[[1]]$status, "success")
  expect_identical(manifest$artifacts$feature_rescoring$status, "empty")
  local_mocked_bindings(rescoreFeatures = rescoring)
  expect_s3_class(rescoreFeatures(fit), "feature_dyad_discovery")
})

test_that("rescoring eligibility counts matched successful seeds per representation", {
  listed <- fxFitTasks()
  listed$built$matrices$mode_id <- "vanilla"
  fit <- structure(list(built = listed$built, tasks = tidyr::expand_grid(
    matrix_id = c("m1", "m2"), seed = 1:3, labels = c("real", "shuffled")) |>
    dplyr::mutate(task_id = paste0("t", dplyr::row_number()), status = "fitted")),
    class = "amr_fit")
  fit$tasks$status[fit$tasks$matrix_id == "m2" & fit$tasks$seed == 3 &
                     fit$tasks$labels == "shuffled"] <- "failed"
  expect_equal(.rescoreScopes(fit, 3L)$matrix_id, "m1")
  fit$tasks$status[fit$tasks$matrix_id == "m1" & fit$tasks$seed == 2 &
                     fit$tasks$labels == "real"] <- "failed"
  expect_equal(nrow(.rescoreScopes(fit, 3L)), 0L)
  expect_error(rescoreFeatures(fit), "at least 3", class = "amrml_invalid_argument")
})

test_that("dominant features cross the cumulative cutoff and seed counts precede filtering", {
  p <- withr::local_tempfile(fileext = ".parquet")
  t <- withr::local_tempfile(fileext = ".parquet")
  perf <- tidyr::expand_grid(seed = 1:3, shuffled = c(FALSE, TRUE)) |>
    dplyr::mutate(species = "Tst", drug_label = "drug", drug_or_class = "AMP",
      feature_type = "gene", feature_subtype = "binary", model = "LR",
      fit_penalty = 0.01, fit_mixture = 1, mcc = ifelse(shuffled, 0.1, 0.9),
      n_feat = 2L, n_feats_returned = 2L)
  top <- tidyr::expand_grid(seed = 1:3, Variable = c("dominant", "minor")) |>
    dplyr::mutate(shuffled = FALSE, species = "Tst", drug_label = "drug",
      drug_or_class = "AMP", feature_type = "gene", feature_subtype = "binary",
      Importance = ifelse(Variable == "dominant", 0.9, 0.1), Sign = "NEG")
  arrow::write_parquet(perf, p)
  arrow::write_parquet(top, t)
  result <- topFeaturesPerDrugOrClass(t, all_performance_parquet = p)
  expect_equal(result$variable, "dominant")
  expect_equal(result$seed_ratio, 1)
  perf$mcc[!perf$shuffled & perf$seed == 3] <- 0
  arrow::write_parquet(perf, p)
  result <- topFeaturesPerDrugOrClass(t, all_performance_parquet = p)
  expect_equal(result$seed_ratio, 2/3)
  expect_equal(nrow(topFeaturesPerDrugOrClass(t, all_performance_parquet = p,
                                             seed_ratio_threshold = 1)), 0L)
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  expect_error(runFeatureDyadDiscovery(t, p, con, MCC_threshold = 1),
               "No model fits passed", class = "amrml_no_rescoring_evidence")
  expect_error(runFeatureDyadDiscovery(t, p, con),
               "No features passed", class = "amrml_no_rescoring_evidence")
})

test_that("magic-ready tables keep every original column and add the per-seed calculations", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  args <- fxWorkflowArgs(fxRescoringOrb(), file.path(withr::local_tempdir(), "run"))
  fit <- do.call(runModelingWorkflow, args)
  paths <- .runPaths(fit$built$out_dir)
  merged <- file.path(paths$merged, c("all_perf.parquet", "all_top_features.parquet"))
  ready <- c(paths$magic_ready_performance, paths$magic_ready_importance)

  perf <- arrow::read_parquet(ready[[1]])
  top <- arrow::read_parquet(ready[[2]])
  original <- lapply(merged, arrow::read_parquet)
  expect_true(all(names(original[[1]]) %in% names(perf)))
  expect_true(all(names(original[[2]]) %in% names(top)))
  expect_equal(nrow(perf), nrow(original[[1]]))
  expect_equal(nrow(top), nrow(original[[2]]))
  expect_true(all(c("shuffled_mcc", "mcc_diff", "above_background", "feat_return_ratio") %in% names(perf)))
  expect_true(all(c("variable_clean", "contribution", "rank_score", "cum_contrib") %in% names(top)))
  manifest <- jsonlite::read_json(file.path(args$out_dir, "manifest.json"))
  expect_identical(manifest$artifacts$feature_rescoring$status, "ready")
  expect_identical(manifest$artifacts$feature_rescoring$magic_ready$performance,
                   "results/magic_ready_performance.parquet")
  expect_identical(manifest$artifacts$feature_rescoring$magic_ready$importance,
                   "results/magic_ready_importance.parquet")
  expect_identical(manifest$artifacts$amRviz_input$results$magic_ready_performance,
                   "results/magic_ready_performance.parquet")
  expect_identical(manifest$artifacts$amRviz_input$results$magic_ready_importance,
                   "results/magic_ready_importance.parquet")
})

test_that("a model with no importance does not stop the magic-ready export", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  args <- fxWorkflowArgs(fxRescoringOrb(), file.path(withr::local_tempdir(), "run"))
  fit <- do.call(runModelingWorkflow, args)
  paths <- .runPaths(fit$built$out_dir)
  top <- arrow::read_parquet(file.path(paths$merged, "all_top_features.parquet"))
  seed <- top$seed[!top$shuffled][[1]]
  zero <- !top$shuffled & top$seed == seed
  top$Importance[zero] <- 0
  arrow::write_parquet(top, file.path(paths$merged, "all_top_features.parquet"))

  ready <- .completeRescoring(file.path(paths$merged, "all_perf.parquet"),
                              file.path(paths$merged, "all_top_features.parquet"))
  expect_equal(nrow(ready$all_top_features), nrow(top))
  expect_true(all(is.na(ready$all_top_features$contribution[zero])))
})

test_that("magic-ready tables are built once and rebuilt only when the models are refitted", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  args <- fxWorkflowArgs(fxRescoringOrb(), file.path(withr::local_tempdir(), "run"))
  fit <- do.call(runModelingWorkflow, args)
  paths <- .runPaths(fit$built$out_dir)
  ready <- c(paths$magic_ready_performance, paths$magic_ready_importance)
  first <- file.mtime(ready)
  expect_true(.magicReadyCurrent(fit$built$out_dir))

  # Rescoring again, with a different threshold, reuses them.
  Sys.sleep(1.1)
  rescoreFeatures(fit, MCC_threshold = 0.4)
  expect_equal(file.mtime(ready), first)

  # Refitted models write newer results, so the tables are rebuilt.
  Sys.sleep(1.1)
  Sys.setFileTime(paths$performance, Sys.time())
  expect_false(.magicReadyCurrent(fit$built$out_dir))
  rescoreFeatures(fit, MCC_threshold = 0.4)
  expect_true(all(file.mtime(ready) > first))
})

test_that("the workflow's rescoring minimum is the one the user passes", {
  listed <- fxFitTasks()
  listed$built$matrices$mode_id <- "vanilla"
  fit <- structure(list(built = listed$built, tasks = tidyr::expand_grid(
    matrix_id = "m1", seed = 1:3, labels = c("real", "shuffled")) |>
    dplyr::mutate(task_id = paste0("t", dplyr::row_number()), status = "fitted")),
    class = "amr_fit")
  manifest <- withr::local_tempfile(fileext = ".json")
  local_mocked_bindings(.manifestStage = function(...) NULL, .manifestArtifact = function(...) NULL,
                        .manifestEvent = function(...) NULL)
  out <- suppressMessages(.rescoreWorkflow(fit, manifest, min_seeds = 4L))
  expect_identical(out$rescoring_status$status, "skipped")
  expect_match(out$rescoring_status$reason, "at least 4")
  expect_error(
    runModelingWorkflow(rescore_min_seeds = 0), class = "amrml_invalid_argument"
  )
})
