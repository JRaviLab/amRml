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
  expect_equal(readRDS(file.path(paths, "discovery.rds")), fit$rescoring)
  expect_true(file.exists(file.path(paths, "features_top.parquet")))
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
  expect_equal(.rescoreScopes(fit)$matrix_id, "m1")
  fit$tasks$status[fit$tasks$matrix_id == "m1" & fit$tasks$seed == 2 &
                     fit$tasks$labels == "real"] <- "failed"
  expect_equal(nrow(.rescoreScopes(fit)), 0L)
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
