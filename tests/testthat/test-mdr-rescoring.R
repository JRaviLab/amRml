# Label-specific MDR evidence must not be pooled or treated as independent fits.
fxMdrDiscovery <- function(env = parent.frame()) {
  perf_path <- withr::local_tempfile(fileext = ".parquet", .local_envir = env)
  top_path <- withr::local_tempfile(fileext = ".parquet", .local_envir = env)
  perf <- tidyr::expand_grid(drug_label = c("mdr", "drug"), seed = 1:3,
      shuffled = c(FALSE, TRUE)) |>
    dplyr::mutate(species = "Tst",
      drug_or_class = ifelse(drug_label == "mdr", "resistant_classes", "AMP"),
      feature_type = "gene", feature_subtype = "binary", model = "LR",
      fit_penalty = .01, fit_mixture = 1, mcc = ifelse(shuffled, .1, .9),
      n_feat = 3L, n_feats_returned = 2L)
  mdr <- tidyr::expand_grid(seed = 1:3, label = c("CLX", "CLX_CLY", "Susceptible"),
      minor = c(FALSE, TRUE)) |>
    dplyr::mutate(species = "Tst", drug_label = "mdr",
      drug_or_class = "resistant_classes", shuffled = FALSE,
      feature_type = "gene", feature_subtype = "binary",
      Variable = ifelse(minor, "minor", ifelse(label == "Susceptible", "sus", "shared")),
      Importance = ifelse(minor, 1, 9), Sign = "POS") |>
    dplyr::select(-minor)
  drug <- tidyr::expand_grid(seed = 1:3, minor = c(FALSE, TRUE)) |>
    dplyr::mutate(species = "Tst", drug_label = "drug", drug_or_class = "AMP",
      label = NA_character_, shuffled = FALSE, feature_type = "gene",
      feature_subtype = "binary", Variable = ifelse(minor, "minor", "shared"),
      Importance = ifelse(minor, 1, 9), Sign = "NEG") |>
    dplyr::select(-minor)
  arrow::write_parquet(perf, perf_path)
  arrow::write_parquet(dplyr::bind_rows(mdr, drug), top_path)
  list(perf = perf_path, top = top_path)
}

test_that("MDR coefficients normalize and recur independently within each label", {
  paths <- fxMdrDiscovery()
  scored <- scoreFeaturesWithinSeed(paths$top, all_performance_parquet = paths$perf)
  sums <- scored |>
    dplyr::group_by(drug_label, label, seed) |>
    dplyr::summarise(total = sum(contribution), .groups = "drop")
  expect_equal(sums$total, rep(1, 12))
  summary <- summariseFeaturesAcrossSeeds(scored,
    arrow::read_parquet(paths$perf) |> dplyr::filter(!shuffled))
  expect_true(all(summary$seed_ratio == 1))
  expect_equal(nrow(summary[summary$drug_label == "mdr" & summary$variable == "shared", ]), 2)
})

test_that("discovery distinguishes shared MDR outcomes from shared MDR and drug evidence", {
  paths <- fxMdrDiscovery()
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(con, "dyad_feature", data.frame(
    dyad = c("p|g", "p2|g2"), feature = c("gene:shared", "gene:sus")))
  result <- runFeatureDyadDiscovery(paths$top, paths$perf, con, rank_score_quantile = 0)
  top <- result$features$top
  expect_equal(nrow(top[top$drug_label == "mdr", ]), 3)
  expect_equal(result$models$summary$n_passing_fits, c(3L, 3L))
  shared <- result$features$mdr_shared
  expect_equal(shared$variable, "shared")
  expect_true(shared$shared_within_mdr)
  expect_true(shared$shared_mdr_single_drug)
  expect_equal(shared$n_mdr_combination_labels, 2L)
  expect_equal(shared$n_model_targets, 2L)
  expect_equal(shared$n_outcomes, 3L)
  expect_equal(result$dyads$mdr_shared$dyad, "p|g")
  expect_equal(nrow(result$features$overlap), 6)
  expect_true(any(result$features$overlap$jaccard == 0))
  expect_true(any(result$features$overlap$jaccard == 1))
  # Binary/counts copies do not count as additional outcomes or memberships.
  duplicate <- dplyr::bind_rows(top, dplyr::mutate(top, feature_subtype = "counts"))
  comparison <- .compareOutcomeEvidence(duplicate, c("feature_type", "variable"))
  expect_equal(comparison$cross_model_sharedness, result$features$cross_model_sharedness)
})

test_that("network normalization and endpoint evidence retain MDR label identity", {
  metadata <- tibble::tibble(task_id = "t1", matrix_id = "m1", unit = "mdr",
    target = "resistant_classes", labels = "real", seed = 1L, encoding = "binary",
    feature_type = "gene")
  importance <- tibble::tibble(task_id = "t1", label = c("A", "B"),
    Variable = "g1", Importance = c(9, 1), Sign = c("POS", "NEG"))
  scored <- .scoreNetworkTasks(importance, metadata, .75)
  expect_equal(scored$contribution, c(1, 1))
  pairs <- tibble::tibble(gene_cluster = "g1", protein_cluster = "p1")
  nodes <- .prepareNetworkNodes(scored, pairs)
  expect_equal(nodes$label, c("A", "B"))
  fitted <- dplyr::bind_rows(
    dplyr::select(nodes, unit, target, label, labels, encoding, feature_type, seed),
    dplyr::mutate(dplyr::select(nodes, unit, target, label, labels, encoding, feature_type, seed),
      seed = 2L))
  summary <- .summarizeNetworkNodes(nodes, fitted)
  expect_equal(summary$n_fitted_seeds, c(2L, 2L))
  expect_equal(summary$n_returned_seeds, c(1L, 1L))
})

test_that("network outcome inventories retain labels with no returned coefficients", {
  metadata <- tibble::tibble(task_id = c("g1", "g2", "b1"),
    matrix_id = c("g", "g", "b"), unit = c("mdr", "mdr", "drug"),
    target = c("resistant_classes", "resistant_classes", "AMP"))
  predictions <- tibble::tibble(task_id = c("g1", "g1", "g2", "g2", "b1"),
    truth = c("A", "B", "A", "B", "Resistant"))
  labels <- .networkOutcomeLabels(metadata, predictions)
  expect_setequal(labels$label[labels$matrix_id == "g"], c("A", "B"))
  expect_true(is.na(labels$label[labels$matrix_id == "b"]))
})

test_that("one-versus-rest label diagnostics preserve model identity and undefined values", {
  listed <- fxFitTasks()
  listed$built$matrices$unit <- "mdr"
  listed$built$matrices$target <- "resistant_classes"
  fit <- structure(list(built = listed$built, tasks = listed$tasks,
    predictions = tibble::tibble(task_id = listed$tasks$task_id[[1L]],
      species = "Tst", truth = c("A", "A", "B", "B"), predicted = "A")),
    class = "amr_fit")
  fit$tasks$status <- "fitted"
  diag <- .mdrLabelDiagnostics(fit)
  expect_equal(nrow(diag), 2L)
  expect_equal(diag$n_truth, c(2L, 2L))
  expect_equal(diag$n_predicted, c(4L, 0L))
  expect_true(all(is.na(diag$label_mcc)))
  expect_equal(diag$sensitivity, c(1, 0))
  expect_true(is.na(diag$precision[diag$label == "B"]))
  expect_equal(dplyr::n_distinct(diag$task_id), 1L)
})

test_that("candidate-aware scoring keeps MDR labels separate and accepts old binary inventories", {
  paths <- fxMdrDiscovery()
  scored <- scoreFeaturesWithinSeed(paths$top, all_performance_parquet = paths$perf)
  keys <- c("species", "drug_label", "drug_or_class", "feature_type",
    "feature_subtype", "seed", "label", "variable")
  inventory <- scored[keys]
  result <- summarizeFeatureEvidence(scored, inventory)
  shared <- result[result$drug_label == "mdr" & result$variable == "shared", ]
  expect_equal(nrow(shared), 2L)
  expect_equal(shared$n_eligible_seeds, c(3L, 3L))
  expect_equal(shared$direct_score, c(.9, .9))
  binary <- scored[scored$drug_label == "drug", ]
  old_inventory <- binary[setdiff(keys, "label")]
  expect_equal(nrow(summarizeFeatureEvidence(binary, old_inventory)), 2L)
  expect_error(summarizeFeatureEvidence(scored, inventory[setdiff(keys, "label")]),
    "Expected per-seed")
})

test_that("MDR diagnostics include false positives for prediction-only labels", {
  listed <- fxFitTasks()
  listed$built$matrices$unit <- "mdr"
  listed$built$matrices$target <- "resistant_classes"
  fit <- structure(list(built = listed$built, tasks = listed$tasks,
    predictions = tibble::tibble(task_id = listed$tasks$task_id[[1L]],
      species = "Tst", truth = c("A", "A", "B", "B"),
      predicted = c("A", "C", "C", "C"))), class = "amr_fit")
  fit$tasks$status <- "fitted"
  diag <- .mdrLabelDiagnostics(fit)
  expect_setequal(diag$label, c("A", "B", "C"))
  expect_equal(sum(diag$n_truth), 4L)
  expect_equal(sum(diag$n_predicted), 4L)
  missing <- diag[diag$label == "C", ]
  expect_equal(missing$n_truth, 0L)
  expect_equal(missing$n_predicted, 3L)
  expect_equal(missing$fp, 3L)
  expect_equal(missing$tn, 1L)
  expect_equal(missing$precision, 0)
  expect_true(is.na(missing$sensitivity))
  expect_true(is.na(missing$label_mcc))
})
