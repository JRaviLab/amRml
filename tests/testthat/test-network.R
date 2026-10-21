testthat::test_that("MAGIC network cores include the crossing feature and conserve importance", {
  importance <- tibble::tibble(
    task_id = rep("fit1", 4),
    Variable = c("g3", "g1", "g4", "g2"),
    Importance = c(1, 4, 0.5, 2.5),
    Sign = "NEG"
  )
  metadata <- tibble::tibble(
    task_id = "fit1", unit = "drug", target = "AMP",
    labels = "real", seed = 31L, feature_type = "gene",
    encoding = "binary"
  )
  scored <- .scoreNetworkTasks(importance, metadata, 0.75)
  testthat::expect_equal(sum(scored$contribution), 1)
  testthat::expect_equal(scored$Variable, c("g1", "g2", "g3", "g4"))
  testthat::expect_equal(scored$Variable[scored$in_core], c("g1", "g2"))
  testthat::expect_true(sum(scored$contribution[scored$in_core]) >= 0.75)
  testthat::expect_error(.scoreNetworkTasks(
    dplyr::bind_rows(importance, importance[1, ]), metadata, 0.75),
    "Invalid returned")
})

testthat::test_that("protein predictors already represent CD-HIT clusters", {
  scored <- tibble::tibble(
    task_id = rep("fit1", 3), unit = "drug", target = "AMP",
    labels = "real", seed = 31L, feature_type = "protein", encoding = "binary",
    Variable = c("fig.1.peg.1", "clusterB", "member1"),
    contribution = c(.4, .35, .25), in_core = c(TRUE, TRUE, FALSE),
    Sign = c("NEG", "NEG", "POS")
  )
  pairs <- tibble::tibble(gene_cluster = "g1",
                          protein_cluster = "fig|1.peg.1")
  result <- .prepareNetworkNodes(scored, pairs)
  testthat::expect_equal(nrow(result), 3L)
  testthat::expect_equal(result$node_id, c("fig|1.peg.1", "clusterB", "member1"))
  testthat::expect_equal(result$contribution, c(.4, .35, .25))
  testthat::expect_equal(result$in_core, c(TRUE, TRUE, FALSE))
  testthat::expect_equal(result$in_orb_graph, c(TRUE, FALSE, FALSE))
  testthat::expect_true(all(result$n_returned_features == 1L))

  duplicate <- dplyr::bind_rows(scored[1, ], scored[1, ])
  duplicate$Variable[[2L]] <- "fig|1.peg.1"
  testthat::expect_error(.prepareNetworkNodes(duplicate, pairs),
                         "same cluster ID")
})

testthat::test_that("network recurrence is separate from direct feature scores", {
  rows <- tibble::tibble(
    task_id = c("fit1", "fit2"), unit = "drug", target = "AMP",
    labels = "real", seed = c(31L, 79L), encoding = "binary",
    feature_type = "gene", node_id = "g1", in_core = c(TRUE, TRUE),
    contribution = c(.4, .2), n_returned_features = 1L
  )
  fitted <- tibble::tibble(
    unit = "drug", target = "AMP", labels = "real",
    encoding = "binary", feature_type = "gene", seed = c(31L, 79L, 51L)
  )
  nodes <- .summarizeNetworkNodes(rows, fitted)
  testthat::expect_identical(nodes$n_returned_seeds, 2L)
  testthat::expect_identical(nodes$n_fitted_seeds, 3L)
  testthat::expect_identical(nodes$n_core_seeds, 2L)
  testthat::expect_true(nodes$reproducibly_important)
  testthat::expect_equal(nodes$mean_returned_contribution, .3)
  testthat::expect_false("direct_score" %in% names(nodes))
})

testthat::test_that("task-based ORB network export preserves unsupported and shuffled evidence", {
  testthat::skip_if_not_installed("arrow")
  testthat::skip_if_not_installed("duckdb")
  root <- tempfile("magic-network-test-")
  run <- file.path(root, "run")
  orb <- file.path(root, "dataset", "orb")
  out <- file.path(root, "network")
  dir.create(file.path(run, "results"), recursive = TRUE)
  dir.create(orb, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = file.path(orb, "fixture_parquet.duckdb"))
  DBI::dbWriteTable(con, "dyads", data.frame(gene = "g1", protein = "fig|1.peg.1"))
  DBI::dbWriteTable(con, "protein_members", data.frame(member = "member1", cluster = "fig|1.peg.1"))
  DBI::dbWriteTable(con, "dyad_ann_pfam", data.frame(cluster = "fig|1.peg.1", feature = "PF.1"))
  for (view in c("dyad_ann_cog", "dyad_ann_amr", "dyad_ann_defense")) {
    DBI::dbWriteTable(con, view,
      data.frame(cluster = character(), feature = character()))
  }
  DBI::dbDisconnect(con)

  matrices <- tibble::tibble(
    matrix_id = c("g", "p", "a"), mode_id = "vanilla", unit = "drug",
    target = "AMP", feature_type = c("gene", "protein", "Pfam"),
    encoding = "binary"
  )
  tasks <- tidyr::crossing(matrix_id = matrices$matrix_id,
                           labels = c("real", "shuffled"), seed = c(31L, 79L)) |>
    dplyr::mutate(task_id = paste(matrix_id, labels, seed, sep = "__"),
                  status = "fitted")
  performance <- tasks |>
    dplyr::transmute(task_id, seed,
      mcc = dplyr::if_else(labels == "real", 0.8, -0.2))
  importance <- tasks |>
    dplyr::left_join(matrices, by = "matrix_id") |>
    dplyr::transmute(task_id,
      Variable = dplyr::case_when(feature_type == "gene" ~ "g1",
                                   feature_type == "protein" ~ "fig.1.peg.1",
                                   TRUE ~ "PF.1"),
      Importance = 1, Sign = "NEG")
  predictions <- tidyr::crossing(task_id = tasks$task_id,
                                 genome_id = c("x", "y")) |>
    dplyr::mutate(truth = dplyr::if_else(genome_id == "x", "Resistant", "Susceptible"),
                  predicted = truth)
  arrow::write_parquet(matrices, file.path(run, "matrices.parquet"))
  arrow::write_parquet(tasks, file.path(run, "tasks.parquet"))
  arrow::write_parquet(performance, file.path(run, "results/performance.parquet"))
  arrow::write_parquet(importance, file.path(run, "results/importance.parquet"))
  arrow::write_parquet(predictions, file.path(run, "results/predictions.parquet"))

  testthat::expect_no_warning(result <- exportNetwork(run, orb, out))
  testthat::expect_equal(result$summary$n_distinct_high_priority_pairs, 1)
  provenance <- arrow::read_parquet(file.path(out, "network_provenance.parquet"))
  testthat::expect_equal(provenance$n_cluster_pairs, 1)
  testthat::expect_equal(provenance$core_threshold, 0.75)
  ids <- arrow::read_parquet(file.path(out, "protein_cluster_id_crosswalk.parquet"))
  testthat::expect_equal(ids$model_feature, "fig.1.peg.1")
  testthat::expect_equal(ids$protein_cluster, "fig|1.peg.1")
  testthat::expect_true(ids$in_orb_graph)
  cluster_rows <- arrow::read_parquet(file.path(out, "gene_protein_seed_evidence.parquet"))
  testthat::expect_true(all(cluster_rows$n_returned_features == 1L))
  testthat::expect_true(all(cluster_rows$node_id[cluster_rows$feature_type == "protein"] == "fig|1.peg.1"))
  testthat::expect_false(file.exists(file.path(out, "gene_protein_projected_seed_evidence.parquet")))
  graph <- arrow::read_parquet(file.path(out, "cluster_pair_network_evidence.parquet"))
  testthat::expect_equal(nrow(graph), 2)
  testthat::expect_true(graph$high_priority[graph$labels == "real"])
  testthat::expect_false(graph$high_priority[graph$labels == "shuffled"])
  testthat::expect_true(all(graph$pair_robust75))
  testthat::expect_true(all(graph$protein_evidence_origin == "direct_protein_cluster_model"))
  neighbors <- arrow::read_parquet(file.path(out, "cluster_annotation_neighbor_evidence.parquet"))
  testthat::expect_equal(nrow(neighbors), 2)
  testthat::expect_true(all(neighbors$returned_nonzero))
  testthat::expect_equal(result$summary$n_high_priority_pair_encoding_records, 1L)
  testthat::expect_equal(result$summary$n_distinct_real_robust_pairs, 1L)
  testthat::expect_setequal(list.files(out), c(
    "network_provenance.parquet", "feature_returned_seed_evidence.parquet",
    "gene_protein_cluster_node_evidence.parquet", "annotation_node_evidence.parquet",
    "cluster_pair_network_evidence.parquet", "cluster_annotation_neighbor_evidence.parquet",
    "model_reliability.parquet", "reliability_mcc_audit.parquet",
    "protein_cluster_id_crosswalk.parquet", "gene_protein_seed_evidence.parquet"
  ))
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = file.path(orb, "fixture_parquet.duckdb"))
  DBI::dbExecute(con, "DELETE FROM dyads")
  DBI::dbWriteTable(con, "dyads", data.frame(
    gene = c("g1", "g2"), protein = c("member1", "fig|1.peg.2")
  ), append = TRUE)
  DBI::dbDisconnect(con)
  testthat::expect_error(exportNetwork(run, orb, file.path(root, "legacy")),
                         "individual protein IDs")
})
