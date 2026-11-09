test_that("MAGIC never combines endpoints from different MDR combination labels", {
  root <- withr::local_tempdir()
  run <- file.path(root, "run")
  orb <- file.path(root, "dataset", "orb")
  dir.create(file.path(run, "results"), recursive = TRUE)
  dir.create(file.path(run, "metadata"))
  dir.create(orb, recursive = TRUE)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = file.path(orb, "test_parquet.duckdb"))
  DBI::dbWriteTable(con, "dyads", data.frame(
    gene = c("g1", "g2"), protein = c("p1", "p2"), dyad_type = "paired"))
  DBI::dbWriteTable(con, "protein_members", data.frame(member = c("p1", "p2"),
    cluster = c("p1", "p2")))
  for (view in c("dyad_ann_pfam", "dyad_ann_cog", "dyad_ann_amr", "dyad_ann_defense")) {
    DBI::dbWriteTable(con, view, data.frame(cluster = character(), feature = character()))
  }
  DBI::dbDisconnect(con, shutdown = TRUE)
  matrices <- tibble::tibble(matrix_id = c("g", "p"), mode_id = "vanilla",
    unit = "mdr", target = "resistant_classes", feature_type = c("gene", "protein"),
    encoding = "binary")
  tasks <- tidyr::expand_grid(matrix_id = c("g", "p"), seed = 1:3,
    labels = c("real", "shuffled")) |>
    dplyr::mutate(task_id = paste(matrix_id, seed, labels, sep = ":"), status = "fitted")
  performance <- tasks |>
    dplyr::transmute(task_id, seed, mcc = ifelse(labels == "real", .8, -.2))
  importance <- tidyr::crossing(tasks, label = c("A", "B"), feature_number = 1:2) |>
    dplyr::mutate(Variable = paste0(matrix_id, feature_number),
      Importance = ifelse(
        (matrix_id == "g" & label == "A" | matrix_id == "p" & label == "B") ==
          (feature_number == 1), 9, 1),
      Sign = "POS") |>
    dplyr::select(task_id, label, Variable, Importance, Sign)
  predictions <- tidyr::crossing(task_id = tasks$task_id, truth = c("A", "B")) |>
    dplyr::mutate(predicted = truth)
  paths <- .runPaths(run)
  .writeParquet(matrices, paths$matrices)
  .writeParquet(tasks, paths$tasks)
  .writeParquet(performance, paths$performance)
  .writeParquet(importance, paths$importance)
  .writeParquet(predictions, paths$predictions)
  out <- file.path(root, "mismatched")
  expect_no_warning(exportNetwork(run, orb, out))
  graph <- arrow::read_parquet(file.path(out, "cluster_pair_network_evidence.parquet"))
  expect_setequal(graph$label, c("A", "B"))
  expect_equal(nrow(graph), 8L)
  expect_false(any(graph$high_priority))
  expect_false(any(graph$pair_robust75))
  # Align each protein with the gene supporting the SAME label.
  importance$Importance[grepl("^p", importance$Variable)] <-
    10 - importance$Importance[grepl("^p", importance$Variable)]
  .writeParquet(importance, paths$importance)
  out <- file.path(root, "matched")
  expect_no_warning(exportNetwork(run, orb, out))
  graph <- arrow::read_parquet(file.path(out, "cluster_pair_network_evidence.parquet"))
  priority <- graph[graph$high_priority, ]
  expect_equal(nrow(priority), 2L)
  expect_equal(priority$gene_cluster[priority$label == "A"], "g1")
  expect_equal(priority$protein_cluster[priority$label == "A"], "p1")
  expect_equal(priority$gene_cluster[priority$label == "B"], "g2")
  # Model reliability stays at 2 representations, not 4 label views.
  reliability <- arrow::read_parquet(file.path(out, "model_reliability.parquet"))
  expect_equal(nrow(reliability), 2L)
})
