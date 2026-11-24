test_that("MAGIC network export excludes new and legacy MDR targets", {
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
  matrices <- tibble::tibble(matrix_id = c("g", "p"), mode_id = "mdr",
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
  for (mode in c("mdr", "vanilla")) {
    matrices$mode_id <- mode
    for (target in c("mdr", "resistant_classes")) {
      matrices$target <- target
      .writeParquet(matrices, paths$matrices)
      expect_error(exportNetwork(run, orb, file.path(root, paste(mode, target))),
                   "Expected distinct fitted vanilla tasks")
    }
  }
})
