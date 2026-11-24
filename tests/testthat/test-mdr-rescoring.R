test_that("MAGIC eligibility excludes both MDR targets without changing drug/class pairs", {
  listed <- fxFitTasks()
  matrices <- listed$built$matrices
  drug <- matrices[1, ]
  drug$mode_id <- "vanilla"
  drug$unit <- "drug"
  drug$target <- "DRA"
  drug_class <- drug
  drug_class$matrix_id <- "class"
  drug_class$unit <- "drug_class"
  drug_class$target <- "CLX"
  mdr <- drug
  mdr$matrix_id <- "mdr"
  mdr$unit <- "mdr"
  mdr$target <- "mdr"
  combination <- mdr
  combination$matrix_id <- "combination"
  combination$target <- "resistant_classes"
  for (mode in c("mdr", "vanilla")) {
    mdr$mode_id <- combination$mode_id <- mode
    listed$built$matrices <- dplyr::bind_rows(drug, drug_class, mdr, combination)
    tasks <- tidyr::expand_grid(matrix_id = listed$built$matrices$matrix_id,
      seed = 1:3, labels = c("real", "shuffled")) |>
      dplyr::mutate(task_id = paste(matrix_id, seed, labels), status = "fitted")
    fit <- structure(list(built = listed$built, tasks = tasks), class = "amr_fit")
    scopes <- .rescoreScopes(fit, min_seeds = 3L)
    expect_setequal(scopes$matrix_id, c(drug$matrix_id, "class"))
    expect_equal(scopes$n_paired_seeds, c(3L, 3L))
    fit$tasks <- tasks[!tasks$matrix_id %in% c("mdr", "combination"), ]
    expect_equal(.rescoreScopes(fit, min_seeds = 3L), scopes)
  }
})
