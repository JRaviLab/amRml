test_that("per-group diagnostics are verbose-gated while progress and file logs remain", {
  local_mocked_bindings(.runTaskGroup = function(job, ...) {
    .saveFitCheckpoint(rep(list(list(error = "fixture")), nrow(job$tasks)), job$checkpoint)
    TRUE
  })
  for (verbose in c(FALSE, TRUE)) {
    tasks <- fxFitTasks(labels = "real")
    messages <- character()
    fit <- withCallingHandlers(
      fitModels(tasks, BPPARAM = BiocParallel::SerialParam(), verbose = verbose),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    expect_identical(any(grepl("^Running group-", messages)), verbose)
    expect_identical(any(grepl("^Diagnostics:", messages)), verbose)
    expect_true(any(grepl("models complete", messages)))
    expect_true(any(grepl("^Model fitting complete:", messages)))
    logs <- list.files(file.path(tasks$built$out_dir, "diagnostics"), full.names = TRUE)
    expect_true(any(grepl("group_start", unlist(lapply(logs, readLines)))))
  }
})

test_that("prepared arrays preserve verbosity and permit a quiet worker override", {
  tasks <- fxFitTasks(labels = "real")
  plan <- fitModels(tasks, execution = "prepare", save_state = TRUE, verbose = TRUE)
  expect_true(plan$verbose)
  local_mocked_bindings(.runTaskGroup = function(job, ...) {
    expect_false(job$verbose)
    .saveFitCheckpoint(rep(list(list(error = "fixture")), nrow(job$tasks)), job$checkpoint)
    TRUE
  })
  messages <- character()
  withCallingHandlers(runTaskChunk(tasks$built$out_dir, 1, verbose = FALSE),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_false(any(grepl("^Running |^Array worker:|^Array chunk ", messages)))
  expect_true(any(grepl("models complete", messages)))
})
