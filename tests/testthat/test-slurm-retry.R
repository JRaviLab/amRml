# Exercise scheduler recovery with controlled command replies.
slurmHelper <- function() {
  ns <- asNamespace("amRml")
  e <- new.env(parent = ns)
  for (name in grep("^([.]slurm|[.]runHPC$)", ls(ns, all.names = TRUE), value = TRUE)) {
    f <- get(name, ns)
    if (!is.primitive(f)) environment(f) <- e
    assign(name, f, e)
  }
  e$.slurmPause <- function(...) NULL
  e
}

slurmState <- function(root, memory = c(3700L, 10000L, 20000L)) {
  dir.create(file.path(root, "diagnostics"), recursive = TRUE)
  list(run_dir = root, runner = "/runner with spaces.sh", memory_mb = memory,
       cpus_per_task = 1L, array_limit = 2L, worker_args = "--account=lab",
       waves = list(list(job_id = "123", chunks = 1:2, memory_mb = memory[[1L]],
                         attempt = 0L, status = "submitted")))
}

slurmReplies <- function(e, replies) {
  e$calls <- list()
  e$.slurmCommand <- function(command, args) {
    index <- length(e$calls) + 1L
    e$calls[[index]] <- list(command = command, args = args)
    if (index > length(replies)) stop("Unexpected extra scheduler command")
    replies[[index]]
  }
}

test_that("R coordinator distinguishes OOM states and refuses incomplete accounting", {
  e <- slurmHelper()
  rows <- c("123_1|COMPLETED|0:0|", "123_2|FAILED|1:0|",
            "123_2.batch|OUT_OF_MEMORY|0:9|", "999_1|OUT_OF_MEMORY|0:9|")
  checked <- e$.slurmAccounting(rows, "123", 1:2)
  expect_identical(checked$oom, 2L)
  expect_identical(checked$states, c(`1` = "COMPLETED", `2` = "FAILED"))
  expect_length(e$.slurmAccounting("123_1|FAILED|137:0|", "123", 1L)$oom, 0L)
  for (rows in list(character(), "123_1|RUNNING|", "123_1|COMPLETED|")) {
    expect_error(e$.slurmAccounting(rows, "123", 1:2), "incomplete")
  }
})

test_that("R monitor uses exact memory tiers and removes only the terminated attempt's locks", {
  e <- slurmHelper()
  root <- withr::local_tempdir()
  state <- slurmState(root)
  locks <- file.path(root, "state", "fit-checkpoints")
  owners <- list(
    `chunk-000002.lock` = list(array_job_id = "123", array_task_id = "2"),
    `group-000002.rds.lock` = list(array_job_id = "123", array_task_id = "2"),
    `chunk-000001.lock` = list(array_job_id = "123", array_task_id = "1"),
    `group-000003.rds.lock` = list(array_job_id = "999", array_task_id = "2"),
    `group-000004.rds.lock` = list(array_job_id = "123"),
    `group-000005.rds.lock` = NULL,
    `group-000006.rds.lock` = "malformed JSON",
    `group-000007.rds.lock` = "scalar JSON"
  )
  for (name in names(owners)) {
    lock <- file.path(locks, name)
    dir.create(lock, recursive = TRUE)
    if (identical(owners[[name]], "malformed JSON")) {
      cat("{", file = file.path(lock, "owner.json"))
    } else if (!is.null(owners[[name]])) {
      jsonlite::write_json(owners[[name]], file.path(lock, "owner.json"), auto_unbox = TRUE)
    }
  }
  path <- file.path(root, "retry.rds")
  saveRDS(state, path)
  slurmReplies(e, list(c("123_1|COMPLETED|", "123_2|OUT_OF_MEMORY|"), character(), "456", "457"))
  e$.slurmWatch(path, 0L)
  saved <- readRDS(path)
  expect_identical(saved$waves[[2L]]$chunks, 2L)
  expect_identical(saved$waves[[2L]]$memory_mb, 10000L)
  expect_true("--mem=10000M" %in% e$calls[[3L]]$args)
  expect_true("--array=2%2" %in% e$calls[[3L]]$args)
  expect_true("--no-requeue" %in% e$calls[[3L]]$args)
  expect_true("--dependency=afterany:456" %in% e$calls[[4L]]$args)
  expect_true("--array" %in% e$calls[[1L]]$args)
  expect_false(any(dir.exists(file.path(locks, names(owners)[1:2]))))
  expect_true(all(dir.exists(file.path(locks, names(owners)[3:8]))))
  expect_error(e$.slurmWatch(path, 0L), "already been handled")
  slurmReplies(e, list("456_2|OUT_OF_MEMORY|", character(), "789", "790"))
  e$.slurmWatch(path, 1L)
  expect_identical(readRDS(path)$waves[[3L]]$memory_mb, 20000L)
  expect_true("--mem=20000M" %in% e$calls[[3L]]$args)
  slurmReplies(e, list("789_2|OUT_OF_MEMORY|", character()))
  e$.slurmWatch(path, 2L)
  expect_length(e$calls, 2L)
  expect_length(readRDS(path)$waves, 3L)
})

test_that("incomplete accounting and queued arrays cause no retry or cleanup", {
  for (queued in c(FALSE, TRUE)) {
    e <- slurmHelper()
    root <- withr::local_tempdir()
    path <- file.path(root, "retry.rds")
    saveRDS(slurmState(root), path)
    replies <- if (queued) rep(list(c("123_1|COMPLETED|", "123_2|OUT_OF_MEMORY|"), "123"), 6L) else
      rep(list(character()), 6L)
    slurmReplies(e, replies)
    e$.slurmCleanLocks <- function(...) stop("cleanup must not run")
    expect_error(e$.slurmWatch(path, 0L), if (queued) "squeue" else "incomplete")
    expect_identical(readRDS(path)$waves[[1L]]$status, "submitted")
    expect_false(dir.exists(paste0(path, ".watch.lock")))
  }
})

test_that("R scheduler wrapper preserves argument boundaries", {
  skip_on_os("windows")
  e <- slurmHelper()
  value <- "space ' quote $HOME ; literal"
  expect_identical(e$.slurmCommand("printf", c("%s", value)), value)
})


test_that("non-OOM failures are recorded without automatic submission", {
  e <- slurmHelper()
  root <- withr::local_tempdir()
  path <- file.path(root, "retry.rds")
  saveRDS(slurmState(root), path)
  slurmReplies(e, list(c("123_1|COMPLETED|", "123_2|FAILED|1:0|"), character()))
  expect_message(e$.slurmWatch(path, 0L), "Chunks failed without confirmed OOM error: 2")
  wave <- readRDS(path)$waves[[1L]]
  expect_identical(wave$failed_chunks, 2L)
  expect_length(wave$oom_chunks, 0L)
  expect_identical(wave$status, "accounted")
  expect_length(e$calls, 2L)
})

test_that("ambiguous scheduler replies preserve submitted array state", {
  for (reply in list("unexpected reply", c("123", "456"))) {
    e <- slurmHelper()
    root <- withr::local_tempdir()
    state <- slurmState(root)
    state$waves <- list()
    path <- file.path(root, "retry.rds")
    slurmReplies(e, list("123", reply))
    expect_error(e$.slurmSubmit(state, path, 1:2, 0L), "Unexpected monitor sbatch reply")
    saved <- readRDS(path)$waves[[1L]]
    expect_identical(saved$job_id, "123")
    expect_identical(saved$status, "submitted")
    expect_null(saved$monitor_job_id)
    expect_length(e$calls, 2L)
  }
})


test_that("submission holds the state guard until monitor registration is saved", {
  e <- slurmHelper()
  root <- withr::local_tempdir()
  state <- slurmState(root)
  state$waves <- list()
  state$chunks <- 1:2
  path <- file.path(root, "submission.rds")
  saveRDS(state, path)
  calls <- 0L
  e$.slurmCommand <- function(command, args) {
    calls <<- calls + 1L
    expect_true(dir.exists(paste0(path, ".watch.lock")))
    if (calls == 2L) expect_error(e$.slurmWatch(path, 0L), "Monitor guard exists")
    as.character(122L + calls)
  }
  e$.runHPC(c(path, "submit"))
  expect_false(dir.exists(paste0(path, ".watch.lock")))
  expect_identical(readRDS(path)$waves[[1L]]$monitor_job_id, "124")
  expect_error(e$.runHPC(c(path, "submit")), "already attempted")
  expect_identical(calls, 2L)
})
