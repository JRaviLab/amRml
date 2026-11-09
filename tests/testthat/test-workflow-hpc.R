test_that("HPC stages preserve workflow settings, fit chunks and finish without refitting", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid, .registerManifest = function(path) path)
  orb <- fxOrb(metadata = fxRows(20, 20), parquets = list(gene_count = fxSignalGenes(20, 20)))
  out <- file.path(withr::local_tempdir(), "run with 'quotes'")
  commands <- list()
  local_mocked_bindings(.slurmCommand = function(command, args) {
    commands[[length(commands) + 1L]] <<- c(command, args)
    switch(command, sbatch = as.character(120L + length(commands)),
      sacct = c("121_1|COMPLETED|0:0|", "121_2|COMPLETED|0:0|"), squeue = character())
  })
  prep <- suppressWarnings(runModelingWorkflowHPC("prep", out_dir = out,
    orb_dir = orb, model_types = "vanilla", feature_types = "gene", seeds = 31L,
    min_genomes = 4, n_fold = 2, n_chunks = 2, save_state = TRUE))
  expect_length(commands, 0L)
  expect_identical(prep$n_chunks, 2L)
  expect_equal(sum(prep$models$n_models), prep$n_models)
  expect_true(is.data.frame(prep$eligibility))
  job <- runModelingWorkflowHPC("run", out_dir = out, memory_mb = c(3700L, 10000L, 20000L),
    scheduler_options = list(account = "account", qos = "qos"))
  expect_length(commands, 0L) # Generating instructions never calls a scheduler.
  expect_equal(job$n_models, prep$n_models)
  expect_setequal(list.files(file.path(out, "hpc")), c("run.sh", "submission.rds"))
  expect_match(job$commands, "submit$")
  expect_error(runModelingWorkflowHPC("run", out_dir = out, memory_mb = 3700L), "already prepared")
  .runHPC(c(job$state_file, "submit"))
  expect_length(commands, 2L)
  expect_true(all(vapply(commands, function(x) "--qos=qos" %in% x, logical(1))))
  expect_true("--mem=3700M" %in% commands[[1L]])
  expect_true("--dependency=afterany:121" %in% commands[[2L]])
  expect_true("--time=00:10:00" %in% commands[[2L]])
  expect_error(.runHPC(c(job$state_file, "submit")), "already attempted")
  expect_error(runModelingWorkflowHPC("finish", out_dir = out), "unfinished")
  for (chunk in 1:2) {
    withr::with_envvar(c(SLURM_ARRAY_TASK_ID = as.character(chunk), SLURM_CPUS_PER_TASK = "1"),
                     suppressWarnings(.runHPC(c(job$state_file, "worker"))))
  }
  .runHPC(c(job$state_file, "watch", "0"))
  local_mocked_bindings(.fitModel = function(...) stop("collection must not refit"))
  collected <- withVisible(runModelingWorkflowHPC("finish", out_dir = out))
  expect_false(collected$visible)
  expect_true(all(collected$value$tasks$status == "fitted"))
  expect_true(file.exists(file.path(out, "results", "merged", "all_perf.parquet")))
})

test_that("resource validation cannot submit jobs or replace saved configuration", {
  tasks <- fxFitTasks()
  fitModels(tasks, execution = "prepare")
  out <- tasks$built$out_dir
  local_mocked_bindings(.slurmCommand = function(...) stop("must not submit"))
  for (memory in list(NULL, numeric(), c(3700, 2000), c(NA, 10000), c(3700, Inf), 0, 3.7)) {
    expect_error(runModelingWorkflowHPC("run", out_dir = out, memory_mb = memory), "memory_mb")
  }
  for (options in list(list("account"), list(mem = "4G"), list(qos = NA_character_),
                       list(qos = c("a", "b")), list(qos = "a", qos = "b"))) {
    expect_error(runModelingWorkflowHPC("run", out_dir = out, memory_mb = 4000,
                                      scheduler_options = options), "scheduler_options")
  }
  expect_error(runModelingWorkflowHPC("run", out_dir = out, memory_mb = 4000,
                                    scheduler = "pbs"), "scheduler")
  expect_error(runModelingWorkflowHPC("run", out_dir = out, memory_mb = 4000,
                                    cpus_per_task = 0), "cpus_per_task")
  expect_error(runModelingWorkflowHPC("run", out_dir = out, memory_mb = 4000,
                                    max_concurrent = 0), "max_concurrent")
  expect_error(runModelingWorkflowHPC("run", out_dir = out, memory_mb = 4000,
                                    n_chunks = 2), "only for stage")
  expect_error(runModelingWorkflowHPC("prep", out_dir = out, execution = "fit"), "execution")
  expect_false(dir.exists(file.path(out, "hpc")))
})

test_that("the single generated runner preserves argument boundaries and preamble", {
  skip_on_os("windows")
  tasks <- fxFitTasks()
  fitModels(tasks, execution = "prepare")
  bin <- withr::local_tempdir()
  fake_r <- file.path(bin, "Rscript")
  writeLines(c("#!/bin/sh", 'printf "%s\\n" "$@"'), fake_r)
  Sys.chmod(fake_r, "0755")
  path <- paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep)
  job <- runModelingWorkflowHPC("run", out_dir = tasks$built$out_dir, memory_mb = 4000L,
                               preamble = paste("export PATH=", shQuote(path), sep = ""))
  expect_identical(system2("bash", c("-n", shQuote(job$runner))), 0L)
  actual <- system2("bash", c(shQuote(job$runner), "worker"), stdout = TRUE)
  expect_identical(actual, c("-e", "amRml:::.runHPC(commandArgs(TRUE))", job$state_file, "worker"))
  expect_length(readRDS(job$state_file)$waves, 0L)
})

test_that("empty HPC plans finish without submission", {
  tasks <- fxFitTasks()
  tasks$tasks <- tasks$tasks[0, ]
  fitModels(tasks, execution = "prepare")
  out <- tasks$built$out_dir
  expect_error(runModelingWorkflowHPC("run", out_dir = out, memory_mb = 4000L), "no chunks")
  fit <- runModelingWorkflowHPC("finish", out_dir = out)
  expect_equal(nrow(fit$tasks), 0L)
})
