# Small args for a quick vanilla run.
wfArgs <- function(orb_dir, ...) {
  c(list(orb_dir = orb_dir, model_types = "vanilla", seeds = 31L, min_genomes = 4, n_fold = 2,
         save_fits = FALSE, BPPARAM = BiocParallel::SerialParam()), list(...))
}

test_that("a missing ORB is chosen from the registry", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  orb <- fxOrb(metadata = fxRows(10, 10), parquets = list(gene_count = fxSignalGenes(10, 10)))
  withr::local_dir(withr::local_tempdir())
  local_mocked_bindings(.selectORB = function(call, argument = "path") orb)
  fit <- suppressWarnings(do.call(runModelingWorkflow, wfArgs(NULL, out_dir = "run")))
  expect_equal(fit$built$out_dir, normalizePath("run"))
})

test_that("a missing ORB errors outside an interactive session", {
  local_mocked_bindings(is_interactive = function() FALSE, .package = "rlang")
  expect_error(runModelingWorkflow(out_dir = withr::local_tempdir()), class = "amrml_invalid_argument")
})

test_that("resume without an ORB uses the saved one, or asks if it has moved", {
  local_mocked_bindings(.tuningGrid = fxSmallGrid)
  orb <- fxOrb(metadata = fxRows(10, 10), parquets = list(gene_count = fxSignalGenes(10, 10)))
  withr::local_dir(withr::local_tempdir())
  args <- wfArgs(orb, out_dir = "run", save_state = TRUE)
  suppressWarnings(do.call(runModelingWorkflow, args))
  local_mocked_bindings(.selectORB = function(call, argument = "path") stop("must not ask"),
                        .fitModel = function(...) stop("completed models must not refit"))
  args$orb_dir <- NULL
  fit <- do.call(runModelingWorkflow, c(args, list(resume = TRUE)))
  expect_true(all(fit$tasks$status == "fitted"))

  # The saved ORB has moved: the user is asked, and the run is refused as a different ORB.
  saved <- readRDS(file.path("run", "state", "workflow.rds"))
  saved$settings$orb_dir <- "/no/such/orb"
  saveRDS(saved, file.path("run", "state", "workflow.rds"))
  asked <- FALSE
  local_mocked_bindings(.selectORB = function(call, argument = "path") { asked <<- TRUE; orb })
  expect_error(do.call(runModelingWorkflow, c(args, list(resume = TRUE))), class = "amrml_invalid_argument")
  expect_true(asked)
})
