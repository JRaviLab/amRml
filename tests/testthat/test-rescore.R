# Rescoring without a fit in memory: reopening completed runs from disk.

# A completed run folder, as runModelingWorkflow() leaves it. Its matrices are relabelled as
# vanilla so the run can be rescored; `workflow.rds` holds the built matrices.
fxCompletedRun <- function(env = parent.frame()) {
  fit <- fxRunFit(fxFitTasks(env = env))
  fit$built$matrices$mode_id <- "vanilla"
  saveRDS(list(settings = list(), built = fit$built), file.path(fit$built$out_dir, "workflow.rds"))
  fit
}

test_that("a completed run folder is reopened as the fit that wrote it", {
  fit <- fxCompletedRun()
  reopened <- .readFitRun(fit$built$out_dir, rlang::current_env())

  expect_s3_class(reopened, "amr_fit")
  expect_equal(reopened$tasks, fit$tasks)
  expect_equal(reopened$performance, fit$performance)
  expect_equal(reopened$importance, fit$importance)
  expect_equal(reopened$built$out_dir, fit$built$out_dir)
  expect_true(.hasRescorableFits(reopened))
})

test_that("a run folder that was moved is read from where it is now", {
  fit <- fxCompletedRun()
  moved <- withr::local_tempdir("moved-run")
  file.copy(list.files(fit$built$out_dir, full.names = TRUE), moved, recursive = TRUE)
  reopened <- .readFitRun(moved, rlang::current_env())

  expect_equal(reopened$built$out_dir, normalizePath(moved, winslash = "/"))
})

test_that("an unfinished run folder is refused with what is missing", {
  fit <- fxCompletedRun()
  unlink(file.path(fit$built$out_dir, "results", "fit-completion.rds"))

  expect_error(
    .readFitRun(fit$built$out_dir, rlang::current_env()),
    class = "amrml_invalid_argument"
  )
  expect_error(
    .readFitRun(file.path(fit$built$out_dir, "nowhere"), rlang::current_env()),
    class = "amrml_invalid_argument"
  )
})

test_that("MDR models do not count toward rescoring", {
  fit <- fxCompletedRun()
  fit$built$matrices$unit <- "mdr"

  expect_false(.hasRescorableFits(fit))
})

test_that("the newest rescorable run in an ORB is found, skipping incomplete ones", {
  orb <- withr::local_tempdir("orb")
  runs <- file.path(orb, "amRml", c("2026-01-01_000000", "2026-01-02_000000", "2026-01-03_000000"))
  fit <- fxCompletedRun()
  dir.create(dirname(runs[[1]]))
  for (run in runs[1:2]) {
    dir.create(run)
    file.copy(list.files(fit$built$out_dir, full.names = TRUE), run, recursive = TRUE)
  }
  # The newest folder is unfinished.
  dir.create(runs[[3]])
  file.create(file.path(runs[[3]], "workflow.rds"))

  found <- .findFitRun(orb, rlang::current_env())
  expect_equal(found$built$out_dir, normalizePath(runs[[2]], winslash = "/"))
})

test_that("an ORB without a rescorable run is an error", {
  orb <- withr::local_tempdir("orb")
  expect_error(.findFitRun(orb, rlang::current_env()), class = "amrml_invalid_argument")

  # A completed run with no vanilla models is not rescorable.
  fit <- fxRunFit(fxFitTasks())
  saveRDS(list(settings = list(), built = fit$built), file.path(fit$built$out_dir, "workflow.rds"))
  dir.create(file.path(orb, "amRml"))
  file.copy(fit$built$out_dir, file.path(orb, "amRml"), recursive = TRUE)
  expect_error(.findFitRun(orb, rlang::current_env()), class = "amrml_invalid_argument")
})

test_that("rescoreFeatures() takes at most one source and needs a fit that can be rescored", {
  fit <- fxCompletedRun()

  expect_error(rescoreFeatures(fit, run_dir = fit$built$out_dir), class = "amrml_invalid_argument")
  expect_error(rescoreFeatures(fit, orb_dir = "x"), class = "amrml_invalid_argument")
  expect_error(rescoreFeatures("not a fit"), class = "amrml_invalid_argument")
  expect_error(rescoreFeatures(run_dir = file.path(fit$built$out_dir, "nowhere")),
               class = "amrml_invalid_argument")
  # Without a path it needs the registry, which a non-interactive session cannot choose from.
  expect_error(rescoreFeatures(), class = "amrml_invalid_argument")
})

# A small result in the shape runFeatureDyadDiscovery() returns.
fxRescored <- function() {
  fit <- function(drug, mcc, shuffled) {
    tibble::tibble(
      species = "Sfl", drug_label = "drug", drug_or_class = drug, seed = 1L,
      feature_type = "AMRFinder", feature_subtype = "counts", model = "LR",
      fit_penalty = 0.1, fit_mixture = 0.5,
      nonshuffled_MCC = mcc, shuffled_MCC = shuffled, MCC_diff = mcc - shuffled
    )
  }
  good <- dplyr::bind_rows(fit("AMP", 0.9, 0.1), fit("CIP", 0.7, 0.2))
  separated <- dplyr::bind_rows(fit("AMP", 0.9, 0.1), fit("CIP", 0.7, 0.2), fit("TMP", 0.3, 0.1))
  feature <- function(drug, variable) {
    tibble::tibble(
      species = "Sfl", drug_label = "drug", drug_or_class = drug, feature_type = "AMRFinder",
      feature_subtype = "counts", variable = variable, median_rank_score = 0.9
    )
  }
  top <- dplyr::bind_rows(feature("AMP", "sul2"), feature("CIP", "sul2"), feature("CIP", "gyrA"))
  sharedness <- tibble::tibble(
    species = "Sfl", drug_label = "drug", feature_type = "AMRFinder",
    variable = c("sul2", "gyrA"), n_drug_or_class = c(2L, 1L), drug_or_class_csv = c("AMP, CIP", "CIP")
  )
  dyad <- function(drug, id) {
    tibble::tibble(species = "Sfl", drug_label = "drug", drug_or_class = drug, dyad = id)
  }
  dyads <- dplyr::bind_rows(dyad("AMP", "d1"), dyad("CIP", "d1"), dyad("CIP", "d2"))
  structure(
    list(
      models = list(
        good_performance = good, good_shuffled_separation = separated, qualified = good,
        summary = tibble::tibble(species = "Sfl", drug_or_class = c("AMP", "CIP"), n_passing_fits = 1L)
      ),
      features = list(
        top = top, consistent = top[1:2, ], high_contribution = top[3, ], all_seeds = top,
        sharedness = sharedness
      ),
      dyads = list(
        top = dyads, high_score = dyads[1, ],
        sharedness = tibble::tibble(
          species = "Sfl", drug_label = "drug", dyad = c("d1", "d2"),
          n_drug_or_class = c(2L, 1L), drug_or_class_csv = c("AMP, CIP", "CIP")
        )
      ),
      parameters = list(MCC_threshold = 0.5, min_models_shared = 2)
    ),
    class = c("feature_dyad_discovery", "list")
  )
}

test_that("a result is written as four tables, with its subsets kept as flags", {
  run <- withr::local_tempdir("run")
  out <- .writeRescored(fxRescored(), run)

  expect_equal(dirname(out), normalizePath(file.path(run, "magic"), winslash = "/"))
  expect_setequal(
    list.files(out),
    c("optimal_models.parquet", "model_summary.parquet", "qualified_features.parquet",
      "qualified_dyads.parquet", "parameters.json")
  )

  perf <- arrow::read_parquet(file.path(out, "optimal_models.parquet"))
  expect_equal(nrow(perf), 3L)
  expect_equal(perf$is_qualified[order(perf$drug_or_class)], c(TRUE, TRUE, FALSE))
  expect_equal(perf$is_separated_from_shuffled, rep(TRUE, 3))

  features <- arrow::read_parquet(file.path(out, "qualified_features.parquet"))
  expect_equal(nrow(features), 3L)
  expect_equal(features$is_consistent, c(TRUE, TRUE, FALSE))
  expect_equal(features$is_high_contribution, c(FALSE, FALSE, TRUE))
  expect_equal(features$is_shared, c(TRUE, TRUE, FALSE))
  expect_equal(features$is_unique, c(FALSE, FALSE, TRUE))

  dyads <- arrow::read_parquet(file.path(out, "qualified_dyads.parquet"))
  expect_equal(dyads$is_high_score, c(TRUE, FALSE, FALSE))
  expect_equal(dyads$is_shared, c(TRUE, TRUE, FALSE))

  info <- jsonlite::read_json(file.path(out, "parameters.json"))
  expect_equal(info$parameters$MCC_threshold, 0.5)
  expect_equal(info$rows$qualified_features, 3L)
})

test_that("every call writes a new folder and keeps the earlier result", {
  run <- withr::local_tempdir("run")
  first <- .writeRescored(fxRescored(), run)
  second <- .writeRescored(fxRescored(), run)

  expect_false(identical(first, second))
  expect_true(file.exists(file.path(first, "qualified_features.parquet")))
  expect_length(list.files(file.path(run, "magic")), 2L)
})
