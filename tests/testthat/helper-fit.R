# Fixtures for fitting and checkpoint integration tests.
fxFitTasks <- function(sparse = FALSE, labels = c("real", "shuffled"),
                       env = parent.frame()) {
  out <- withr::local_tempdir("fit-run", .local_envir = env)
  dir.create(file.path(out, "features"))
  genomes <- sprintf("g%03d", seq_len(40))
  phenotype <- rep(c("Resistant", "Susceptible"), 20)
  values <- withr::with_seed(817, {
    p <- if (sparse) 24L else 5L
    x <- if (sparse) matrix(stats::rbinom(40 * p, 1, 0.025), 40) else
      matrix(stats::rpois(40 * p, 3), 40)
    x[, 1] <- as.numeric(phenotype == "Resistant") * 3
    x[, 2] <- rep(c(0, 1, 2, 0), 10)
    x
  })
  store <- list(values = as(Matrix::Matrix(values, sparse = TRUE), "dgCMatrix"),
                genomes = genomes, features = paste0("feature", seq_len(ncol(values))))
  saveRDS(store, file.path(out, "features", "0001.rds"))
  rows <- lapply(list(33:36, 37:40), function(test) {
    i <- c(1:32, test)
    tibble::tibble(genome_id = genomes[i], role = c(rep("train", 32), rep("test", 4)),
                   phenotype = phenotype[i])
  })
  matrices <- tibble::tibble(
    matrix_id = c("m1", "m2"), mode_id = "cross_drug", unit = "DRA", target = "drug",
    train_group = "DRA", test_group = c("DRB", "DRC"), feature_type = "gene",
    encoding = "counts", n_train = 32L, n_test = 4L, status = "built",
    rule_id = NA_character_, n_features = ncol(values), file = "features/0001.rds"
  )
  saveRDS(list(rows = rows, matrix_scopes = 1:2,
               columns = rep(list(seq_len(ncol(values))), 2)),
          file.path(out, "matrix_index.rds"))
  arrow::write_parquet(matrices, file.path(out, "matrices.parquet"))
  orb <- list(directory = "/synthetic-orb", dataset_id = "Tst",
              dataset_label = "Synthetic", manifest_path = "/synthetic-orb/manifest.json",
              producer_run_id = "fixture", finished_at = "2026-01-01 00:00:00")
  eligibility <- list(profile = list(orb = orb), settings = list(n_fold = 2L, holdout = 0.2),
                      modes = .MODES[.MODES$mode_id == "cross_drug", ])
  built <- structure(list(eligibility = eligibility, out_dir = normalizePath(out),
                          matrices = matrices, storage = "shared_sparse"), class = "amr_matrices")
  expandTasks(built, seeds = 31L, labels = labels)
}

fxSmallGrid <- function() tibble::tibble(penalty = c(0.01, 0.1), mixture = c(1, 1))

fxComparableFit <- function(fit) {
  fit$performance$run_time_sec <- NULL
  fit[c("tasks", "performance", "importance", "predictions")]
}

fxFitOptions <- function() {
  if ("save_fits" %in% names(formals(fitModels))) list(save_fits = TRUE) else list()
}

fxRunFit <- function(tasks, ...) {
  suppressWarnings(do.call(fitModels, c(list(tasks = tasks), list(...), fxFitOptions())))
}
