test_that("sparse storage keeps a labelled genome with no observed features", {
  metadata <- fxRows(10, 10)
  genes <- fxSignalGenes(10, 10)
  genes <- genes[genes$genome_id != "g020", ]
  eligibility <- fxEligibility(metadata, modes = "vanilla", min_genomes = 4,
                               n_fold = 2, parquets = list(gene_count = genes))
  out <- withr::local_tempdir()
  built <- buildMatrices(planMatrices(eligibility), out)
  store <- readRDS(file.path(built$out_dir, unique(built$matrices$file)[[1]]))
  row <- match("g020", store$genomes)
  expect_false(is.na(row))
  expect_equal(as.numeric(store$values[row, ]), rep(0, ncol(store$values)))
  index <- readRDS(file.path(built$out_dir, "matrix_index.rds"))
  expect_true(all(vapply(index$rows, function(rows) "g020" %in% rows$genome_id, logical(1))))
})

test_that("feature presets select available tables and explicit names stay strict", {
  eligibility <- fxEligibility(fxRows(4, 4), modes = "vanilla", min_genomes = 4, n_fold = 2)
  for (preset in c("quick", "standard", "full")) {
    plan <- planMatrices(eligibility, feature_types = preset)
    expect_equal(unique(plan$matrices$feature_type), "gene")
  }
  explicit <- planMatrices(eligibility, feature_types = "widget")
  expect_equal(unique(explicit$matrices$feature_type), "widget")
  expect_error(planMatrices(eligibility, feature_types = "typo"), class = "amrml_invalid_argument")
})
