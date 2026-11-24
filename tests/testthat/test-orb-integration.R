# Tests for reading and profiling a real amRdata ORB, set by AMRML_TEST_ORB.

test_that("a real ORB reads and profiles consistently without changing its files", {
  path <- Sys.getenv("AMRML_TEST_ORB")
  skip_if_not(nzchar(path), "AMRML_TEST_ORB is not set")
  expect_true(dir.exists(path), info = "AMRML_TEST_ORB must point to an existing directory")

  # Top-level files only; genomes/ and panaroo_out_* are not read.
  entries <- list.files(path, all.files = TRUE, no.. = TRUE)
  files <- file.path(path, entries)
  files <- files[!dir.exists(files)]
  before <- tools::md5sum(files)

  orb <- readORB(path)
  profile <- profileORB(orb)

  expect_gt(nrow(orb$files), 0)
  expect_gt(nrow(orb$feature_tables), 0)

  targets <- profile$targets
  metadata_genomes <- unique(arrow::read_parquet(orb$metadata_parquet)$genome.genome_id)
  expect_equal(targets$n_resistant + targets$n_susceptible, targets$n_genomes)
  expect_true(all(profile$phenotypes$genome_id %in% profile$genomes$genome_id))
  expect_setequal(
    metadata_genomes,
    c(profile$genomes$genome_id, profile$excluded_genomes$genome_id)
  )

  expect_equal(list.files(path, all.files = TRUE, no.. = TRUE), entries)
  expect_equal(tools::md5sum(files), before)
})
