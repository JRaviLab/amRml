# Tests for readORB(), in the order it reads an ORB.

# The path argument

test_that("readORB() requires a single directory path", {
  expect_error(readORB(), class = "amrml_invalid_argument")
  expect_error(readORB(c("a", "b")), class = "amrml_invalid_argument")
  expect_error(readORB(""), class = "amrml_invalid_argument")
})

test_that("readORB() tells a missing directory from one with no manifest", {
  expect_error(
    readORB(file.path(tempdir(), "definitely-not-here")),
    class = "amrml_orb_dir_not_found"
  )
  expect_error(readORB(withr::local_tempdir()), class = "amrml_manifest_not_found")
})

# Finding manifests

test_that("only top-level manifest files are found", {
  dir <- withr::local_tempdir()
  writeLines("42", file.path(dir, "settings.json"))
  nested <- file.path(dir, "nested")
  dir.create(nested)
  fxWriteManifest(nested, fxManifest())

  expect_error(readORB(dir), class = "amrml_manifest_not_found")
})

# The checked record

test_that("declared files keep each path with its size, and modified_at is optional", {
  run <- fxManifest()$runs[[1]]
  run$stages[[1]]$outputs[[2]]$modified_at <- NULL

  expect_equal(.declaredFiles(run), tibble::tibble(
    path = c("/orb/metadata.parquet", "/orb/gene_count.parquet", "/orb/Tst_parquet.duckdb"),
    name = c("metadata.parquet", "gene_count.parquet", "Tst_parquet.duckdb"),
    size_bytes = c(100, 200, 0),
    modified_at = c("2026-01-01 00:00:01", NA_character_, "2026-01-01 00:00:01")
  ))
})

test_that("declared file sizes must be non-negative numbers", {
  for (size in list("100", Inf, -1)) {
    run <- fxManifest()$runs[[1]]
    run$stages[[1]]$outputs[[1]]$size_bytes <- size
    expect_error(.declaredFiles(run), "invalid path or size", fixed = TRUE)
  }
})

test_that("the record keeps the fields the reader uses", {
  expect_silent(record <- .manifestRecord("manifest.json", fxManifest()))

  expect_equal(record$manifest_path, "manifest.json")
  expect_equal(record$dataset_id, "Tst")
  expect_equal(record$directory, "/orb")
  expect_equal(record$metadata_parquet, "/orb/metadata.parquet")
  expect_equal(record$producer, "amRdata")
  expect_equal(record$producer_run_id, "run_fixture_producer")
  expect_equal(record$finished_at, "2026-01-01 00:00:01")
})

test_that("the dataset label lists what amRdata was asked for, or the dataset ID", {
  manifest <- fxManifest()
  manifest$dataset$selection$user_bacs <- list("Species one", "Species two")
  expect_equal(.manifestRecord("m.json", manifest)$dataset_label, "Species one, Species two")

  manifest$dataset$selection$user_bacs <- NULL
  expect_equal(.manifestRecord("m.json", manifest)$dataset_label, "Tst")
})

test_that("the record follows the artifact's producer run and its export stage", {
  manifest <- fxManifest()
  producer <- manifest$runs[[1]]
  other_run <- producer
  other_run$run_id <- "unrelated"
  other_run$finished_at <- "2030-01-01 00:00:00"
  last_run <- other_run
  last_run$run_id <- "also_unrelated"
  producer$stages <- c(
    list(list(name = "before_export", outputs = list())), producer$stages,
    list(list(name = "after_export", outputs = list()))
  )
  # The producer run and its export stage are neither first nor last.
  manifest$runs <- list(other_run, producer, last_run)

  record <- .manifestRecord("m.json", manifest)
  expect_equal(record$producer_run_id, "run_fixture_producer")
  expect_equal(record$finished_at, "2026-01-01 00:00:01")
  expect_equal(record$files$name, c("metadata.parquet", "gene_count.parquet", "Tst_parquet.duckdb"))
})

test_that("an unready artifact or failed producer run is returned as an error, not thrown", {
  expect_silent(error <- .manifestRecord("m.json", fxManifest(artifact_status = "building")))
  expect_s3_class(error, "amrml_orb_not_ready")
  expect_equal(error$observed$status, "building")

  manifest <- fxManifest()
  manifest$artifacts <- list()
  expect_equal(.manifestRecord("m.json", manifest)$observed$status, "absent")

  expect_silent(error <- .manifestRecord("m.json", fxManifest(producer_status = "failed")))
  expect_s3_class(error, "amrml_producer_run_invalid")
  expect_equal(error$observed, list(producer_run_id = "run_fixture_producer", status = "failed"))

  manifest <- fxManifest()
  manifest$artifacts$amRml_input$producer_run_id <- "missing_run"
  expect_equal(.manifestRecord("m.json", manifest)$observed$status, "missing")
})

test_that("the metadata file must be among the declared files", {
  manifest <- fxManifest()
  manifest$artifacts$amRml_input$metadata_parquet <- "/orb/undeclared.parquet"

  expect_error(.manifestRecord("m.json", manifest), "undeclared.parquet", fixed = TRUE)
})

# Checking each manifest

test_that("unreadable JSON and invalid contents are told apart", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, fxManifestName())

  writeLines("{ not json", path)
  expect_silent(error <- .manifestReadiness(path))
  expect_s3_class(error, "amrml_manifest_unreadable")

  writeLines("42", path)
  expect_silent(error <- .manifestReadiness(path))
  expect_s3_class(error, "amrml_manifest_invalid")
})

test_that("a malformed field is recorded as manifest_invalid, naming the field", {
  broken <- list(
    "amRdata v1" = function(m) {
      m$schema_version <- 2L
      m
    },
    "metadata_parquet" = function(m) {
      m$artifacts$amRml_input$metadata_parquet <- NULL
      m
    },
    "finished_at" = function(m) {
      m$runs[[1]]$finished_at <- "invalid"
      m
    },
    "clean_metadata_and_export" = function(m) {
      m$runs[[1]]$stages <- list()
      m
    },
    # `$` would have matched path_original as path; the field must be named exactly.
    "invalid path or size" = function(m) {
      output <- m$runs[[1]]$stages[[1]]$outputs[[1]]
      names(output)[names(output) == "path"] <- "path_original"
      m$runs[[1]]$stages[[1]]$outputs[[1]] <- output
      m
    }
  )
  dir <- withr::local_tempdir()
  path <- file.path(dir, fxManifestName())

  for (reason in names(broken)) {
    fxWriteManifest(dir, broken[[reason]](fxManifest()))
    expect_silent(error <- .manifestReadiness(path))
    expect_s3_class(error, "amrml_manifest_invalid")
    expect_match(error$observed$error, reason, fixed = TRUE)
  }

  # A shape the checks don't expect is still manifest_invalid, not an R error.
  manifest <- fxManifest()
  manifest$artifacts$amRml_input <- 42
  fxWriteManifest(dir, manifest)
  expect_s3_class(.manifestReadiness(path), "amrml_manifest_invalid")
})

# Choosing the newest manifest

test_that(".newestReady() picks the latest finish, then the first file name", {
  older <- list(manifest_path = "a.json", finished_at = "2026-01-01 00:00:01.5")
  newer <- list(manifest_path = "z.json", finished_at = "2026-01-01 00:00:02")
  expect_identical(.newestReady(list(older, newer)), newer)
  expect_identical(.newestReady(list(newer, older)), newer)

  tied <- list(manifest_path = "b.json", finished_at = "2026-01-01 00:00:02")
  expect_identical(.newestReady(list(newer, tied)), tied)
  expect_identical(.newestReady(list(tied, newer)), tied)

  # Fractional seconds sort correctly as text.
  quarter <- list(manifest_path = "q.json", finished_at = "2026-01-01 00:00:01.25")
  expect_identical(.newestReady(list(older, quarter)), older)
})

test_that("the newest producer run wins, not the newest file, with a message", {
  dir <- fxOrb()

  # An older run in a newer file.
  older <- fxReadManifest(dir)
  older$runs[[1]]$finished_at <- "2020-01-01 00:00:00"
  fxWriteManifest(dir, older, "manifest_run_older.json")
  Sys.setFileTime(file.path(dir, "manifest_run_older.json"), Sys.time() + 3600)

  expect_message(
    orb <- readORB(dir),
    "manifest_run_fixture.json, the newest of 2 ready manifests"
  )
  expect_equal(basename(orb$manifest_path), fxManifestName())

  newer <- fxReadManifest(dir)
  newer$runs[[1]]$finished_at <- "2030-01-01 00:00:00"
  fxWriteManifest(dir, newer, "manifest_run_newer.json")

  expect_message(orb <- readORB(dir), "newest of 3 ready manifests")
  expect_equal(basename(orb$manifest_path), "manifest_run_newer.json")
})

test_that("a building manifest beside a ready one does not stop the ready one", {
  dir <- fxOrb()
  building <- fxReadManifest(dir)
  building$artifacts$amRml_input$status <- "building"
  fxWriteManifest(dir, building, "manifest_run_failed.json")

  expect_equal(basename(readORB(dir)$manifest_path), fxManifestName())
})

test_that("when no manifest is usable, every manifest's reason is reported", {
  dir <- fxOrb(artifact_status = "building")
  broken <- fxReadManifest(dir)
  broken$artifacts$amRml_input$status <- "ready"
  broken$runs[[1]]$finished_at <- "invalid"
  fxWriteManifest(dir, broken, "manifest_run_broken.json")

  cnd <- rlang::catch_cnd(readORB(dir), classes = "amrml_error")
  expect_s3_class(cnd, "amrml_orb_no_usable_manifest")
  expect_s3_class(cnd$observed[[fxManifestName()]], "amrml_orb_not_ready")
  expect_s3_class(cnd$observed$manifest_run_broken.json, "amrml_manifest_invalid")

  message <- conditionMessage(cnd)
  for (text in c(fxManifestName(), "building", "manifest_run_broken.json", "finished_at")) {
    expect_match(message, text, fixed = TRUE)
  }
})

# Checking the ORB's files

test_that("a moved ORB is refused as moved, not as missing files", {
  dir <- fxOrb()
  destination <- file.path(withr::local_tempdir(), "moved")
  expect_true(file.rename(dir, destination))

  cnd <- rlang::catch_cnd(readORB(destination), classes = "amrml_error")
  expect_s3_class(cnd, "amrml_orb_moved")
  expect_equal(cnd$observed$actual, normalizePath(destination))
  expect_match(conditionMessage(cnd), basename(dir), fixed = TRUE)
})

test_that("a missing declared file is refused", {
  dir <- fxOrb()
  unlink(file.path(dir, "weird.parquet"))

  cnd <- rlang::catch_cnd(readORB(dir), classes = "amrml_error")
  expect_s3_class(cnd, "amrml_orb_files_missing")
  expect_equal(basename(cnd$observed$missing), "weird.parquet")
  expect_equal(rlang::call_name(cnd$call), "readORB")
})

test_that("an unreadable declared file counts as missing", {
  skip_on_os("windows")
  skip_if(Sys.info()[["user"]] == "root", "root can read any file")

  dir <- fxOrb()
  path <- file.path(dir, "weird.parquet")
  Sys.chmod(path, "000")
  withr::defer(Sys.chmod(path, "644"))

  cnd <- rlang::catch_cnd(readORB(dir), classes = "amrml_error")
  expect_s3_class(cnd, "amrml_orb_files_missing")
  expect_equal(basename(cnd$observed$missing), "weird.parquet")
})

test_that("the file check reports only files whose size changed, with both sizes", {
  dir <- normalizePath(withr::local_tempdir())
  paths <- file.path(dir, c("first", "second", "empty"))
  writeBin(charToRaw("abc"), paths[[1]])
  writeBin(charToRaw("12345"), paths[[2]])
  file.create(paths[[3]])
  record <- list(directory = dir, files = tibble::tibble(
    path = paths, name = c("first", "second", "empty"), size_bytes = c(3, 5, 0)
  ))
  expect_no_error(.checkOrbFiles(record, dir, call = NULL))

  record$files$size_bytes <- c(1, 10, 0)
  error <- rlang::catch_cnd(.checkOrbFiles(record, dir, call = NULL), classes = "error")
  expect_s3_class(error, "amrml_orb_file_changed")
  expect_equal(error$observed, list(
    file = paths[1:2], size_bytes = c(1, 10), actual_bytes = c(3, 5)
  ))
})

test_that("a failed file check never falls back to an older usable manifest", {
  dir <- fxOrb()

  # The older manifest doesn't declare weird.parquet, so on its own it reads.
  newest <- fxReadManifest(dir)
  older <- newest
  older$runs[[1]]$finished_at <- "2020-01-01 00:00:00"
  older$runs[[1]]$stages[[1]]$outputs <- Filter(
    function(o) basename(o$path) != "weird.parquet", older$runs[[1]]$stages[[1]]$outputs
  )
  fxWriteManifest(dir, older)
  unlink(file.path(dir, "weird.parquet"))
  expect_equal(basename(readORB(dir)$manifest_path), fxManifestName())

  fxWriteManifest(dir, newest, "manifest_run_newer.json")
  expect_error(
    expect_message(readORB(dir), "newest of 2 ready manifests"),
    class = "amrml_orb_files_missing"
  )
})

# Feature tables

test_that("feature tables are found by their columns, whatever their order or file name", {
  dir <- withr::local_tempdir()
  tables <- list(
    a_widget = tibble::tibble(value = 1, widget = "w", genome_id = "g"),
    z_gene = tibble::tibble(gene = "a", genome_id = "g", value = 1),
    no_genome = tibble::tibble(other = "g", gene = "a", value = 1),
    no_value = tibble::tibble(genome_id = "g", gene = "a", other = 1),
    extra = tibble::tibble(genome_id = "g", gene = "a", value = 1, extra = 2)
  )
  paths <- file.path(dir, paste0(names(tables), ".parquet"))
  for (i in seq_along(paths)) {
    arrow::write_parquet(tables[[i]], paths[[i]])
  }
  files <- tibble::tibble(path = paths, name = basename(paths))

  expect_equal(.featureTables(files, call = NULL), tibble::tibble(
    feature_type = c("gene", "widget"), path = paths[c(2, 1)]
  ))
})

test_that("feature tables sort the same way in any locale", {
  # Tests run in the C locale; in en_US a plain sort puts "gene" before "Pfam".
  skip_on_os("windows")
  skip_if_not("en_US.UTF-8" %in% system2("locale", "-a", stdout = TRUE))
  withr::local_collate("en_US.UTF-8")
  skip_if_not(identical(sort(c("Pfam", "gene")), c("gene", "Pfam")))

  dir <- withr::local_tempdir()
  paths <- file.path(dir, c("gene.parquet", "pfam.parquet"))
  arrow::write_parquet(tibble::tibble(genome_id = "g", gene = "a", value = 1), paths[[1]])
  arrow::write_parquet(tibble::tibble(genome_id = "g", Pfam = "a", value = 1), paths[[2]])
  files <- tibble::tibble(path = paths, name = basename(paths))

  expect_equal(.featureTables(files, call = NULL)$feature_type, c("Pfam", "gene"))
})

test_that("undeclared parquet files are ignored", {
  dir <- fxOrb()
  writeLines("not Parquet", file.path(dir, "undeclared.parquet"))

  orb <- readORB(dir)
  expect_equal(orb$feature_tables$feature_type, c("gene", "widget"))
  expect_false("undeclared.parquet" %in% orb$files$name)
})

test_that("an unreadable declared parquet is an error, not a smaller feature set", {
  dir <- fxOrb()
  writeLines("this is not parquet", file.path(dir, "weird.parquet"))
  fxRecordSizes(dir)

  cnd <- rlang::catch_cnd(readORB(dir), classes = "amrml_error")
  expect_s3_class(cnd, "amrml_parquet_unreadable")
  expect_equal(basename(cnd$observed$file), "weird.parquet")
})

test_that("an ORB with no feature tables is refused", {
  dir <- fxOrb(parquets = fxDefaultParquets()["decoy_count"])

  expect_error(readORB(dir), class = "amrml_orb_no_feature_tables")
})

test_that("two feature tables sharing a feature ID column are refused", {
  parquets <- fxDefaultParquets()
  parquets$gene_again <- parquets$gene_count
  dir <- fxOrb(parquets = parquets)

  cnd <- rlang::catch_cnd(readORB(dir), classes = "amrml_error")
  expect_s3_class(cnd, "amrml_orb_duplicate_feature_type")
  expect_equal(cnd$observed$feature_type, "gene")
})

# The ORB object

test_that("readORB() opens a valid ORB and describes it", {
  dir <- fxOrb()
  expect_no_message(orb <- readORB(file.path(dir, ".")))

  expect_s3_class(orb, "amr_orb")
  expect_setequal(names(orb), c(
    "dataset_id", "dataset_label", "directory", "manifest_path", "metadata_parquet",
    "producer", "producer_run_id", "finished_at", "files", "feature_tables"
  ))
  expect_equal(orb$dataset_id, "Tst")
  expect_equal(orb$dataset_label, "Testus fixtus")
  expect_equal(orb$directory, normalizePath(dir))
  expect_equal(orb$manifest_path, file.path(normalizePath(dir), fxManifestName()))
  expect_equal(orb$metadata_parquet, file.path(dir, "metadata.parquet"))
  expect_equal(orb$producer, "amRdata")
  expect_equal(orb$producer_run_id, "run_fixture_producer")
  expect_equal(orb$finished_at, "2026-01-01 00:00:01")
  expect_setequal(orb$files$name, c(
    "metadata.parquet", "gene_count.parquet", "weird.parquet",
    "decoy_count.parquet", "Tst_parquet.duckdb"
  ))
  expect_equal(orb$feature_tables, tibble::tibble(
    feature_type = c("gene", "widget"),
    path = file.path(dir, c("gene_count.parquet", "weird.parquet"))
  ))
})

test_that("readORB() does not write to the ORB", {
  dir <- fxOrb()
  before <- fxOrbState(dir)

  readORB(dir)

  expect_equal(fxOrbState(dir), before)
})

test_that("printing an ORB shows its summary and returns it invisibly", {
  orb <- readORB(fxOrb())
  output <- paste(capture.output(result <- withVisible(print(orb))), collapse = "\n")

  for (text in c(
    "Tst", "Testus fixtus", orb$directory, fxManifestName(),
    "2026-01-01 00:00:01", "5 files", "gene, widget"
  )) {
    expect_match(output, text, fixed = TRUE)
  }
  expect_identical(result$value, orb)
  expect_false(result$visible)
})
