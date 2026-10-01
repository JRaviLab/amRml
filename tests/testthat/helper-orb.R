# Test fixtures: synthetic amRdata manifests and ORBs

fxGenomes <- function() c("g1", "g2", "g3", "g4", "g5")

# Column names follow amRdata's metadata.parquet.
fxMetadata <- function() {
  tibble::tibble(
    `genome.genome_id` = rep(fxGenomes(), each = 2),
    `genome_drug.antibiotic` = rep(c("antibiotic_a", "antibiotic_b"), 5),
    `genome_drug.resistant_phenotype` = rep(c("Resistant", "Susceptible"), 5),
    drug_class = "class_x",
    drug_abbr = rep(c("DRA", "DRB"), 5),
    class_abbr = "CLX",
    country_abbr = "CX",
    year_bin = "2000-2004",
    resistant_classes = "CLX",
    num_resistant_classes = 1L
  )
}

# Metadata with one row per genome and drug; the other columns are recycled.
fxMetadataRows <- function(genome, drug, phenotype, class = "CLX",
                           country = "CX", year_bin = "2000-2004") {
  tibble::tibble(
    `genome.genome_id` = genome,
    `genome_drug.resistant_phenotype` = phenotype,
    drug_abbr = drug,
    class_abbr = class,
    country_abbr = country,
    year_bin = year_bin
  )
}

# weird is a feature table with an odd name; decoy_count is named like one but isn't.
fxDefaultParquets <- function() {
  genomes <- fxGenomes()

  list(
    gene_count = tibble::tibble(
      genome_id = rep(genomes, each = 2),
      gene = rep(c("geneA", "geneB"), times = length(genomes)),
      value = seq_len(2 * length(genomes)) %% 3
    ),
    weird = tibble::tibble(
      genome_id = rep(genomes, each = 2),
      widget = rep(c("w1", "w2"), times = length(genomes)),
      value = seq_len(2 * length(genomes)) %% 5
    ),
    decoy_count = tibble::tibble(genome_id = genomes, a = 1, b = 2, c = 3)
  )
}

fxManifestName <- function() "manifest_run_fixture.json"

# A ready ORB in a temporary directory, removed when `env` exits, declaring every file.
fxOrb <- function(env = parent.frame(),
                  metadata = fxMetadata(),
                  parquets = fxDefaultParquets(),
                  artifact_status = "ready",
                  producer_status = "success") {
  dir <- withr::local_tempdir("orb", .local_envir = env)

  metadata_path <- file.path(dir, "metadata.parquet")
  arrow::write_parquet(metadata, metadata_path)

  table_paths <- vapply(names(parquets), function(name) {
    path <- file.path(dir, paste0(name, ".parquet"))
    arrow::write_parquet(parquets[[name]], path)
    path
  }, character(1), USE.NAMES = FALSE)

  duckdb_path <- file.path(dir, "Tst_parquet.duckdb")
  file.create(duckdb_path)

  declared <- c(metadata_path, table_paths, duckdb_path)
  manifest <- fxManifest(
    dir,
    paths = declared, sizes = as.numeric(file.info(declared)$size),
    artifact_status = artifact_status, producer_status = producer_status
  )
  fxWriteManifest(dir, manifest)
  dir
}

# An amRdata manifest with fixed IDs, timestamps and sizes; its paths need not exist.
fxManifest <- function(dir = "/orb",
                       paths = file.path(dir, c(
                         "metadata.parquet", "gene_count.parquet", "Tst_parquet.duckdb"
                       )),
                       sizes = c(100, 200, 0),
                       artifact_status = "ready",
                       producer_status = "success") {
  run_id <- "run_fixture_producer"
  metadata_path <- file.path(dir, "metadata.parquet")
  duckdb_path <- file.path(dir, "Tst_parquet.duckdb")

  list(
    schema_version = 1L,
    manifest_type = "amR_dataset",
    manifest_id = "manifest_run_fixture",
    manifest_created_at = "2026-01-01 00:00:00",
    manifest_updated_at = "2026-01-01 00:00:01",
    dataset_id = "Tst",
    dataset = list(
      duckdb = duckdb_path,
      selection = list(user_bacs = "Testus fixtus")
    ),
    artifacts = list(
      amRml_input = list(
        status = artifact_status,
        updated_at = "2026-01-01 00:00:01",
        producer = "amRdata",
        producer_run_id = run_id,
        directory = dir,
        parquet_duckdb = duckdb_path,
        metadata_parquet = metadata_path
      )
    ),
    runs = list(
      list(
        run_id = run_id,
        status = producer_status,
        started_at = "2026-01-01 00:00:00",
        finished_at = "2026-01-01 00:00:01",
        stages = list(
          list(
            name = "clean_metadata_and_export",
            status = "success",
            outputs = lapply(seq_along(paths), function(i) {
              list(
                path = paths[[i]],
                exists = TRUE,
                size_bytes = sizes[[i]],
                modified_at = "2026-01-01 00:00:01"
              )
            })
          )
        ),
        events = list()
      )
    )
  )
}

fxReadManifest <- function(dir, name = fxManifestName()) {
  jsonlite::read_json(file.path(dir, name), simplifyVector = FALSE)
}

fxWriteManifest <- function(dir, manifest, name = fxManifestName()) {
  jsonlite::write_json(
    manifest,
    file.path(dir, name),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  invisible(dir)
}

# Apply `fn` to the parsed manifest and write it back.
fxEditManifest <- function(dir, fn, name = fxManifestName()) {
  fxWriteManifest(dir, fn(fxReadManifest(dir, name)), name)
}

# Re-record declared file sizes, so a changed file doesn't trip the size check.
fxRecordSizes <- function(dir) {
  fxEditManifest(dir, function(m) {
    outputs <- m$runs[[1]]$stages[[1]]$outputs
    m$runs[[1]]$stages[[1]]$outputs <- lapply(outputs, function(o) {
      o$size_bytes <- as.numeric(file.info(o$path)$size)
      o
    })
    m
  })
}

# Every entry in `dir`, hidden and empty ones included, with each file's checksum.
fxOrbState <- function(dir) {
  entries <- sort(list.files(dir, recursive = TRUE, all.files = TRUE, include.dirs = TRUE))
  paths <- file.path(dir, entries)
  files <- paths[!dir.exists(paths)]
  list(entries = entries, checksums = tools::md5sum(files))
}
