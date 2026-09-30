# Tests for profileORB(), in the order it profiles an ORB.

fxProfile <- function(metadata = fxMetadata(), parquets = fxDefaultParquets(),
                      env = parent.frame()) {
  profileORB(readORB(fxOrb(env = env, metadata = metadata, parquets = parquets)))
}

# Whole-dataset counts for one target.
fxTarget <- function(profile, unit, target) {
  t <- profile$targets
  t[t$unit == unit & t$target == target & is.na(t$stratum), , drop = FALSE]
}

# Parquets whose feature tables cover only `genomes`.
fxParquetsFor <- function(genomes) {
  lapply(fxDefaultParquets(), function(x) x[x$genome_id %in% genomes, , drop = FALSE])
}

# The argument

test_that("anything that isn't an ORB is rejected, naming profileORB()", {
  cnd <- rlang::catch_cnd(profileORB(list()), classes = "error")

  expect_s3_class(cnd, "amrml_invalid_argument")
  expect_equal(rlang::call_name(cnd$call), "profileORB")
})

# Reading the metadata

test_that("metadata without a required column is rejected", {
  metadata <- fxMetadata()
  metadata$class_abbr <- NULL

  expect_error(fxProfile(metadata), class = "amrml_metadata_columns_missing")
})

# Rows that give a genome a phenotype

test_that("only Resistant and Susceptible rows count; genomes without one are excluded", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g2", "g3", "g4"),
    drug = "DRA",
    phenotype = c("Resistant", "Intermediate", "", NA)
  )
  expect_message(profile <- fxProfile(metadata), "3 genomes are excluded")

  expect_equal(unique(profile$phenotypes$genome_id), "g1")
  expect_equal(profile$excluded_genomes$genome_id, c("g2", "g3", "g4"))
  expect_equal(unique(profile$excluded_genomes$reason), "no_usable_phenotype")
})

# Drug phenotypes and conflicts

test_that("each genome gets its own phenotype per drug", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g1", "g2", "g2"),
    drug = c("DRA", "DRB", "DRA", "DRB"),
    phenotype = c("Resistant", "Susceptible", "Susceptible", "Resistant")
  )
  profile <- fxProfile(metadata)
  drugs <- profile$phenotypes[profile$phenotypes$unit == "drug", ]

  expect_equal(drugs, tibble::tibble(
    genome_id = c("g1", "g1", "g2", "g2"),
    unit = "drug",
    target = c("DRA", "DRB", "DRA", "DRB"),
    phenotype = c("Resistant", "Susceptible", "Susceptible", "Resistant")
  ))
})

test_that("a genome recorded as both Resistant and Susceptible to a drug gets no phenotype", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g1", "g1", "g2"),
    drug = c("DRA", "DRA", "DRB", "DRA"),
    phenotype = c("Resistant", "Susceptible", "Susceptible", "Susceptible")
  )
  expect_message(profile <- fxProfile(metadata), "1 genome-drug pairs")

  expect_equal(profile$conflicts, tibble::tibble(genome_id = "g1", drug_abbr = "DRA"))
  expect_equal(fxTarget(profile, "drug", "DRA")$n_genomes, 1)
  expect_equal(fxTarget(profile, "drug", "DRB")$n_genomes, 1)
})

# Genomes missing from a feature table

test_that("genomes missing from a feature table are excluded, with the tables lacking them", {
  parquets <- fxDefaultParquets()
  parquets$gene_count <- parquets$gene_count[parquets$gene_count$genome_id != "g5", ]

  expect_message(profile <- fxProfile(parquets = parquets), "missing_features: 1")

  expect_equal(profile$genomes$genome_id, c("g1", "g2", "g3", "g4"))
  expect_false("g5" %in% profile$phenotypes$genome_id)
  expect_equal(profile$excluded_genomes$reason, "missing_features")
  expect_equal(profile$excluded_genomes$missing_from, "gene")
})

test_that("conflicts are recorded for genomes excluded for missing features", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g2", "g2"),
    drug = "DRA",
    phenotype = c("Resistant", "Resistant", "Susceptible")
  )

  expect_message(
    expect_message(profile <- fxProfile(metadata, fxParquetsFor("g1")), "missing_features: 1"),
    "1 genome-drug pairs"
  )
  expect_equal(profile$conflicts$genome_id, "g2")
})

# Class phenotypes

test_that("a class is Resistant if any member drug is, Susceptible only if all are", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g1", "g2", "g2", "g3"),
    drug = c("DRA", "DRB", "DRA", "DRB", "DRA"),
    phenotype = c("Resistant", "Susceptible", "Susceptible", "Susceptible", "Susceptible")
  )
  profile <- fxProfile(metadata)
  classes <- profile$phenotypes[profile$phenotypes$unit == "drug_class", ]

  expect_equal(classes$genome_id, c("g1", "g2", "g3"))
  expect_equal(classes$phenotype, c("Resistant", "Susceptible", "Susceptible"))
})

test_that("a contradictory member drug stops a class being Susceptible", {
  # g1's DRA might have been Resistant; g2 is Resistant through DRB whatever DRA was.
  metadata <- fxMetadataRows(
    genome = c("g1", "g1", "g1", "g2", "g2", "g2"),
    drug = c("DRA", "DRA", "DRB", "DRA", "DRA", "DRB"),
    phenotype = c(
      "Resistant", "Susceptible", "Susceptible",
      "Resistant", "Susceptible", "Resistant"
    )
  )
  expect_message(profile <- fxProfile(metadata), "2 genome-drug pairs")
  classes <- profile$phenotypes[profile$phenotypes$unit == "drug_class", ]

  expect_equal(classes$genome_id, "g2")
  expect_equal(classes$phenotype, "Resistant")
})

test_that("a drug with no class keeps its drug phenotype but adds no class phenotype", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g1"),
    drug = c("DRA", "DRZ"),
    phenotype = "Resistant",
    class = c("CLX", NA)
  )
  profile <- fxProfile(metadata)

  expect_equal(fxTarget(profile, "drug", "DRZ")$n_resistant, 1)
  expect_equal(unique(profile$targets$target[profile$targets$unit == "drug_class"]), "CLX")
})

# Excluded genomes

test_that("every genome is either profiled or excluded with one reason", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g2", "g3", "g3", "g4"),
    drug = "DRA",
    phenotype = c("Resistant", "Intermediate", "Resistant", "Susceptible", "Susceptible")
  )

  expect_message(
    expect_message(
      profile <- fxProfile(metadata, fxParquetsFor(c("g1", "g2", "g3"))),
      "3 genomes are excluded"
    ),
    "both Resistant and Susceptible"
  )

  excluded <- profile$excluded_genomes
  expect_equal(profile$genomes$genome_id, "g1")
  expect_equal(excluded$genome_id, c("g2", "g3", "g4"))
  expect_equal(excluded$reason, c("no_usable_phenotype", "contradictory", "missing_features"))
  expect_equal(excluded$missing_from, c(NA, NA, "gene, widget"))
})

test_that("an ORB with no genome left to profile is an error that lists the reasons", {
  only_intermediate <- fxMetadataRows(genome = "g1", drug = "DRA", phenotype = "Intermediate")
  err <- expect_error(fxProfile(only_intermediate), class = "amrml_no_labelled_genomes")
  expect_equal(err$observed$excluded$reason, "no_usable_phenotype")

  only_contradictory <- fxMetadataRows(
    genome = "g1", drug = "DRA", phenotype = c("Resistant", "Susceptible")
  )
  err <- expect_error(fxProfile(only_contradictory), class = "amrml_no_labelled_genomes")
  expect_equal(err$observed$excluded$reason, "contradictory")
})

# Genomes and strata

test_that("a genome with two values for a stratum is rejected", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g1"),
    drug = c("DRA", "DRB"),
    phenotype = "Resistant",
    country = c("CX", "CY")
  )

  err <- expect_error(fxProfile(metadata), class = "amrml_stratum_inconsistent")
  expect_equal(err$observed$genomes, "g1")
})

test_that("profiling counts by whatever strata the registry lists", {
  local_mocked_bindings(.STRATA = c(.STRATA, host = "genome.host_common_name"))
  metadata <- fxMetadataRows(genome = c("g1", "g2"), drug = "DRA", phenotype = "Resistant")
  metadata$genome.host_common_name <- c("Human", "Cow")
  profile <- fxProfile(metadata)

  by_host <- profile$targets[profile$targets$stratum %in% "host", ]
  expect_equal(profile$genomes$host, c("Human", "Cow"))
  expect_setequal(by_host$stratum_value, c("Human", "Cow"))
})

test_that("a stratum column the metadata lacks is left out, not an error", {
  metadata <- fxMetadata()
  metadata$year_bin <- NULL
  profile <- fxProfile(metadata)

  expect_false("year" %in% names(profile$genomes))
  expect_equal(unique(stats::na.omit(profile$targets$stratum)), "country")
})

# Counts

test_that("strata count observed values; genomes without one count in the whole dataset only", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g2", "g3", "g4"),
    drug = "DRA",
    phenotype = "Resistant",
    country = c("CX", "CY", NA, ""),
    year_bin = factor("2000-2004", levels = c("1995-1999", "2000-2004", "2005-2009"))
  )
  profile <- fxProfile(metadata)
  targets <- profile$targets
  by_country <- targets[targets$stratum %in% "country" & targets$unit == "drug", ]

  expect_type(profile$genomes$year, "character")
  expect_equal(unique(targets$stratum_value[targets$stratum %in% "year"]), "2000-2004")
  expect_equal(by_country$stratum_value, c("CX", "CY"))
  expect_equal(by_country$n_genomes, c(1, 1))
  expect_equal(fxTarget(profile, "drug", "DRA")$n_genomes, 4)
})

# The profile object

test_that("the profile has its parts and prints a summary", {
  profile <- fxProfile()

  expect_s3_class(profile, "amr_orb_profile")
  expect_named(
    profile,
    c("orb", "genomes", "phenotypes", "targets", "conflicts", "excluded_genomes")
  )
  expect_output(print(profile), "5 profiled, 0 excluded")
  expect_output(print(profile), "country \\(values: 1, missing: 0\\)")
  expect_output(print(profile), "DRA:5/0, DRB:0/5")
  expect_output(print(profile), "CLX:5/0")
})

test_that("profiling does not write to the ORB", {
  dir <- fxOrb()
  orb <- readORB(dir)
  before <- fxOrbState(dir)

  profileORB(orb)

  expect_equal(fxOrbState(dir), before)
})
