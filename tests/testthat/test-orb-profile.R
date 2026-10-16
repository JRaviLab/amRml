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

# Drug phenotypes

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

test_that("repeated rows for a genome and drug count once if they agree", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g1", "g2"),
    drug = "DRA",
    phenotype = c("Resistant", "Resistant", "Susceptible")
  )
  profile <- fxProfile(metadata)

  expect_equal(fxTarget(profile, "drug", "DRA")$n_genomes, 2)
  expect_equal(fxTarget(profile, "drug", "DRA")$n_resistant, 1)
})

test_that("a genome recorded as both Resistant and Susceptible for a drug is an error", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g1", "g2"),
    drug = "DRA",
    phenotype = c("Resistant", "Susceptible", "Susceptible")
  )

  err <- expect_error(fxProfile(metadata), class = "amrml_phenotype_conflict")
  expect_equal(err$observed$conflicts, tibble::tibble(genome_id = "g1", drug_abbr = "DRA"))
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

# Excluded genomes

test_that("every genome is either profiled or excluded with one reason", {
  metadata <- fxMetadataRows(
    genome = c("g1", "g2", "g3"),
    drug = "DRA",
    phenotype = c("Resistant", "Intermediate", "Susceptible")
  )

  expect_message(
    profile <- fxProfile(metadata, fxParquetsFor(c("g1", "g2"))),
    "2 genomes are excluded"
  )

  excluded <- profile$excluded_genomes
  expect_equal(profile$genomes$genome_id, "g1")
  expect_equal(excluded$genome_id, c("g2", "g3"))
  expect_equal(excluded$reason, c("no_usable_phenotype", "missing_features"))
  expect_equal(excluded$missing_from, c(NA, "gene, widget"))
})

test_that("an ORB with no genome left to profile is an error that lists the reasons", {
  metadata <- fxMetadataRows(genome = "g1", drug = "DRA", phenotype = "Intermediate")

  err <- expect_error(fxProfile(metadata), class = "amrml_no_labelled_genomes")
  expect_equal(err$observed$excluded$reason, "no_usable_phenotype")
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
  expect_equal(profile$drug_classes, tibble::tibble(drug = "DRA", class = "CLX"))
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

# MDR

# g1: Resistant in 3 classes. g2: Resistant in 2. g3: tested in 3, Resistant in 1.
# g4: tested in 1 class only. g5: tested in 2 drugs of one class only.
fxMdrMetadata <- function() {
  fxMetadataRows(
    genome = c("g1", "g1", "g1", "g2", "g2", "g2", "g3", "g3", "g3", "g4", "g5", "g5"),
    drug = c("DRA", "DRB", "DRC", "DRA", "DRB", "DRC", "DRA", "DRB", "DRC", "DRA", "DRA", "DRD"),
    phenotype = c(
      "Resistant", "Resistant", "Resistant",
      "Resistant", "Resistant", "Susceptible",
      "Resistant", "Susceptible", "Susceptible",
      "Resistant",
      "Resistant", "Resistant"
    ),
    class = c(
      "CLX", "CLY", "CLZ", "CLX", "CLY", "CLZ", "CLX", "CLY", "CLZ", "CLX", "CLX", "CLX"
    )
  )
}

test_that("MDR status counts classes tested and resistant per genome", {
  profile <- fxProfile(fxMdrMetadata(), parquets = fxDefaultParquets())

  expect_equal(profile$mdr, tibble::tibble(
    genome_id = c("g1", "g2", "g3", "g4", "g5"),
    n_classes_tested = c(3L, 3L, 3L, 1L, 1L),
    n_classes_resistant = c(3L, 2L, 1L, 1L, 1L),
    mdr = c("mdr", "mdr", "not_mdr", "undetermined", "undetermined")
  ))
})

test_that("undetermined genomes get no mdr phenotype; the rest are Resistant or Susceptible", {
  profile <- fxProfile(fxMdrMetadata(), parquets = fxDefaultParquets())
  mdr <- profile$phenotypes[profile$phenotypes$unit == "mdr", ]

  expect_equal(mdr, tibble::tibble(
    genome_id = c("g1", "g2", "g3"),
    unit = "mdr",
    target = "mdr",
    phenotype = c("Resistant", "Resistant", "Susceptible")
  ))
  expect_equal(fxTarget(profile, "mdr", "mdr")$n_resistant, 2)
  expect_equal(fxTarget(profile, "mdr", "mdr")$n_susceptible, 1)
  expect_output(print(profile), "2 MDR, 1 not MDR, 2 undetermined")
})

test_that("a higher mdr_min_classes needs more classes tested and resistant", {
  orb <- readORB(fxOrb(metadata = fxMdrMetadata()))
  profile <- profileORB(orb, mdr_min_classes = 3)

  expect_equal(profile$mdr$mdr, c("mdr", "not_mdr", "not_mdr", "undetermined", "undetermined"))
  expect_equal(profile$settings$mdr_min_classes, 3)
})

test_that("mdr_min_classes must be a whole number of at least 2", {
  orb <- readORB(fxOrb())

  expect_error(profileORB(orb, mdr_min_classes = 1), class = "amrml_invalid_argument")
  expect_error(profileORB(orb, mdr_min_classes = 2.5), class = "amrml_invalid_argument")
})

# The profile object

test_that("the profile has its parts and prints a summary", {
  profile <- fxProfile()

  expect_s3_class(profile, "amr_orb_profile")
  expect_named(
    profile,
    c(
      "orb", "genomes", "phenotypes", "drug_classes", "mdr", "settings", "targets",
      "excluded_genomes"
    )
  )
  expect_output(print(profile), "5 profiled, 0 excluded")
  expect_output(print(profile), "country \\(values: 1, missing: 0\\)")
  expect_output(print(profile), "DRA:5/0, DRB:0/5")
  expect_output(print(profile), "CLX:5/0")
  expect_output(print(profile), "MDR \\(2\\+ classes\\): 0 MDR, 0 not MDR, 5 undetermined")
})

test_that("profiling does not write to the ORB", {
  dir <- fxOrb()
  orb <- readORB(dir)
  before <- fxOrbState(dir)

  profileORB(orb)

  expect_equal(fxOrbState(dir), before)
})
