# Tests for eligibleScopes(), in the order it decides what can be modelled.

# Scopes of one mode, without the mode column.
fxScopes <- function(eligibility, mode_id) {
  scopes <- eligibility$scopes
  scopes[scopes$mode_id == mode_id, -1]
}

# Genome IDs kept for one eligible scope and role.
fxMemberIds <- function(eligibility, mode_id, target, train_group, test_group, role) {
  m <- eligibility$members
  keep <- m$mode_id == mode_id & m$target == target & m$role == role &
    m$train_group %in% train_group & m$test_group %in% test_group
  sort(m$genome_id[keep])
}

# The arguments

test_that("anything that isn't a profile is rejected, naming eligibleScopes()", {
  cnd <- rlang::catch_cnd(eligibleScopes(list()), classes = "error")

  expect_s3_class(cnd, "amrml_invalid_argument")
  expect_equal(rlang::call_name(cnd$call), "eligibleScopes")
})

test_that("unknown modes and thresholds that aren't whole numbers are rejected", {
  profile <- profileORB(readORB(fxOrb()))

  expect_error(eligibleScopes(profile, modes = "nope"), class = "amrml_invalid_argument")
  expect_error(eligibleScopes(profile, n_fold = 1), class = "amrml_invalid_argument")
  expect_error(eligibleScopes(profile, min_genomes = 2.5), class = "amrml_invalid_argument")
  expect_error(eligibleScopes(profile, min_genomes = Inf), class = "amrml_invalid_argument")
  expect_error(eligibleScopes(profile, modes = list("vanilla")), class = "amrml_invalid_argument")
})

# Whether each mode can run

test_that("a mode that isn't supported yet doesn't run", {
  eligibility <- fxEligibility(modes = "lodo")

  expect_equal(eligibility$modes$rule_id, "mode_not_supported")
  expect_equal(nrow(eligibility$scopes), 0)
})

test_that("grouped modes need 2 groups, and leave-one-out needs 3", {
  metadata <- fxRows(2, 2, year_bin = c("Y1", "Y1", "Y2", "Y2"))
  eligibility <- fxEligibility(metadata, modes = c("stratified_year", "cross_year", "loto"))
  modes <- eligibility$modes

  expect_equal(modes$feasible, c(TRUE, TRUE, FALSE))
  expect_equal(modes$rule_id, c(NA, NA, "too_few_groups"))
  expect_equal(modes$n_groups, c(2, 2, 2))
  # Feasible modes can still have no eligible scope: Y1 is all Resistant, Y2 all Susceptible.
  expect_equal(modes$n_scopes, c(4, 4, NA))
  expect_equal(modes$n_eligible, c(0, 0, NA))

  one_drug <- fxEligibility(metadata, modes = "cross_drug")$modes
  expect_equal(one_drug$rule_id, "too_few_groups")
})

test_that("without a year column, the year modes have no groups", {
  metadata <- fxRows(10, 10)
  metadata$year_bin <- NULL
  modes <- fxEligibility(metadata, modes = c("stratified_year", "cross_year", "loto"))$modes

  expect_equal(modes$rule_id, rep("too_few_groups", 3))
  expect_equal(modes$n_groups, c(0, 0, 0))
})

test_that("every supported mode finds an eligible scope on an ORB with enough data", {
  eligibility <- fxEligibility(fxCrossed(), min_genomes = 4, n_fold = 2)
  eligible <- eligibility$scopes[eligibility$scopes$eligible & eligibility$scopes$unit == "drug", ]

  expect_setequal(unique(eligible$mode_id), c(
    "vanilla", "stratified_year", "stratified_country", "cross_year", "cross_country",
    "cross_drug", "loto", "logo"
  ))

  modes <- eligibility$modes[eligibility$modes$feasible, ]
  scopes <- eligibility$scopes
  expect_equal(modes$n_scopes, as.vector(table(scopes$mode_id)[modes$mode_id]))
  expect_equal(modes$n_eligible, as.vector(table(scopes$mode_id[scopes$eligible])[modes$mode_id]))
})

# Whether each drug or class can run in a mode

test_that("a drug or class tested in too few of a mode's groups gets one decision, not scopes", {
  # DRA is tested in all three years; DRB, in its own class, only in Y1 and Y2.
  metadata <- rbind(
    fxRows(3, 3, year_bin = rep(c("Y1", "Y2", "Y3"), 2)),
    fxRows(2, 2, drug = "DRB", from = 7, class = "CLY", year_bin = c("Y1", "Y2", "Y1", "Y2"))
  )
  scopes <- fxScopes(fxEligibility(metadata, modes = "loto"), "loto")
  drb <- scopes[scopes$target == "DRB", ]

  expect_equal(nrow(scopes[scopes$target == "DRA", ]), 3)
  expect_equal(nrow(drb), 1)
  expect_equal(drb$rule_id, "too_few_groups_tested")
  expect_equal(drb$n_genomes, 4)
  expect_true(is.na(drb$test_group))
})

test_that("a drug tested only on genomes without a group still gets a decision", {
  # DRB's genomes have no year, so it is in none of the year modes' groups.
  metadata <- rbind(
    fxRows(3, 3, year_bin = rep(c("Y1", "Y2"), 3)),
    fxRows(2, 2, drug = "DRB", from = 7, class = "CLY", year_bin = NA)
  )
  scopes <- fxScopes(fxEligibility(metadata, modes = "stratified_year"), "stratified_year")
  drb <- scopes[scopes$target == "DRB", ]

  expect_equal(nrow(drb), 1)
  expect_equal(drb$rule_id, "too_few_groups_tested")
  expect_equal(drb$n_genomes, 0)
})

test_that("a drug tested in one group gets no stratified or cross scope of its own", {
  # DRB is tested in Y1 only, so its stratified scope would repeat vanilla.
  metadata <- rbind(
    fxRows(2, 2, year_bin = c("Y1", "Y2", "Y1", "Y2")),
    fxRows(2, 2, drug = "DRB", from = 5, class = "CLY", year_bin = "Y1")
  )
  eligibility <- fxEligibility(metadata, modes = c("stratified_year", "cross_year"))
  drb <- eligibility$scopes[eligibility$scopes$target == "DRB", ]

  expect_equal(drb$mode_id, c("stratified_year", "cross_year"))
  expect_equal(drb$rule_id, rep("too_few_groups_tested", 2))
})

# Scopes without a separate test set

test_that("vanilla scopes are the profile's whole-dataset counts", {
  profile <- profileORB(readORB(fxOrb()))
  scopes <- fxScopes(eligibleScopes(profile, modes = "vanilla"), "vanilla")
  whole <- profile$targets[is.na(profile$targets$stratum), ]

  expect_equal(scopes$target, whole$target)
  expect_equal(scopes$n_genomes, whole$n_genomes)
  expect_equal(scopes$n_resistant, whole$n_resistant)
  expect_true(all(is.na(scopes$test_n_genomes)))
})

test_that("stratified scopes are per value, leaving out genomes without one", {
  metadata <- fxRows(2, 2, country = c("CX", "CX", "CY", NA))
  scopes <- fxScopes(fxEligibility(metadata, modes = "stratified_country"), "stratified_country")
  drug <- scopes[scopes$unit == "drug", ]

  expect_equal(drug$train_group, c("CX", "CY"))
  expect_equal(drug$n_genomes, c(2, 1))
})

# Scopes with a separate test set

test_that("cross scopes train on one group and test on another, both ways round", {
  metadata <- fxRows(2, 3, year_bin = c("Y1", "Y2", "Y1", "Y2", "Y2"))
  scopes <- fxScopes(fxEligibility(metadata, modes = "cross_year"), "cross_year")
  drug <- scopes[scopes$unit == "drug", ]

  expect_equal(drug$train_group, c("Y1", "Y2"))
  expect_equal(drug$test_group, c("Y2", "Y1"))
  expect_equal(drug$n_genomes, c(2, 3))
  expect_equal(drug$test_n_genomes, c(3, 2))
  expect_equal(drug$test_n_resistant, c(1, 1))
})

test_that("cross-drug test sets leave out genomes the training drug was trained on", {
  # DRA covers g001-g004 and DRB covers g003-g006, so each tests on the other's two others.
  metadata <- rbind(fxRows(2, 2), fxRows(2, 2, drug = "DRB", from = 3, class = "CLY"))
  scopes <- fxScopes(fxEligibility(metadata, modes = "cross_drug"), "cross_drug")

  expect_equal(scopes$train_group, c("DRA", "DRB"))
  expect_equal(scopes$test_group, c("DRB", "DRA"))
  expect_equal(scopes$n_genomes, c(4, 4))
  expect_equal(scopes$test_n_genomes, c(2, 2))
})

test_that("a cross-drug test set left empty counts zero genomes and fails the test-set rule", {
  # DRB's genomes are all DRA's, so DRA has nothing left to test on.
  metadata <- rbind(fxRows(5, 5), fxRows(2, 2, drug = "DRB", class = "CLY"))
  eligibility <- fxEligibility(metadata, modes = "cross_drug", min_genomes = 2, n_fold = 2)
  scope <- fxScopes(eligibility, "cross_drug")
  scope <- scope[scope$train_group == "DRA", ]

  expect_equal(scope$test_n_genomes, 0)
  expect_equal(scope$rule_id, "test_too_small")
})

test_that("leave-one-out trains on every other group and tests on the held-out one", {
  metadata <- fxRows(3, 3, year_bin = c("Y1", "Y2", "Y3", "Y1", "Y2", "Y3"))
  scopes <- fxScopes(fxEligibility(metadata, modes = "loto"), "loto")
  drug <- scopes[scopes$unit == "drug", ]

  expect_equal(drug$test_group, c("Y1", "Y2", "Y3"))
  expect_true(all(is.na(drug$train_group)))
  expect_equal(drug$n_genomes, c(4, 4, 4))
  expect_equal(drug$test_n_genomes, c(2, 2, 2))
})

test_that("each grouped mode splits genomes by its own stratum", {
  metadata <- fxCrossed()
  eligibility <- fxEligibility(metadata, min_genomes = 4, n_fold = 2)
  dra <- metadata[metadata$drug_abbr == "DRA", ]
  ids <- function(keep) sort(dra$genome.genome_id[keep])

  expect_equal(
    fxMemberIds(eligibility, "stratified_year", "DRA", "Y1", NA, "train"),
    ids(dra$year_bin == "Y1")
  )
  expect_equal(
    fxMemberIds(eligibility, "cross_country", "DRA", "C1", "C2", "train"),
    ids(dra$country_abbr == "C1")
  )
  expect_equal(
    fxMemberIds(eligibility, "cross_country", "DRA", "C1", "C2", "test"),
    ids(dra$country_abbr == "C2")
  )
  expect_equal(
    fxMemberIds(eligibility, "logo", "DRA", NA, "C2", "train"),
    ids(dra$country_abbr != "C2")
  )
  expect_equal(
    fxMemberIds(eligibility, "logo", "DRA", NA, "C2", "test"),
    ids(dra$country_abbr == "C2")
  )
})

# The rules

test_that("a class with data for one drug is a duplicate; one with two drugs is not", {
  # CLX holds DRA only. CLY holds DRB and DRC on different genomes. CLZ holds DRD and DRE
  # on the same genomes with the same labels, so its labels match DRD's but it has two drugs.
  metadata <- rbind(
    fxRows(3, 3),
    fxRows(3, 3, drug = "DRB", class = "CLY"),
    fxRows(3, 3, drug = "DRC", from = 7, class = "CLY"),
    fxRows(3, 3, drug = "DRD", class = "CLZ"),
    fxRows(3, 3, drug = "DRE", class = "CLZ")
  )
  scopes <- fxScopes(
    fxEligibility(metadata, modes = "vanilla", min_genomes = 2, n_fold = 2),
    "vanilla"
  )
  classes <- scopes[scopes$unit == "drug_class", ]

  expect_equal(classes$target, c("CLX", "CLY", "CLZ"))
  expect_equal(classes$rule_id, c("duplicate_class", NA, NA))
})

test_that("a class's drugs are counted within each scope's genomes", {
  # Only DRB has data in Y1, so CLY repeats it there. In Y2, CLY combines DRB and DRC; it
  # then fails only because it has no second eligible year.
  metadata <- rbind(
    fxRows(3, 3, drug = "DRB", class = "CLY", year_bin = "Y1"),
    fxRows(3, 3, drug = "DRB", from = 7, class = "CLY", year_bin = "Y2"),
    fxRows(3, 3, drug = "DRC", from = 7, class = "CLY", year_bin = "Y2")
  )
  scopes <- fxScopes(
    fxEligibility(metadata, modes = "stratified_year", min_genomes = 2, n_fold = 2),
    "stratified_year"
  )
  classes <- scopes[scopes$unit == "drug_class", ]

  expect_equal(classes$train_group, c("Y1", "Y2"))
  expect_equal(classes$rule_id, c("duplicate_class", "too_few_eligible_groups"))
})

test_that("only the first rule a scope fails is recorded", {
  # 3 Resistant genomes: the drug has one phenotype and too few genomes, and the phenotype
  # rule comes first. Its class holds one drug, a rule checked before both.
  scopes <- fxScopes(fxEligibility(fxRows(3, 0), modes = "vanilla"), "vanilla")

  expect_equal(scopes$rule_id[scopes$unit == "drug"], "single_phenotype")
  expect_equal(scopes$rule_id[scopes$unit == "drug_class"], "duplicate_class")
})

test_that("a scope with fewer training genomes than min_genomes is skipped", {
  scope <- function(min_genomes) {
    eligibility <- fxEligibility(fxRows(10, 10), modes = "vanilla", min_genomes = min_genomes, n_fold = 2)
    fxScopes(eligibility, "vanilla")[1, ]
  }

  expect_equal(scope(21)$rule_id, "too_few_genomes")
  expect_true(scope(20)$eligible)
})

test_that("the fold check counts the 20% holdout only when there's no test set", {
  # Internal: the split keeps floor(0.8 x n) of each phenotype, so 6 Resistant leave 4 for
  # 5 folds, too few; 7 leave 5. With 4 folds, 5 Resistant leave exactly 4, enough.
  internal <- function(n_r, n_fold = 5) {
    eligibility <- fxEligibility(fxRows(n_r, 30), modes = "vanilla", min_genomes = 2, n_fold = n_fold)
    fxScopes(eligibility, "vanilla")[1, ]
  }
  expect_equal(internal(5)$rule_id, "too_few_for_cv")
  expect_equal(internal(6)$rule_id, "too_few_for_cv")
  expect_true(internal(7)$eligible)
  expect_true(internal(5, n_fold = 4)$eligible)

  # Cross: no holdout, so 5 Resistant training genomes are enough for 5 folds; 4 are not.
  cross <- function(n_r) {
    metadata <- rbind(
      fxRows(n_r, 30, year_bin = "Y1"),
      fxRows(2, 3, from = 36, year_bin = "Y2")
    )
    scopes <- fxScopes(fxEligibility(metadata, modes = "cross_year", min_genomes = 2), "cross_year")
    scopes[scopes$unit == "drug" & scopes$train_group == "Y1", ]
  }
  expect_true(cross(5)$eligible)
  expect_equal(cross(4)$rule_id, "too_few_for_cv")
})

test_that("a test set needs enough genomes and enough of its rarer phenotype", {
  testSetRule <- function(n_r, n_s) {
    metadata <- rbind(
      fxRows(10, 10, year_bin = "Y1"),
      fxRows(n_r, n_s, from = 21, year_bin = "Y2")
    )
    scopes <- fxScopes(
      fxEligibility(metadata, modes = "cross_year", min_genomes = 2, n_fold = 2),
      "cross_year"
    )
    scopes$rule_id[scopes$unit == "drug" & scopes$train_group == "Y1"]
  }

  # 1 of the rarer phenotype is below 2; 4 genomes are below 5.
  expect_equal(testSetRule(1, 4), "test_too_small")
  expect_equal(testSetRule(2, 2), "test_too_small")
  expect_true(is.na(testSetRule(2, 3)))

  # The same rule for a leave-one-out held-out group: Y3 has 2 genomes.
  metadata <- rbind(
    fxRows(10, 10, year_bin = "Y1"),
    fxRows(10, 10, from = 21, year_bin = "Y2"),
    fxRows(1, 1, from = 41, year_bin = "Y3")
  )
  scopes <- fxScopes(fxEligibility(metadata, modes = "loto", min_genomes = 2, n_fold = 2), "loto")
  drug <- scopes[scopes$unit == "drug", ]
  expect_equal(drug$rule_id, c(NA, NA, "test_too_small"))
})

test_that("a stratified drug needs at least two groups with enough data", {
  # Y1 has enough genomes; Y2 is too small, so Y1 would be the only stratified model.
  metadata <- rbind(
    fxRows(10, 10, year_bin = "Y1"),
    fxRows(1, 1, from = 21, year_bin = "Y2")
  )
  scopes <- fxScopes(
    fxEligibility(metadata, modes = "stratified_year", min_genomes = 4, n_fold = 2),
    "stratified_year"
  )
  drug <- scopes[scopes$unit == "drug", ]

  expect_equal(drug$train_group, c("Y1", "Y2"))
  expect_equal(drug$rule_id, c("too_few_eligible_groups", "too_few_genomes"))
})

# The result

test_that("the result keeps the profile and settings, and every rule is an error code", {
  eligibility <- fxEligibility(n_fold = 3)
  rules <- c(eligibility$modes$rule_id, eligibility$scopes$rule_id)

  expect_s3_class(eligibility, "amr_eligibility")
  expect_s3_class(eligibility$profile, "amr_orb_profile")
  expect_equal(eligibility$settings$n_fold, 3)
  expect_equal(eligibility$settings$holdout, 0.2)
  expect_true(all(stats::na.omit(rules) %in% .RULES$rule_id))
})

test_that("when no mode runs, scopes and members are empty but keep their columns", {
  ran <- fxEligibility(modes = "vanilla")
  none <- fxEligibility(modes = "lodo")

  expect_equal(nrow(none$scopes), 0)
  expect_equal(nrow(none$members), 0)
  expect_identical(none$scopes, ran$scopes[0, ])
  expect_identical(none$members, ran$members[0, ])
})

test_that("the genomes of each eligible scope are kept, with train and test kept apart", {
  # Leave-one-out: hold out Y1 for DRA.
  crossed <- fxCrossed()
  loto <- fxEligibility(crossed, modes = "loto", min_genomes = 4, n_fold = 2)
  dra <- crossed[crossed$drug_abbr == "DRA", ]
  train <- fxMemberIds(loto, "loto", "DRA", NA, "Y1", "train")
  test <- fxMemberIds(loto, "loto", "DRA", NA, "Y1", "test")

  expect_equal(train, sort(dra$genome.genome_id[dra$year_bin != "Y1"]))
  expect_equal(test, sort(dra$genome.genome_id[dra$year_bin == "Y1"]))
  expect_length(intersect(train, test), 0)

  # Cross-drug: DRA covers g001-g004 and DRB g003-g006, so DRA -> DRB tests on g005 and g006.
  metadata <- rbind(
    fxRows(2, 2),
    fxMetadataRows(
      genome = c("g003", "g004", "g005", "g006"), drug = "DRB", class = "CLY",
      phenotype = c("Resistant", "Susceptible", "Resistant", "Susceptible")
    )
  )
  cross <- fxEligibility(
    metadata,
    modes = "cross_drug", min_genomes = 2, n_fold = 2, min_test_genomes = 2, min_test_minority = 1
  )

  expect_equal(fxMemberIds(cross, "cross_drug", "DRA", "DRA", "DRB", "train"), sprintf("g%03d", 1:4))
  expect_equal(fxMemberIds(cross, "cross_drug", "DRA", "DRA", "DRB", "test"), c("g005", "g006"))
  # DRB -> DRA would test on g001 and g002, both Resistant, so it fails the test-set rule
  # and keeps no genomes.
  scopes <- fxScopes(cross, "cross_drug")
  expect_equal(scopes$rule_id[scopes$train_group == "DRB"], "test_too_small")
  expect_false(any(cross$members$train_group == "DRB"))
})

test_that("the print shows each mode's result", {
  eligibility <- fxEligibility(modes = c("vanilla", "stratified_year", "lodo"))

  expect_output(print(eligibility), "40 training genomes, 5 folds")
  expect_output(print(eligibility), "vanilla +: 0 of 3 eligible")
  expect_output(print(eligibility), "stratified_year +: not run \\(too_few_groups\\)")
  expect_output(print(eligibility), "lodo +: not run \\(mode_not_supported\\)")
})
