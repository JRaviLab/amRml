# Tests for fitting the MDR combination target, whose labels aren't Resistant and Susceptible.

# `n` genomes per label, tested on DRA (class CLX) and DRB (class CLY): resistant to CLX only,
# to both, or to neither. geneX marks CLX resistance, geneY CLY resistance, and geneN alternates
# with no signal.
fxLabelled <- function(n) {
  genomes <- sprintf("g%03d", seq_len(3 * n))
  label <- rep(c("x", "xy", "s"), each = n)
  metadata <- fxMetadataRows(
    genome = rep(genomes, each = 2),
    drug = rep(c("DRA", "DRB"), 3 * n),
    class = rep(c("CLX", "CLY"), 3 * n),
    phenotype = ifelse(
      c(rbind(label %in% c("x", "xy"), label == "xy")), "Resistant", "Susceptible"
    )
  )
  genes <- tibble::tibble(
    genome_id = rep(genomes, 3),
    gene = rep(c("geneX", "geneY", "geneN"), each = 3 * n),
    value = c(
      as.integer(label %in% c("x", "xy")), as.integer(label == "xy"), rep_len(0:1, 3 * n)
    )
  )
  list(metadata = metadata, parquets = list(gene_count = genes))
}

# The combination target's matrices from `fxLabelled(n)`, built and fitted.
fxCombinationFit <- function(n = 20, env = parent.frame()) {
  data <- fxLabelled(n)
  eligibility <- fxEligibility(
    data$metadata, modes = "vanilla", parquets = data$parquets,
    min_genomes = n, n_fold = 2, min_label_genomes = n, env = env
  )
  eligibility$scopes <- eligibility$scopes[eligibility$scopes$target == "resistant_classes", ]
  eligibility$members <- eligibility$members[eligibility$members$target == "resistant_classes", ]
  plan <- planMatrices(eligibility)
  built <- buildMatrices(plan, withr::local_tempdir(.local_envir = env))
  fitModels(expandTasks(built), BPPARAM = BiocParallel::SerialParam())
}

test_that("folds are assigned for every label, not only Resistant and Susceptible", {
  phenotype <- factor(rep(c("CLX", "CLX_CLY", "Susceptible"), each = 4))
  fold <- withr::with_seed(1, .stratifiedFolds(phenotype, 2))

  expect_true(all(fold %in% 1:2))
  expect_equal(as.vector(table(phenotype, fold)), rep(2, 6))
})

test_that("Resistant and Susceptible fold assignment is unchanged", {
  phenotype <- factor(rep(c("Resistant", "Susceptible"), c(5, 7)), levels = .PHENOTYPES)
  fold <- withr::with_seed(1, .stratifiedFolds(phenotype, 3))

  expected <- withr::with_seed(1, {
    expected <- integer(12)
    expected[1:5] <- sample(rep_len(1:3, 5))
    expected[6:12] <- sample(rep_len(1:3, 7))
    expected
  })
  expect_equal(fold, expected)
})

test_that("the combination target is fitted by multinomial regression, one task per matrix", {
  fit <- fxCombinationFit()

  expect_equal(fit$tasks$status, rep("fitted", nrow(fit$tasks)))
  expect_true(all(is.na(fit$tasks$rule_id)))
  expect_equal(unique(fit$performance$n_labels), 3)
  expect_false(anyNA(fit$performance$mcc))
  expect_true(all(is.na(fit$performance$res_prop)))
  expect_true(all(is.na(fit$performance$auprc)))
})

test_that("importance is per label, and predictions carry the predicted label's probability", {
  fit <- fxCombinationFit()

  expect_true(all(fit$importance$label %in% c("CLX", "CLX_CLY", "Susceptible")))
  expect_true("geneX" %in% fit$importance$Variable)
  expect_true(all(is.na(fit$predictions$prob_resistant)))
  expect_true(all(fit$predictions$prob_predicted >= 1 / 3 & fit$predictions$prob_predicted <= 1))
  expect_true(all(fit$predictions$truth %in% c("CLX", "CLX_CLY", "Susceptible")))
})
