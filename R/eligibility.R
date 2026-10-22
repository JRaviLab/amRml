# Deciding which modes can run on an ORB, and which drugs and classes in each have enough data.

# Share held out for testing when a mode has no test set; kept in settings for fitting.
.CV_HOLDOUT <- 0.2
# The columns that identify a scope.
.SCOPE_KEYS <- c("unit", "target", "train_group", "test_group")


# Fewest groups a mode needs: leave-one-out needs a held-out group and two to train on.
.groupsNeeded <- function(loo) {
  ifelse(loo, 3L, 2L)
}

# For each mode: whether it can run on this ORB, with the number of groups found.
.modeFeasibility <- function(modes, profile) {
  drugs <- profile$phenotypes$target[profile$phenotypes$unit == "drug"]

  n_groups <- vapply(modes$grouping, function(grouping) {
    if (grouping == "none") {
      NA_integer_
    } else if (grouping == "drug") {
      dplyr::n_distinct(drugs)
    } else {
      dplyr::n_distinct(profile$genomes[[grouping]], na.rm = TRUE)
    }
  }, integer(1))

  rule_id <- rep(NA_character_, nrow(modes))
  rule_id[modes$grouping != "none" & n_groups < .groupsNeeded(modes$LOO)] <- "too_few_groups"
  rule_id[!modes$supported] <- "mode_not_supported"

  tibble::tibble(
    mode_id = modes$mode_id,
    feasible = is.na(rule_id),
    rule_id = rule_id,
    n_groups = unname(n_groups)
  )
}

# Each phenotype row with its genome's value of `grouping`; rows without one are dropped.
.withGroup <- function(phenotypes, genomes, grouping) {
  if (grouping == "none") {
    phenotypes$group <- NA_character_
    return(phenotypes)
  }

  phenotypes$group <- genomes[[grouping]][match(phenotypes$genome_id, genomes$genome_id)]
  phenotypes[!is.na(phenotypes$group), , drop = FALSE]
}

# A builder's training and test genomes, stacked, with each row's role.
.asMembers <- function(train, test = NULL) {
  columns <- c(.SCOPE_KEYS, "genome_id", "phenotype")
  members <- dplyr::bind_rows(
    tibble::tibble(train[columns], role = "train"),
    if (!is.null(test)) tibble::tibble(test[columns], role = "test")
  )
  members[c(.SCOPE_KEYS, "role", "genome_id", "phenotype")]
}

# Internal scopes: each target within each group, trained and tested by cross-validation.
.internalScopes <- function(grouped) {
  grouped$train_group <- grouped$group
  grouped$test_group <- NA_character_
  .asMembers(grouped)
}

# Cross scopes: each target, trained on one group and tested on another, both ways round.
.crossScopes <- function(grouped) {
  groups <- dplyr::distinct(grouped[c("unit", "target", "group")])
  pairs <- dplyr::inner_join(
    groups, groups,
    by = c("unit", "target"), suffix = c("_train", "_test"), relationship = "many-to-many"
  )
  pairs <- pairs[pairs$group_train != pairs$group_test, , drop = FALSE]
  names(pairs) <- c("unit", "target", "train_group", "test_group")

  train <- dplyr::inner_join(
    pairs, grouped,
    by = c("unit", "target", train_group = "group"), relationship = "many-to-many"
  )
  test <- dplyr::inner_join(
    pairs, grouped,
    by = c("unit", "target", test_group = "group"), relationship = "many-to-many"
  )
  .asMembers(train, test)
}

# Cross-drug scopes: trained on one drug, tested on another's genomes the first never saw.
.crossDrugScopes <- function(phenotypes) {
  drugs <- phenotypes[phenotypes$unit == "drug", c("target", "genome_id", "phenotype")]
  drug_ids <- unique(drugs$target)
  pairs <- expand.grid(train_group = drug_ids, test_group = drug_ids, stringsAsFactors = FALSE)
  pairs <- tibble::as_tibble(pairs[pairs$train_group != pairs$test_group, , drop = FALSE])
  pairs$unit <- "drug"
  pairs$target <- pairs$train_group

  train <- dplyr::inner_join(pairs, drugs, by = "target", relationship = "many-to-many")
  test <- dplyr::inner_join(
    pairs, drugs,
    by = c(test_group = "target"), relationship = "many-to-many"
  )
  test <- dplyr::anti_join(test, drugs, by = c("target", "genome_id"))
  .asMembers(train, test)
}

# Leave-one-out scopes: trained on every other group (`train_group` NA), tested on the held-out one.
.looScopes <- function(grouped) {
  held_out <- dplyr::distinct(grouped[c("unit", "target", "group")])
  names(held_out)[[3]] <- "test_group"
  held_out$train_group <- NA_character_

  train <- dplyr::inner_join(
    held_out, grouped,
    by = c("unit", "target"), relationship = "many-to-many"
  )
  train <- train[train$group != train$test_group, , drop = FALSE]
  test <- dplyr::inner_join(
    held_out, grouped,
    by = c("unit", "target", test_group = "group"), relationship = "many-to-many"
  )
  .asMembers(train, test)
}

# Each target tested in fewer of the mode's groups than it needs, with its genome counts in
# them. A target whose genomes all lack a group has none, and is listed too.
.targetsInTooFewGroups <- function(phenotypes, grouped, loo) {
  counts <- dplyr::summarise(
    dplyr::group_by(grouped, .data$unit, .data$target),
    n_groups = dplyr::n_distinct(.data$group),
    n_genomes = dplyr::n(),
    n_resistant = sum(.data$phenotype == "Resistant"),
    n_susceptible = sum(.data$phenotype == "Susceptible"),
    .groups = "drop"
  )
  targets <- dplyr::left_join(
    dplyr::distinct(phenotypes[c("unit", "target")]), counts,
    by = c("unit", "target")
  )
  columns <- c("n_groups", "n_genomes", "n_resistant", "n_susceptible")
  targets[columns] <- lapply(targets[columns], function(n) replace(n, is.na(n), 0L))

  targets[targets$n_groups < .groupsNeeded(loo), , drop = FALSE]
}

# A mode's `members`, one row per scope, genome and role; `too_few`, the targets without scopes.
.buildScopes <- function(mode, profile) {
  if (mode$grouping == "drug") {
    if (!mode$cross_test) {
      rlang::abort(
        paste0("No scope builder for mode '", mode$mode_id, "'."),
        class = c("amrml_internal_error", "amrml_error")
      )
    }
    return(list(members = .crossDrugScopes(profile$phenotypes), too_few = NULL))
  }

  grouped <- .withGroup(profile$phenotypes, profile$genomes, mode$grouping)
  too_few <- NULL

  if (mode$grouping != "none") {
    too_few <- .targetsInTooFewGroups(profile$phenotypes, grouped, mode$LOO)
    grouped <- dplyr::anti_join(grouped, too_few, by = c("unit", "target"))
  }

  members <- if (mode$LOO) {
    .looScopes(grouped)
  } else if (mode$cross_test) {
    .crossScopes(grouped)
  } else {
    .internalScopes(grouped)
  }

  list(members = members, too_few = too_few)
}

# Per scope, the training and test genome counts.
.scopeCounts <- function(members, has_test) {
  counts <- dplyr::summarise(
    dplyr::group_by(members, dplyr::across(dplyr::all_of(c(.SCOPE_KEYS, "role")))),
    n_genomes = dplyr::n(),
    n_resistant = sum(.data$phenotype == "Resistant"),
    n_susceptible = sum(.data$phenotype == "Susceptible"),
    .groups = "drop"
  )

  sides <- c("n_genomes", "n_resistant", "n_susceptible")
  train <- counts[counts$role == "train", c(.SCOPE_KEYS, sides)]
  test <- counts[counts$role == "test", c(.SCOPE_KEYS, sides)]
  names(test) <- c(.SCOPE_KEYS, paste0("test_", sides))
  scopes <- dplyr::full_join(train, test, by = .SCOPE_KEYS)

  # A cross-drug test set is empty when the training drug saw all the test drug's genomes.
  if (has_test) {
    test_counts <- paste0("test_", sides)
    scopes[test_counts] <- lapply(scopes[test_counts], function(n) replace(n, is.na(n), 0L))
  }

  scopes
}

# For each class scope, how many of the class's drugs have data among its genomes.
.classDrugCounts <- function(members, profile) {
  drugs <- profile$phenotypes[profile$phenotypes$unit == "drug", , drop = FALSE]
  drugs <- tibble::tibble(
    genome_id = drugs$genome_id,
    drug = drugs$target,
    class = profile$drug_classes$class[match(drugs$target, profile$drug_classes$drug)]
  )

  tested <- dplyr::inner_join(
    members[members$unit == "drug_class", c(.SCOPE_KEYS, "genome_id")], drugs,
    by = c(target = "class", "genome_id"), relationship = "many-to-many"
  )
  dplyr::summarise(
    dplyr::group_by(tested, dplyr::across(dplyr::all_of(.SCOPE_KEYS))),
    n_drugs = dplyr::n_distinct(.data$drug),
    .groups = "drop"
  )
}

# Each scope's first failed rule, or NA if it passes them all.
.firstFailedRule <- function(scopes, mode, settings, has_test) {
  # A class with data for only one of its drugs is that drug's model again.
  duplicate <- scopes$unit == "drug_class" & scopes$n_drugs %in% 1

  rarer <- pmin(scopes$n_resistant, scopes$n_susceptible)
  enough_for_cv <- if (has_test) {
    rarer >= settings$n_fold
  } else {
    rarer * (1 - settings$holdout) >= settings$n_fold
  }
  test_ok <- !has_test | (
    scopes$test_n_genomes >= settings$min_test_genomes &
      pmin(scopes$test_n_resistant, scopes$test_n_susceptible) >= settings$min_test_minority
  )

  # Check the rules in order and keep each scope's first failure.
  failed <- list(
    duplicate_class = duplicate,
    single_phenotype = rarer == 0,
    too_few_genomes = scopes$n_genomes < settings$min_genomes,
    too_few_for_cv = !enough_for_cv,
    test_too_small = !test_ok
  )
  rule_id <- rep(NA_character_, nrow(scopes))
  for (rule in names(failed)) {
    rule_id[is.na(rule_id) & failed[[rule]]] <- rule
  }

  # A stratified model only means something next to another group's, so a target keeps
  # its eligible stratified scopes only if at least two of its groups have one.
  if (!has_test && mode$grouping != "none") {
    eligible <- is.na(rule_id)
    n_eligible <- stats::ave(as.integer(eligible), scopes$unit, scopes$target, FUN = sum)
    rule_id[eligible & n_eligible < 2] <- "too_few_eligible_groups"
  }

  rule_id
}

# One decision per scope: its genome counts, and the first rule it fails. Targets tested
# in too few groups get one decision each.
.judgeScopes <- function(built, mode, profile, settings) {
  counts <- c(
    "n_genomes", "n_resistant", "n_susceptible",
    "test_n_genomes", "test_n_resistant", "test_n_susceptible"
  )
  has_test <- mode$LOO || mode$cross_test
  scopes <- .scopeCounts(built$members, has_test)
  scopes <- dplyr::left_join(scopes, .classDrugCounts(built$members, profile), by = .SCOPE_KEYS)
  scopes$rule_id <- .firstFailedRule(scopes, mode, settings, has_test)
  decisions <- scopes[c(.SCOPE_KEYS, counts, "rule_id")]

  too_few <- built$too_few
  if (NROW(too_few)) {
    too_few[c("train_group", "test_group")] <- NA_character_
    too_few$rule_id <- "too_few_groups_tested"
    decisions <- dplyr::bind_rows(decisions, too_few[names(too_few) != "n_groups"])
  }
  decisions <- decisions[do.call(order, unname(decisions[.SCOPE_KEYS])), , drop = FALSE]

  tibble::tibble(
    mode_id = mode$mode_id,
    decisions[c(.SCOPE_KEYS, counts)],
    eligible = is.na(decisions$rule_id),
    rule_id = decisions$rule_id
  )
}

#' Decide what can be modeled
#'
#' Compares the modeling modes with a profiled ORB. For each mode, records
#' whether it can run on this ORB at all; for each drug or class in a mode that
#' can, records whether it has enough data to model, or why not. Nothing is
#' built or fitted.
#'
#' @details
#' Counts are of genomes, each with one phenotype. A drug or class needs 2 of a
#' mode's groups (3 for leave-one-out) to get scopes, and a stratified one needs
#' 2 eligible groups. A scope fails, in order, if it's a class with data for only
#' one of its drugs, has one phenotype, has fewer than `min_genomes` genomes,
#' can't fill `n_fold` folds with its rarer phenotype (after a 20% holdout when
#' there's no test set), or has too small a test set.
#'
#' @section Errors:
#' All errors have class `amrml_error`, plus `amrml_invalid_argument`:
#' `profile` is not from [profileORB()], `modes` names an unknown mode, or a
#' threshold is not a finite whole number of at least 1 (2 for `n_fold`).
#'
#' @param profile An `amr_orb_profile` from [profileORB()].
#' @param modes Mode IDs from [modelingModes()], or `NULL` for all of them.
#' @param min_genomes Fewest genomes a scope needs, before any holdout.
#' @param n_fold Number of cross-validation folds.
#' @param min_test_genomes Fewest genomes a separate test set needs.
#' @param min_test_minority Fewest of the rarer phenotype a test set needs.
#'
#' @return An `amr_eligibility` list:
#'   * `profile` and `settings` (the thresholds, and the `holdout` share).
#'   * `modes`: per mode, `feasible`, `rule_id`, `n_groups`, and how many
#'     scopes it has and how many are eligible (`n_scopes`, `n_eligible`).
#'   * `scopes`: per scope, training and test counts (`NA` without a test set),
#'     `eligible` and `rule_id`.
#'   * `members`: each eligible scope's genomes, by `role` (train or test).
#'
#' @examples
#' \dontrun{
#' eligibility <- eligibleScopes(profileORB(readORB("path/to/amRdata/output")))
#' eligibility$scopes
#' }
#' @export
eligibleScopes <- function(profile,
                           modes = NULL,
                           min_genomes = 40,
                           n_fold = 5,
                           min_test_genomes = 5,
                           min_test_minority = 2) {
  call <- rlang::current_env()

  # Check the arguments.
  .checkArgClass(profile, "profile", "amr_orb_profile", "profileORB()", call = call)
  .checkArgCount(min_genomes, "min_genomes", call = call)
  .checkArgCount(n_fold, "n_fold", min = 2, call = call)
  .checkArgCount(min_test_genomes, "min_test_genomes", call = call)
  .checkArgCount(min_test_minority, "min_test_minority", call = call)
  modes <- modes %||% .MODES$mode_id
  .checkArgModes(modes, call = call)

  settings <- list(
    min_genomes = min_genomes,
    n_fold = n_fold,
    min_test_genomes = min_test_genomes,
    min_test_minority = min_test_minority,
    holdout = .CV_HOLDOUT
  )

  # Decide whether each mode can run on this ORB at all.
  selected <- .MODES[.MODES$mode_id %in% modes, , drop = FALSE]
  feasibility <- .modeFeasibility(selected, profile)

  # Build and judge the scopes of each mode that can run, keeping the genomes of the
  # eligible ones for matrix building.
  runnable <- selected[feasibility$feasible, , drop = FALSE]
  judged <- lapply(seq_len(nrow(runnable)), function(i) {
    mode <- runnable[i, , drop = FALSE]
    built <- .buildScopes(mode, profile)
    decisions <- .judgeScopes(built, mode, profile, settings)
    eligible <- decisions[decisions$eligible, .SCOPE_KEYS]
    list(
      decisions = decisions,
      members = tibble::tibble(
        mode_id = mode$mode_id,
        dplyr::semi_join(built$members, eligible, by = .SCOPE_KEYS)
      )
    )
  })
  scopes <- dplyr::bind_rows(lapply(judged, `[[`, "decisions"))
  members <- dplyr::bind_rows(lapply(judged, `[[`, "members"))

  # If no mode ran, return empty tables that still have every column.
  if (!length(judged)) {
    scopes <- tibble::tibble(
      mode_id = character(), unit = character(), target = character(),
      train_group = character(), test_group = character(),
      n_genomes = integer(), n_resistant = integer(), n_susceptible = integer(),
      test_n_genomes = integer(), test_n_resistant = integer(),
      test_n_susceptible = integer(), eligible = logical(), rule_id = character()
    )
    members <- tibble::tibble(
      mode_id = character(), unit = character(), target = character(),
      train_group = character(), test_group = character(),
      role = character(), genome_id = character(), phenotype = character()
    )
  }

  # Each mode's scope counts, so a mode that ran with nothing eligible is visible; NA if not run.
  per_mode <- match(scopes$mode_id, feasibility$mode_id)
  feasibility$n_scopes <- tabulate(per_mode, nrow(feasibility))
  feasibility$n_eligible <- tabulate(per_mode[scopes$eligible], nrow(feasibility))
  feasibility[!feasibility$feasible, c("n_scopes", "n_eligible")] <- NA_integer_

  structure(
    list(
      profile = profile,
      settings = settings,
      modes = feasibility,
      scopes = scopes,
      members = members
    ),
    class = "amr_eligibility"
  )
}


#' @export
print.amr_eligibility <- function(x, ...) {
  orb <- x$profile$orb
  s <- x$settings

  cat("<amr_eligibility>", orb$dataset_id, "-", orb$dataset_label, "\n")
  cat(
    "  settings :", s$min_genomes, "training genomes,", s$n_fold, "folds; test sets",
    s$min_test_genomes, "genomes,", s$min_test_minority, "of the rarer phenotype\n"
  )

  labels <- format(x$modes$mode_id)
  for (i in seq_len(nrow(x$modes))) {
    mode <- x$modes[i, , drop = FALSE]

    if (!mode$feasible) {
      cat("  ", labels[[i]], " : not run (", mode$rule_id, ")\n", sep = "")
      next
    }

    reasons <- table(x$scopes$rule_id[x$scopes$mode_id == mode$mode_id])
    cat(
      "  ", labels[[i]], " : ", mode$n_eligible, " of ", mode$n_scopes, " eligible",
      if (length(reasons)) paste0("; ", paste(names(reasons), reasons, collapse = ", ")),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}
