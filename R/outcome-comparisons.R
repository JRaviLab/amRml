# Compare cross-target identities
.compareOutcomeEvidence <- function(top, identity) {
  if (!"label" %in% names(top)) top$label <- NA_character_
  outcomes <- dplyr::distinct(top, species, drug_label, drug_or_class, label)
  outcomes$outcome_id <- vapply(seq_len(nrow(outcomes)), function(i) {
    as.character(jsonlite::toJSON(as.list(outcomes[i, c("drug_label", "drug_or_class", "label")]),
      auto_unbox = TRUE, na = "null"))
  }, character(1))
  membership <- top |>
    dplyr::inner_join(outcomes,
      by = c("species", "drug_label", "drug_or_class", "label"),
      relationship = "many-to-one") |>
    dplyr::distinct(dplyr::across(dplyr::all_of(
      c("species", identity, "drug_label", "drug_or_class", "label", "outcome_id"))))
  sharedness <- membership |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("species", identity)))) |>
    dplyr::summarise(
      n_outcomes = dplyr::n_distinct(outcome_id),
      n_model_targets = dplyr::n_distinct(drug_label, drug_or_class),
      n_drug_targets = dplyr::n_distinct(drug_or_class[drug_label == "drug"]),
      n_class_targets = dplyr::n_distinct(drug_or_class[drug_label == "drug_class"]),
      n_mdr_outcomes = dplyr::n_distinct(outcome_id[drug_label == "mdr"]),
      n_mdr_combination_labels = dplyr::n_distinct(
        label[drug_label == "mdr" & drug_or_class == .MDR_COMBINATION]),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      shared_within_mdr = n_mdr_combination_labels >= 2L,
      shared_mdr_single_drug = n_mdr_outcomes > 0L & n_drug_targets > 0L
    )
  sizes <- membership |>
    dplyr::count(species, outcome_id, name = "n_features")
  left <- dplyr::rename(sizes, outcome_a = outcome_id, n_a = n_features)
  right <- dplyr::rename(sizes, outcome_b = outcome_id, n_b = n_features)
  pairs <- dplyr::inner_join(left, right, by = "species", relationship = "many-to-many") |>
    dplyr::filter(outcome_a < outcome_b)
  shared <- dplyr::inner_join(membership, membership,
      by = c("species", identity), relationship = "many-to-many",
      suffix = c("_a", "_b")) |>
    dplyr::filter(outcome_id_a < outcome_id_b) |>
    dplyr::count(species, outcome_a = outcome_id_a, outcome_b = outcome_id_b,
      name = "n_shared")
  overlap <- pairs |>
    dplyr::left_join(shared, by = c("species", "outcome_a", "outcome_b"),
      relationship = "one-to-one") |>
    dplyr::mutate(n_shared = dplyr::coalesce(n_shared, 0L),
      n_union = n_a + n_b - n_shared,
      jaccard = n_shared / n_union,
      fraction_a_shared = n_shared / n_a, fraction_b_shared = n_shared / n_b) |>
    dplyr::left_join(dplyr::rename(outcomes, outcome_a = outcome_id,
        drug_label_a = drug_label, drug_or_class_a = drug_or_class, label_a = label),
      by = c("species", "outcome_a"), relationship = "many-to-one") |>
    dplyr::left_join(dplyr::rename(outcomes, outcome_b = outcome_id,
        drug_label_b = drug_label, drug_or_class_b = drug_or_class, label_b = label),
      by = c("species", "outcome_b"), relationship = "many-to-one")
  list(membership = membership, outcomes = outcomes,
    cross_model_sharedness = sharedness,
    mdr_shared = dplyr::filter(sharedness, shared_within_mdr | shared_mdr_single_drug),
    overlap = overlap)
}
