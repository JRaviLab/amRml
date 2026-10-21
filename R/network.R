# Build network evidence from completed amRml modeling tasks and the existing ORB.
# These helpers summarize fitted models and connect their features to ORB clusters
# and annotations, without refitting models or modifying the source data.

# Calculate each returned feature's share of importance within its fitted task.
# Sort features by contribution (and feature ID to break ties), then mark the
# smallest set whose combined importance reaches the requested core threshold.
# Keep all returned features so later steps can compare core and non-core evidence.
.scoreNetworkTasks <- function(importance, metadata, core_threshold) {
  stopifnot(core_threshold > 0, core_threshold <= 1)
  r <- importance |>
    dplyr::inner_join(metadata, by = "task_id", relationship = "many-to-one")
  if (anyNA(r[c("Variable", "Importance", "Sign")]) ||
      any(!is.finite(r$Importance) | r$Importance <= 0) ||
      anyDuplicated(r[c("task_id", "Variable")])) {
    stop("Invalid returned feature importance.", call. = FALSE)
  }
  r |>
    dplyr::group_by(task_id) |>
    dplyr::mutate(contribution = Importance / sum(Importance)) |>
    dplyr::arrange(dplyr::desc(contribution), Variable, .by_group = TRUE) |>
    dplyr::mutate(in_core = dplyr::row_number() <=
      which(cumsum(contribution) >= core_threshold - 1e-12)[[1L]]) |>
    dplyr::ungroup()
}

# Prepare gene and protein cluster evidence from the scored model features.
# Both model types use cluster-level predictors, so each feature maps directly
# to its corresponding ORB cluster. Restore the "fig|" prefix in CD-HIT protein
# representative IDs when glmnet has written it as "fig.", and record whether
# each cluster appears in the fixed ORB dyad graph. Do not combine coefficients.
.prepareNetworkNodes <- function(scored, pairs) {
  rows <- scored |>
    dplyr::filter(feature_type %in% c("gene", "protein")) |>
    dplyr::transmute(task_id, unit, target, labels, seed, encoding, feature_type,
                     node_id = dplyr::if_else(feature_type == "protein",
                       sub("^fig\\.", "fig|", Variable), Variable),
                     contribution, in_core, sign = Sign,
                     n_returned_features = 1L) |>
    dplyr::mutate(in_orb_graph = dplyr::if_else(feature_type == "protein",
      node_id %in% pairs$protein_cluster, node_id %in% pairs$gene_cluster))
  if (anyDuplicated(rows[c("task_id", "feature_type", "node_id")])) {
    stop("Multiple model features map to the same cluster ID within a task.", call. = FALSE)
  }
  rows
}

# Summarize how often each cluster or annotation was returned across fitted seeds.
# Record its total and average returned importance and how often it entered the
# core group. Use all successful fitted tasks to count seeds, including those
# where a feature was not returned. These tables cannot distinguish an eligible
# zero coefficient from a feature that was unavailable to the model; the final
# candidate-aware MAGIC score requires a separate candidate inventory.
.summarizeNetworkNodes <- function(rows, fitted) {
  groups <- c("unit", "target", "labels", "encoding", "feature_type")
  fits <- fitted |>
    dplyr::count(dplyr::across(dplyr::all_of(groups)), name = "n_fitted_seeds")
  rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(groups, "node_id")))) |>
    dplyr::summarise(
      n_returned_seeds = dplyr::n_distinct(seed),
      n_core_seeds = sum(in_core),
      total_returned_contribution = sum(contribution),
      mean_returned_contribution = mean(contribution),
      n_returned_features = sum(n_returned_features),
      .groups = "drop"
    ) |>
    dplyr::left_join(fits, by = groups, relationship = "many-to-one") |>
    dplyr::mutate(
      # Require core importance in at least two fitted seeds. This measures
      # recurrence among returned coefficients, not candidate-aware importance.
      reproducibly_important = n_core_seeds >= 2L & n_fitted_seeds >= 2L,
      observed_in_fitted_fraction = n_returned_seeds / n_fitted_seeds
    )
}

# Fetch only the requested cluster or feature mappings from one ORB view.
# Query keys in batches so large annotation tables stay in DuckDB instead of
# being loaded into R. Avoid the virtual dyad_feature relation, which expands
# into a billion features that the cluster network does not need.
.queryNetworkKeys <- function(con, view, column, keys, other) {
  if (!length(keys)) return(data.frame())
  pieces <- split(unique(keys), ceiling(seq_along(unique(keys)) / 500L))
  dplyr::bind_rows(lapply(pieces, function(part) {
    quoted <- as.character(DBI::dbQuoteString(con, part))
    DBI::dbGetQuery(con, paste0(
      "SELECT DISTINCT CAST(", other, " AS VARCHAR) AS value, ",
      "CAST(", column, " AS VARCHAR) AS lookup_key FROM ", view,
      " WHERE ", column, " IN (", paste(quoted, collapse = ", "), ")"
    ))
  }))
}

#' Export MAGIC gene-protein cluster networks for amRviz
#'
#' Build a fixed ORB cluster-pair graph with real and shuffled feature evidence,
#' model reliability, and separate connected functional annotations. Does not
#' refit models. Gene and protein models already use cluster predictors;
#' protein cluster IDs are named by CD-HIT representatives. Functional
#' annotations do not increase either endpoint's feature score.
#'
#' @param run_dir Completed modeling run with tasks.parquet, matrices.parquet,
#'   results/performance.parquet, results/importance.parquet, and predictions.
#' @param orb_dir Parquet-backed ORB with cluster dyads, protein_members, and annotation
#'   views. Opened read-only.
#' @param out_dir Output directory outside the modeling run and ORB.
#' @param core_threshold Share of returned coefficient importance for the core.
#' @return Invisibly, a list containing the per-target summary and output directory.
#' @details Undefined MCC is provisionally replaced by zero only when the saved
#'   predictions show constant class predictions and both true classes in the
#'   test set. The raw modeling performance Parquet remains unchanged. Model
#'   reliability is assessed by `assessModelReliability()`. Candidate-aware direct
#'   feature scores are not approximated from the matrix-wide feature count;
#'   node tables here report observed returned-coefficient evidence instead.
#'   Writes Parquet evidence tables only. Does not create review or Cytoscape CSVs.
#' @export
exportNetwork <- function(run_dir, orb_dir, out_dir,
                         core_threshold = 0.75) {
  if (length(core_threshold) != 1L || !is.finite(core_threshold) ||
      core_threshold <= 0 || core_threshold > 1) {
    stop("core_threshold must be in (0, 1].", call. = FALSE)
  }
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  orb_dir <- normalizePath(orb_dir, winslash = "/", mustWork = TRUE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  if (out_dir %in% c(run_dir, orb_dir) ||
      startsWith(out_dir, paste0(run_dir, "/")) ||
      startsWith(out_dir, paste0(orb_dir, "/"))) {
    stop("Output must be outside the modeling run and ORB.", call. = FALSE)
  }
  inputs <- c("tasks.parquet", "matrices.parquet", "results/performance.parquet",
              "results/importance.parquet", "results/predictions.parquet")
  if (!all(file.exists(file.path(run_dir, inputs)))) {
    stop("Missing run output: ", paste(inputs[!file.exists(file.path(run_dir, inputs))],
                                      collapse = ", "), call. = FALSE)
  }
  tasks <- arrow::read_parquet(file.path(run_dir, inputs[[1L]]))
  matrices <- arrow::read_parquet(file.path(run_dir, inputs[[2L]]))
  performance <- arrow::read_parquet(file.path(run_dir, inputs[[3L]]))
  importance <- arrow::read_parquet(file.path(run_dir, inputs[[4L]]))
  predictions <- arrow::read_parquet(file.path(run_dir, inputs[[5L]]))
  required_tasks <- c("task_id", "matrix_id", "labels", "seed", "status")
  required_matrices <- c("matrix_id", "mode_id", "unit", "target", "feature_type", "encoding")
  if (!all(required_tasks %in% names(tasks)) ||
      !all(required_matrices %in% names(matrices)) ||
      !all(c("task_id", "mcc") %in% names(performance)) ||
      !all(c("task_id", "Variable", "Importance", "Sign") %in% names(importance)) ||
      !all(c("task_id", "truth", "predicted") %in% names(predictions)) ||
      anyDuplicated(tasks$task_id) || anyDuplicated(matrices$matrix_id) ||
      anyDuplicated(performance$task_id)) {
    stop("Missing columns or duplicate identities in modeling output.", call. = FALSE)
  }
  metadata <- tasks |>
    dplyr::filter(status == "fitted", labels %in% c("real", "shuffled")) |>
    dplyr::left_join(matrices, by = "matrix_id", relationship = "many-to-one") |>
    dplyr::filter(mode_id == "vanilla") |>
    dplyr::select(task_id, matrix_id, unit, target, labels, seed, feature_type, encoding)
  if (!nrow(metadata) || anyNA(metadata) ||
      anyDuplicated(metadata[c("unit", "target", "labels", "seed", "feature_type", "encoding")])) {
    stop("Expected distinct fitted vanilla tasks with complete scope metadata.", call. = FALSE)
  }
  if ("seed" %in% names(performance) &&
      any(metadata$seed != performance$seed[match(metadata$task_id, performance$task_id)],
          na.rm = TRUE)) {
    stop("Inconsistent fitted-task seeds.", call. = FALSE)
  }
  performance <- metadata |>
    dplyr::inner_join(dplyr::select(performance, -dplyr::any_of("seed")),
                      by = "task_id", relationship = "one-to-one")
  if (nrow(performance) != nrow(metadata)) {
    stop("A fitted task has no performance record.", call. = FALSE)
  }

  scored <- .scoreNetworkTasks(importance, metadata, core_threshold)
  if (dplyr::n_distinct(scored$task_id) != nrow(metadata)) {
    stop("A fitted task has no usable returned coefficient importance.", call. = FALSE)
  }
  if (any(abs((scored |>
    dplyr::group_by(task_id) |>
    dplyr::summarise(total = sum(contribution), .groups = "drop"))$total - 1) > 1e-7)) {
    stop("Normalized contribution is not conserved within fitted tasks.", call. = FALSE)
  }
  # Keep one record for every successful task when counting fitted seeds.
  # A model can fit successfully while returning only a few nonzero features,
  # so returned coefficient rows alone cannot tell us how many seeds were fitted.
  fitted <- metadata |>
    dplyr::select(unit, target, labels, seed, encoding, feature_type)

  audit <- predictions |>
    dplyr::group_by(task_id) |>
    dplyr::summarise(
      n_missing = sum(is.na(truth) | is.na(predicted)),
      n_true = dplyr::n_distinct(truth, na.rm = TRUE),
      n_predicted = dplyr::n_distinct(predicted, na.rm = TRUE),
      .groups = "drop"
    )
  performance <- performance |>
    dplyr::left_join(audit, by = "task_id", relationship = "one-to-one") |>
    dplyr::mutate(
      mcc_adjusted = is.na(mcc) & !is.na(n_missing) & n_missing == 0L &
        n_true == 2L & n_predicted == 1L,
      mcc_for_reliability = dplyr::if_else(mcc_adjusted, 0, mcc)
    )
  species <- basename(dirname(orb_dir))
  perf_reliability <- performance |>
    dplyr::transmute(species = species, drug_label = unit,
                     drug_or_class = target, feature_type,
                     feature_subtype = encoding, seed,
                     shuffled = labels == "shuffled", mcc = mcc_for_reliability)
  reliability <- local({
    path <- tempfile(fileext = ".parquet")
    on.exit(unlink(path), add = TRUE)
    arrow::write_parquet(perf_reliability, path)
    assessModelReliability(path)
  })
  representation <- reliability$by_representation |>
    dplyr::rename(unit = drug_label, target = drug_or_class,
                  encoding = feature_subtype) |>
    dplyr::select(unit, target, feature_type, encoding,
                  evidence_status, supported, n_supported_seeds,
                  n_paired_assessable, n_shuffled_defined)

  # Read the fixed gene-cluster/protein-cluster pairs from the ORB rather than
  # rebuilding edges from model results. The ORB views use relative Parquet paths,
  # so open the database from its directory and restore the working directory
  # and database connection when the query finishes.
  orb_results <- local({
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(orb_dir)
    databases <- list.files(".", pattern = "_parquet[.]duckdb$")
    if (length(databases) != 1L) {
      stop("Expected one *_parquet.duckdb ORB database.", call. = FALSE)
    }
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = databases[[1L]], read_only = TRUE)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    needed <- c("dyads", "protein_members")
    if (!all(needed %in% DBI::dbListTables(con))) {
      stop("ORB lacks the cluster dyads or protein_members relation.", call. = FALSE)
    }
    pairs <- DBI::dbGetQuery(con,
      "SELECT DISTINCT CAST(gene AS VARCHAR) AS gene_cluster, CAST(protein AS VARCHAR) AS protein_cluster FROM dyads")
    if (!nrow(pairs) || anyNA(pairs) ||
        anyDuplicated(pairs[c("gene_cluster", "protein_cluster")])) {
      stop("Invalid ORB cluster-pair graph.", call. = FALSE)
    }
    # Check that each dyad's protein ID is a known CD-HIT cluster representative.
    # Representatives can have "fig|" IDs just like protein members, so the ID
    # prefix alone cannot tell us whether the ORB contains the correct dyads.
    unknown <- DBI::dbGetQuery(con,
      "SELECT COUNT(*) AS n FROM (
         SELECT DISTINCT CAST(protein AS VARCHAR) AS cluster FROM dyads
         EXCEPT
         SELECT DISTINCT CAST(cluster AS VARCHAR) AS cluster FROM protein_members
       )")$n[[1L]]
    if (unknown > 0L) {
      stop("ORB dyads contain protein IDs absent from CD-HIT cluster representatives (possibly individual protein IDs). Rebuild the ORB cluster dyads.",
           call. = FALSE)
    }
    list(pairs = pairs)
  })
  pairs <- orb_results$pairs
  cluster_rows <- .prepareNetworkNodes(scored, pairs)
  nodes <- .summarizeNetworkNodes(cluster_rows, fitted) |>
    dplyr::mutate(evidence_origin = dplyr::if_else(
      feature_type == "protein", "direct_protein_cluster_model", "direct_gene_model"))
  annotations <- scored |>
    dplyr::filter(feature_type %in% c("Pfam", "COG", "AMRFinder", "DefenseCas")) |>
    dplyr::transmute(unit, target, labels, seed, encoding, feature_type,
                     node_id = Variable, contribution, in_core,
                     n_returned_features = 1L) |>
    .summarizeNetworkNodes(fitted)
  annotation_evidence <- annotations |>
    dplyr::left_join(representation, by = c("unit", "target", "feature_type", "encoding"),
                     relationship = "many-to-one")

  scopes <- metadata |>
    dplyr::distinct(unit, target, labels, encoding)
  graph <- dplyr::cross_join(scopes, pairs)
  graph$gene_degree <- as.integer(table(pairs$gene_cluster)[graph$gene_cluster])
  graph$protein_degree <- as.integer(table(pairs$protein_cluster)[graph$protein_cluster])
  for (typ in c("gene", "protein")) {
    n <- nodes |>
      dplyr::filter(feature_type == typ) |>
      dplyr::select(unit, target, labels, encoding, node_id,
                    reproducibly_important, n_core_seeds, n_fitted_seeds,
                    mean_returned_contribution)
    prefix <- if (typ == "gene") "gene" else "protein"
    names(n)[names(n) == "node_id"] <- paste0(prefix, "_cluster")
    names(n)[names(n) == "reproducibly_important"] <- paste0(prefix, "_robust75")
    names(n)[names(n) == "n_core_seeds"] <- paste0(prefix, "_core_seeds")
    names(n)[names(n) == "n_fitted_seeds"] <- paste0(prefix, "_fitted_seeds")
    names(n)[names(n) == "mean_returned_contribution"] <- paste0(prefix, "_mean_returned_contribution")
    graph <- dplyr::left_join(graph, n,
      by = c("unit", "target", "labels", "encoding", paste0(prefix, "_cluster")),
      relationship = "many-to-one")
    col <- paste0(prefix, "_robust75")
    graph[[col]] <- dplyr::coalesce(graph[[col]], FALSE)
  }
  for (typ in c("gene", "protein")) {
    rep <- representation |>
      dplyr::filter(feature_type == typ) |>
      dplyr::select(unit, target, encoding, supported, evidence_status)
    names(rep)[names(rep) == "supported"] <- paste0(typ, "_representation_supported")
    names(rep)[names(rep) == "evidence_status"] <- paste0(typ, "_representation_status")
    graph <- dplyr::left_join(graph, rep, by = c("unit", "target", "encoding"),
                              relationship = "many-to-one")
  }
  graph <- graph |>
    dplyr::mutate(
      gene_representation_supported = dplyr::coalesce(gene_representation_supported, FALSE),
      protein_representation_supported = dplyr::coalesce(protein_representation_supported, FALSE),
      pair_robust75 = gene_robust75 & protein_robust75,
      pair_reliability_supported = gene_representation_supported & protein_representation_supported,
      priority_under_reliability = pair_robust75 & pair_reliability_supported,
      high_priority = labels == "real" & priority_under_reliability,
      protein_evidence_origin = "direct_protein_cluster_model"
    )
  background <- graph |>
    dplyr::filter(labels == "shuffled") |>
    dplyr::select(unit, target, encoding, gene_cluster, protein_cluster,
                  shuffled_pair_robust75 = pair_robust75,
                  shuffled_pair_reliability_supported = priority_under_reliability)
  real <- graph |>
    dplyr::filter(labels == "real") |>
    dplyr::left_join(background,
      by = c("unit", "target", "encoding", "gene_cluster", "protein_cluster"),
      relationship = "one-to-one") |>
    dplyr::mutate(
      shuffled_pair_robust75 = dplyr::coalesce(shuffled_pair_robust75, FALSE),
      shuffled_pair_reliability_supported = dplyr::coalesce(shuffled_pair_reliability_supported, FALSE)
    )
  # Find the protein clusters adjacent to reproducibly important gene/protein
  # pairs. Include robust real-label pairs even if their models did not pass the
  # reliability check, so annotation context remains available for review.
  robust_pairs <- real |>
    dplyr::filter(pair_robust75)

  # Retrieve ORB annotations for protein clusters connected to robust pairs.
  # Keep Pfam, COG, AMRFinder, and DefenseCas links separate from the endpoint
  # scores: a cluster with many annotations does not gain extra importance.
  neighbor_clusters <- unique(robust_pairs$protein_cluster)
  neighbor_links <- local({
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(orb_dir)
    db <- list.files(".", pattern = "_parquet[.]duckdb$")
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db[[1L]], read_only = TRUE)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    views <- c(Pfam = "dyad_ann_pfam", COG = "dyad_ann_cog",
               AMRFinder = "dyad_ann_amr", DefenseCas = "dyad_ann_defense")
    available <- DBI::dbListTables(con)
    dplyr::bind_rows(lapply(names(views), function(typ) {
      if (!views[[typ]] %in% available) {
        warning("Missing ORB annotation view: ", views[[typ]], call. = FALSE)
        return(data.frame())
      }
      found <- .queryNetworkKeys(con, views[[typ]], "cluster",
                                      neighbor_clusters, "feature")
      if (!nrow(found)) return(data.frame())
      data.frame(feature_type = typ, protein_cluster = found$lookup_key,
                 annotation_key = found$value)
    }))
  })
  if (!nrow(neighbor_links)) {
    neighbor_links <- data.frame(feature_type = character(), protein_cluster = character(),
                                 annotation_key = character())
  }
  # Match ORB annotation IDs to model feature IDs without dropping annotations
  # that were not returned by a fitted model. Prefer an exact feature-ID match;
  # use the older AMRFinder spelling only when an exact match is unavailable.
  # Unmatched ORB annotations still provide context, but have no model score.
  annotation_ids <- annotation_evidence |>
    dplyr::distinct(feature_type, node_id) |>
    dplyr::mutate(legacy_key = dplyr::if_else(feature_type == "AMRFinder",
      sub("^X", "", gsub("\\.NCBIFAM", "", node_id)), node_id))
  raw_match <- annotation_ids |>
    dplyr::transmute(feature_type, annotation_key = node_id, node_id) |>
    dplyr::semi_join(neighbor_links, by = c("feature_type", "annotation_key"))
  legacy_match <- annotation_ids |>
    dplyr::filter(node_id != legacy_key) |>
    dplyr::anti_join(raw_match, by = c("feature_type", "node_id")) |>
    dplyr::transmute(feature_type, annotation_key = legacy_key, node_id) |>
    dplyr::semi_join(neighbor_links, by = c("feature_type", "annotation_key"))
  mapping <- dplyr::bind_rows(raw_match, legacy_match) |>
    dplyr::distinct()
  if (anyDuplicated(mapping[c("feature_type", "annotation_key")])) {
    stop("An ORB annotation key matches multiple model feature identities.", call. = FALSE)
  }
  matched <- neighbor_links |>
    dplyr::left_join(mapping, by = c("feature_type", "annotation_key"),
                     relationship = "many-to-one") |>
    dplyr::mutate(node_id = dplyr::coalesce(node_id, annotation_key))
  neighbors <- robust_pairs |>
    dplyr::distinct(unit, target, encoding, protein_cluster) |>
    dplyr::left_join(matched, by = "protein_cluster", relationship = "many-to-many") |>
    dplyr::filter(!is.na(feature_type)) |>
    dplyr::cross_join(tibble::tibble(labels = c("real", "shuffled"))) |>
    dplyr::left_join(annotation_evidence,
      by = c("unit", "target", "encoding", "labels", "feature_type", "node_id"),
      relationship = "many-to-one") |>
    dplyr::mutate(returned_nonzero = !is.na(n_returned_seeds),
                  robust75 = dplyr::coalesce(reproducibly_important, FALSE),
                  annotation_evidence_status = dplyr::if_else(returned_nonzero,
                    "returned", "not_returned_candidate_unknown"))

  by_target <- real |>
    dplyr::group_by(unit, target) |>
    dplyr::summarise(
      n_fixed_pair_encoding_records = dplyr::n(),
      n_real_robust_pair_encoding_records = sum(pair_robust75),
      n_high_priority_pair_encoding_records = sum(high_priority),
      n_robust_pair_encoding_records_also_shuffled =
        sum(pair_robust75 & shuffled_pair_robust75),
      n_distinct_real_robust_pairs = dplyr::n_distinct(
        gene_cluster[pair_robust75], protein_cluster[pair_robust75]),
      n_distinct_high_priority_pairs = dplyr::n_distinct(
        gene_cluster[high_priority], protein_cluster[high_priority]),
      .groups = "drop"
    )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  provenance <- data.frame(
    modeling_run = run_dir,
    orb_directory = orb_dir,
    n_cluster_pairs = nrow(pairs),
    n_fitted_tasks = nrow(metadata),
    core_threshold = core_threshold,
    min_supported_seeds = 2L,
    shuffled_reference = "mean_plus_one_sd",
    constant_prediction_mcc = "zero_when_both_truth_classes_present"
  )
  arrow::write_parquet(provenance, file.path(out_dir, "network_provenance.parquet"))
  arrow::write_parquet(scored, file.path(out_dir, "feature_returned_seed_evidence.parquet"))
  arrow::write_parquet(nodes, file.path(out_dir, "gene_protein_cluster_node_evidence.parquet"))
  arrow::write_parquet(annotation_evidence, file.path(out_dir, "annotation_node_evidence.parquet"))
  arrow::write_parquet(graph, file.path(out_dir, "cluster_pair_network_evidence.parquet"))
  arrow::write_parquet(neighbors, file.path(out_dir, "cluster_annotation_neighbor_evidence.parquet"))
  arrow::write_parquet(reliability$by_representation, file.path(out_dir, "model_reliability.parquet"))
  arrow::write_parquet(performance |>
    dplyr::select(task_id, unit, target, labels, seed, feature_type, encoding,
                  mcc, mcc_adjusted, mcc_for_reliability),
    file.path(out_dir, "reliability_mcc_audit.parquet"))
  protein_ids <- scored |>
    dplyr::filter(feature_type == "protein") |>
    dplyr::distinct(Variable) |>
    dplyr::transmute(model_feature = Variable,
                     protein_cluster = sub("^fig\\.", "fig|", Variable),
                     in_orb_graph = protein_cluster %in% pairs$protein_cluster)
  arrow::write_parquet(protein_ids,
    file.path(out_dir, "protein_cluster_id_crosswalk.parquet"))
  if (!any(cluster_rows$feature_type == "gene" & cluster_rows$in_orb_graph) ||
      !any(cluster_rows$feature_type == "protein" & cluster_rows$in_orb_graph)) {
    warning("No matching gene or protein clusters in ORB dyads; check feature ID normalization and ORB version.",
            call. = FALSE)
  }
  arrow::write_parquet(cluster_rows,
    file.path(out_dir, "gene_protein_seed_evidence.parquet"))

  message("MAGIC network export: ", nrow(pairs), " fixed cluster pairs; ",
          sum(real$high_priority), " high-priority pair/encoding records (",
          sum(by_target$n_distinct_high_priority_pairs), " distinct pairs).\n",
          "Outputs: ", out_dir)
  invisible(list(summary = by_target, out_dir = out_dir))
}
