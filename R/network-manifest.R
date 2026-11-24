# Network evidence remains separate from the completed modeling bundle.
.writeNetworkManifest <- function(out_dir, run_dir, orb_dir, run, core_threshold,
                                   metadata, pairs) {
  tables <- c(
    dyads = "cluster_pair_network_evidence.parquet",
    nodes = "gene_protein_cluster_node_evidence.parquet",
    annotations = "annotation_node_evidence.parquet",
    annotation_edges = "cluster_annotation_neighbor_evidence.parquet",
    feature_seeds = "feature_returned_seed_evidence.parquet",
    cluster_seeds = "gene_protein_seed_evidence.parquet",
    fitted_tasks = "fitted_tasks.parquet",
    model_reliability = "model_reliability.parquet",
    model_reliability_seeds = "model_reliability_by_seed.parquet",
    mcc_audit = "reliability_mcc_audit.parquet",
    protein_crosswalk = "protein_cluster_id_crosswalk.parquet",
    provenance = "network_provenance.parquet"
  )
  run$status <- "success"
  run$finished_at <- as.character(Sys.time())
  run$stages <- list(list(
    name = "export_network", status = "success", started_at = run$started_at,
    finished_at = run$finished_at,
    inputs = .manifestFileInfo(c(
      .runPaths(run_dir)$manifest, .runPaths(run_dir)$matrices,
      .runPaths(run_dir)$tasks, .runPaths(run_dir)$performance,
      .runPaths(run_dir)$importance, .runPaths(run_dir)$predictions,
      list.files(orb_dir, pattern = "(_parquet[.]duckdb$|^manifest_.*[.]json$)",
                 full.names = TRUE)
    ), out_dir),
    outputs = .manifestFileInfo(file.path(out_dir, tables), out_dir),
    metrics = list(n_dyads = nrow(pairs), n_fitted_tasks = nrow(metadata))
  ))
  .manifestWrite(list(
    schema_version = 1L, manifest_type = "amR_ml_network",
    manifest_id = paste0("manifest_", run$run_id),
    manifest_created_at = run$started_at, runs = list(run),
    sources = list(modeling_manifest = .runPaths(run_dir)$manifest,
                   modeling_directory = run_dir, orb_directory = orb_dir),
    artifacts = list(amRviz_network = list(
      status = "ready", producer_run_id = run$run_id,
      directory = ".", path_mode = "relative_to_manifest", tables = as.list(tables),
      parquet_compression = "zstd",
      network = list(
        default_view = list(labels = "real", high_priority = TRUE),
        defaults = list(core_threshold = core_threshold, min_core_seeds = 2L,
                        min_supported_seeds = 2L, min_shuffled_seeds = 2L,
                        shuffled_sd_multiplier = 1, min_real_mcc = 0),
        available = list(targets = I(unique(metadata$target)),
                         units = I(unique(metadata$unit)), seeds = I(unique(metadata$seed)),
                         labels = I(unique(metadata$labels)),
                         encodings = I(unique(metadata$encoding)),
                         feature_types = I(unique(metadata$feature_type)),
                         dyad_types = I(unique(pairs$dyad_type))),
        rules = list(
          core = "Smallest prefix reaching core_threshold after sorting contribution descending, then feature ID; crossing feature included.",
          contribution = "Importance divided by total returned importance within task; means use returned seeds only.",
          robust = "Core in at least min_core_seeds successfully fitted seeds.",
          reliability = "Real MCC strictly above zero, matched shuffled MCC, and shuffled mean plus one sample SD; at least two defined shuffled seeds.",
          dyads = "Every available endpoint must be robust and its representation supported; absent orphan endpoints impose no requirement.",
          priority = "Real-label dyad satisfying robustness and reliability.",
          annotations = "Context only; annotation counts do not increase endpoint scores.",
          undefined_mcc = "Zero only for constant predictions with both truth classes and no missing predictions; raw MCC retained in audit."
        ),
        filtering = list(
          direct = I(c("unit", "target", "encoding", "labels", "dyad_type")),
          recompute = list(
            core_threshold = "Recompute per-task in_core from feature_seeds, then cluster/node/dyad summaries.",
            min_core_seeds = "Recompute robustness from cluster_seeds and fitted_tasks.",
            model_reliability = "Recompute from model_reliability_seeds and mcc_audit, then dyad priority.",
            seed_subset = "Filter fitted_tasks and per-seed evidence, then recompute all summaries."
          )
        ),
        coverage = list(
          all_dyads_retained = TRUE, all_returned_features_retained = TRUE,
          annotation_context = "All protein endpoints, independent of default priority cutoffs.",
          limitations = "Candidate availability and zero coefficients are not inventoried; candidate-aware MAGIC scores cannot be reconstructed.",
          flag_names = "Fields ending in robust75 retain their historical names; the actual core cutoff is defaults.core_threshold."
        )
      )
    ))
  ), file.path(out_dir, "manifest.json"))
}
