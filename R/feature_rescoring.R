#' Score features within each seed
#'
#' @param all_top_features_parquet The path to the Parquet file containing all top features with their importance scores.
#'
#' @returns a tibble of scored top features with their contribution, rank, and rank score within each seed.
#' within a seed
#' contribution is calculated as the importance of a feature divided by the sum of importance scores for all features.
#' Rank is assigned based on the descending order of contribution, and
#' rank score is calculated as (n_features - rank) / (n_features - 1), where n_features is the total number of features within the same seed.
#' rank score is ranged between 0 and 1, with higher values indicating higher importance.
#'
#' @export
#' @examples
#' scoreFeaturesWithinSeed(all_top_features.parquet)
#'
scoreFeaturesWithinSeed <- function(all_top_features_parquet) {

  stopifnot(file.exists(all_top_features_parquet))
  scored_top_features <- arrow::read_parquet(normalizePath(all_top_features_parquet)) |>
    dplyr::filter(!shuffled) |>
    dplyr::select(
      species, drug_label, drug_or_class, seed,
      feature_type, feature_subtype, variable=Variable,
      importance=Importance, sign=Sign
    ) |>
    dplyr::mutate(
      variable = dplyr::case_when(
        feature_type == "domains" ~ sub("_.+$", "", variable),
        feature_type == "proteins" ~ sub("fig.", "fig|", variable, fixed = TRUE),
        feature_type == "args" ~ sub(
          "^X", "",
          gsub("\\.NCBIFAM", "", variable)
        ),
        TRUE ~ variable
      )
    )|>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      feature_type,
      feature_subtype,
      seed
    ) |>
    dplyr::mutate(
      contribution = importance / sum(importance, na.rm = TRUE),
      rank = dplyr::dense_rank(dplyr::desc(contribution)),
      n_features = dplyr::n(),
      rank_score = dplyr::if_else(
        n_features > 1,
        (n_features - rank) / (n_features - 1),
        1
      )
    ) |>
    dplyr::ungroup()

  return(scored_top_features)
}

#' Summarize a variable across seeds
#    "which molecular features are consistent?"
#'
#' @param scored_features The tibble of scored top features with their contribution, rank, and rank score within each seed generated from `scoreFeaturesWithinSeed()`
#'
#' @returns a tibble of summarized features across seeds, including the number of seeds,
#' mean and median rank,
#' mean and median contribution,
#' mean rank score (scaled between 0 and 1 with higher values indicating higher importance),
#' best rank,
#' rank consistency (TRUE if the rank is same across seeds), rank standard deviation (how much the rank varies across seeds),
#' sign consistency (TRUE if the sign is same across seeds), and sign (the sign of the feature; POS, NEG or MIXED).
#'
#' @export
#' @examples
#' summariseFeaturesAcrossSeeds(scoreFeaturesWithinSeed(all_top_features.parquet))
summariseFeaturesAcrossSeeds <- function(scored_features = scoreFeaturesWithinSeed(all_top_features_parquet)) {
 feature_summary <- scored_features |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      feature_type,
      feature_subtype,
      variable
    ) |>
    dplyr::summarise(
    n_seeds = dplyr::n_distinct(seed),
    mean_rank = mean(rank, na.rm = TRUE),
    mean_contribution = mean(contribution, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    median_contribution = median(contribution, na.rm = TRUE),
    mean_rank_score = mean(rank_score, na.rm = TRUE),
    best_rank = min(rank, na.rm = TRUE),
    rank_consistent = dplyr::n_distinct(rank) == 1,
    rank_sd = sd(rank, na.rm = TRUE),
    sign_consistent = dplyr::n_distinct(sign) == 1,
    sign = if (sign_consistent) dplyr::first(sign) else "MIXED",
    .groups = "drop"
  )

  return(feature_summary)
}

#' Build a per-drug top-feature table with cutoffs
#' This returns the top features for each drug/class.
#'
#' @param feature_summary The tibble of summarized features across seeds generated from `summariseFeaturesAcrossSeeds()`
#' @param rank_score_quantile The quantile threshold for filtering features based on their mean rank scores
#' @param threshold_sd_rank The threshold for filtering features based on the standard deviation of ranks across seeds
#'
#' @returns a tibble of top features for each drug/class,
#' filtered based on the specified rank score quantile and standard deviation threshold.
#' The resulting filtered table includes species, drug label, drug or class, feature type, feature subtype, variable,
#' mean rank score, mean rank, best rank, median rank, mean contribution, median contribution, number of seeds,
#' rank standard deviation, sign consistency, and sign.
#' The features are filtered to include only those with a negative sign, low rank standard deviation,
#' and are grouped by drug label and drug or class, retaining only the features with the
#' maximum number of seeds and mean rank scores above the specified quantile threshold.
#' Every drug/class may not have variables from all feature types.
#'
#' @export
topFeaturesPerDrugOrClass <- function(feature_summary = summariseFeaturesAcrossSeeds(scored_features),
                                        rank_score_quantile = 0.95,
                                        threshold_sd_rank = 1
                                        ) {

top_features <- feature_summary |>
  dplyr::filter(
    sign == "NEG",
    rank_sd <= threshold_sd_rank
  ) |>
  dplyr::group_by(species, drug_label, drug_or_class, feature_type, variable) |>
    dplyr::mutate( n_subtype = dplyr::n_distinct(feature_subtype),
  subtype_csv = paste(sort(unique(feature_subtype)), collapse = ",") ) |>
  dplyr::ungroup() |>
  dplyr::group_by(drug_label, drug_or_class) |>
  dplyr::filter(
  #   n_seeds == max(n_seeds),
    mean_rank_score >= quantile(mean_rank_score, rank_score_quantile)
  ) |>
  # dplyr::slice_max(
  #   order_by = mean_rank_score,
  #   n = 10,
  #   with_ties = FALSE
  # ) |>
  dplyr::ungroup() |>
    dplyr::distinct(
      species, drug_label, drug_or_class,
      feature_type, feature_subtype, variable, n_subtype, subtype_csv,
      mean_rank_score, mean_rank, best_rank,
      median_rank, mean_contribution, median_contribution,
      n_seeds, rank_sd,
      sign_consistent, sign
    ) |>
dplyr::filter(sign_consistent)

  return(top_features)
}

#' Aggregate mapped features to protein clusters
#'
#' @param top_features The tibble of top features for each drug/class generated from `topFeaturesPerDrugOrClass()`
#' @param cluster_feature_parquet The path to the Parquet file containing the mapping of features to protein clusters
#'
#' @returns a tibble of summarized clusters, including species, drug label, drug or class, cluster, frequency (number of features in the cluster),
#' number of distinct variables, number of distinct feature types, feature types as a comma-separated string,
#' cluster mean rank score (mean of mean_rank_score), cluster median rank score (median of mean_rank_score),
#' cluster max rank score (max of mean_rank_score), cluster rank score standard deviation, and cluster best rank.
#'
#' @export
summariseClusters <- function(top_features = topFeaturesPerDrugOrClass(feature_summary),
                               cluster_feature_parquet
                              ) {
  stopifnot(file.exists(cluster_feature_parquet))
  cluster_feature <- arrow::read_parquet(normalizePath(cluster_feature_parquet))

  top_clusters <- top_features |>
    dplyr::left_join(cluster_feature, by = dplyr::join_by(variable == feature)) |>
    dplyr::group_by(
      species,
      drug_label,
      drug_or_class,
      cluster
    ) |>
    dplyr::summarise(
      frequency = dplyr::n(),
      n_variables = dplyr::n_distinct(variable),
      variables_csv = paste(sort(unique(variable)), collapse = ","),
      n_feature_types = dplyr::n_distinct(feature_type),
      feature_types_csv = paste(sort(unique(feature_type)), collapse = ","),
      cluster_mean_rank_score = mean(mean_rank_score, na.rm = TRUE),
      cluster_median_rank_score = median(mean_rank_score, na.rm = TRUE),
      cluster_max_rank_score = max(mean_rank_score, na.rm = TRUE),
      cluster_rank_score_sd = sd(mean_rank_score, na.rm = TRUE),
      cluster_best_rank = min(best_rank, na.rm = TRUE),
      .groups = "drop"
    ) |>
  dplyr::arrange(dplyr::desc(frequency), dplyr::desc(n_feature_types))

  return(top_clusters)
}

#' Find clusters that appear across multiple drugs/classes
#'
#' @param top_clusters The tibble of summarized clusters generated from `summariseClusters()`
#' @param label The label to filter clusters by, either "drug" or "drug_class"
#' @param  min_drugs_or_classes The minimum number of distinct drugs or classes required for a cluster to be considered shared
#'
#' @returns a tibble of shared clusters, including cluster, number of distinct drugs or classes, and a comma-separated string of the distinct drugs or classes.
#'
#' @export
findSharedClusters <- function(top_clusters = summariseClusters(top_features, cluster_feature_parquet),
                                label = "drug",
                                 min_drugs_or_classes = 2
                                ) {
  shared_clusters <- top_clusters |>
    dplyr::filter(!is.na(cluster), drug_label == label) |>
    dplyr::group_by(cluster) |>
    dplyr::mutate(
      n_drug_or_class = dplyr::n_distinct(drug_or_class),
      drug_or_class_csv = paste(sort(unique(drug_or_class)), collapse = ", ")
    ) |>
    dplyr::filter(n_drug_or_class >= min_drugs_or_classes) |>
    dplyr::ungroup() |>
    dplyr::select(cluster, n_drug_or_class, drug_or_class_csv) |>
    dplyr::arrange(dplyr::desc(n_drug_or_class))

  return(shared_clusters)
}

#' find the clusters that are unique to a single drug/class
#'
#' @param top_clusters The tibble of summarized clusters generated from `summariseClusters()`
#' @param label The label to filter clusters by, either "drug" or "drug_class"
#' @param protein_names_parquet The path to the Parquet file containing the annotations to protein cluster names
#'
#' @returns
#'
#' @export
#' @examples
findUniqueClusters <- function(top_clusters = summariseClusters(top_features, cluster_feature_parquet),
                            label = "drug",
                            protein_names_parquet
) {

  protein_names <- arrow::read_parquet(normalizePath(protein_names_parquet)) |>
    dplyr::distinct()

  unique_clusters <- top_clusters |>
    dplyr::filter(!is.na(cluster), drug_label == label) |>
    dplyr::group_by(cluster) |>
    dplyr::mutate(
      n_drug_or_class = dplyr::n_distinct(drug_or_class),
      drug_or_class_csv = paste(sort(unique(drug_or_class)), collapse = ", ")
    ) |>
    dplyr::filter(n_drug_or_class == 1) |>
    dplyr::ungroup() |>
    dplyr::select(cluster, drug_or_class_csv, cluster_mean_rank_score) |>
    dplyr::arrange(dplyr::desc(cluster_mean_rank_score)) |>
dplyr::rename(drug_or_class = drug_or_class_csv) |>
    dplyr::left_join(protein_names, by = dplyr::join_by(cluster == proteinID)) |>
    dplyr::rename(cluster_name = proteinName) |>
    dplyr::select(drug_or_class, cluster, cluster_name, cluster_mean_rank_score)

  return(unique_clusters)
}

#' Build a feature network from selected top features and top clusters
#'
#' @param top_features Output of topFeaturesPerDrugOrClass().
#' @param top_clusters Output of summariseClusters().
#' @param cluster_feature_parquet Path to the Parquet file mapping variables to clusters.
#' @param protein_names_parquet Path to the Parquet file with cluster name annotations.
#'
#' @return A list with feature_table, cluster_table, nodes, edges, and graph.
#' @export
buildFeatureNetwork <- function(top_features,
                               top_clusters,
                               cluster_feature_parquet,
                               protein_names_parquet
                              ) {
  stopifnot(is.data.frame(top_features))
  stopifnot(is.data.frame(top_clusters))
  stopifnot(file.exists(cluster_feature_parquet))
  stopifnot(file.exists(protein_names_parquet))

  required_feature_cols <- c(
    "species", "drug_label", "drug_or_class",
    "feature_type", "variable",
    "mean_rank_score"
  )
  required_cluster_cols <- c(
    "species", "drug_label", "drug_or_class",
    "cluster", "cluster_mean_rank_score"
  )

  missing_feature_cols <- setdiff(required_feature_cols, names(top_features))
  missing_cluster_cols <- setdiff(required_cluster_cols, names(top_clusters))

  if (length(missing_feature_cols) > 0) {
    stop("top_features is missing required columns: ",
         paste(missing_feature_cols, collapse = ", "))
  }
  if (length(missing_cluster_cols) > 0) {
    stop("top_clusters is missing required columns: ",
         paste(missing_cluster_cols, collapse = ", "))
  }

  cluster_feature <- arrow::read_parquet(normalizePath(cluster_feature_parquet)) |>
    dplyr::distinct()

  protein_names <- arrow::read_parquet(normalizePath(protein_names_parquet)) |>
    dplyr::distinct()

  make_model_id <- function(drug_label, drug_or_class) {
    paste(drug_label, drug_or_class, sep = ".")
  }

  feature_table <- top_features |>
    dplyr::mutate(model_id = make_model_id(drug_label, drug_or_class)) |>
    dplyr::group_by(species, model_id, feature_type, variable) |>
    dplyr::summarise(
      mean_score = mean(mean_rank_score, na.rm = TRUE),
      .groups = "drop"
    )

  cluster_table <- top_clusters |>
    dplyr::mutate(model_id = make_model_id(drug_label, drug_or_class)) |>
    dplyr::left_join(
      protein_names,
      by = dplyr::join_by(cluster == proteinID)
    ) |>
    dplyr::group_by(species, model_id, cluster, proteinName) |>
    dplyr::summarise(
      mean_score = mean(cluster_mean_rank_score, na.rm = TRUE),
      .groups = "drop"
    )

  model_nodes <- dplyr::bind_rows(
    feature_table |>
      dplyr::distinct(species, model_id),
    cluster_table |>
      dplyr::distinct(species, model_id)
  ) |>
    dplyr::distinct(species, model_id) |>
    dplyr::transmute(
      name = model_id,
      label = model_id,
      node_type = "model",
      species = species,
      score = NA_real_,
      breadth = NA_real_,
      node_size = 4
    )

  feature_nodes <- feature_table |>
    dplyr::group_by(species, variable) |>
    dplyr::summarise(
      score = mean(mean_score, na.rm = TRUE),
      breadth = dplyr::n_distinct(model_id),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      name = variable,
      label = variable,
      node_type = "feature",
      species = species,
      score = score,
      breadth = breadth,
      node_size = pmax(3, pmin(10, breadth + 2))
    )

  cluster_nodes <- cluster_table |>
    dplyr::group_by(species, cluster, proteinName) |>
    dplyr::summarise(
      score = mean(mean_score, na.rm = TRUE),
      breadth = dplyr::n_distinct(model_id),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      name = cluster,
      label = dplyr::if_else(
        is.na(proteinName) | proteinName == "",
        cluster,
        proteinName
      ),
      node_type = "cluster",
      species = species,
      score = score,
      breadth = breadth,
      node_size = pmax(3, pmin(10, breadth + 2))
    )

  nodes <- dplyr::bind_rows(model_nodes, feature_nodes, cluster_nodes) |>
    dplyr::distinct(name, .keep_all = TRUE)

  feature_edges <- feature_table |>
    dplyr::transmute(
      from = model_id,
      to = variable,
      weight = mean_score,
      edge_type = "model_feature"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  cluster_edges <- cluster_table |>
    dplyr::transmute(
      from = model_id,
      to = cluster,
      weight = mean_score,
      edge_type = "model_cluster"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  feature_cluster_edges <- feature_table |>
    dplyr::left_join(cluster_feature, by = dplyr::join_by(variable == feature)) |>
    dplyr::filter(!is.na(cluster)) |>
    dplyr::transmute(
      from = variable,
      to = cluster,
      weight = 1,
      edge_type = "feature_cluster"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  edges <- dplyr::bind_rows(feature_edges, feature_cluster_edges)

    edges <- dplyr::bind_rows(edges, cluster_edges)


  missing_vertices <- setdiff(unique(c(edges$from, edges$to)), nodes$name)
  if (length(missing_vertices) > 0) {
    extra_nodes <- tibble::tibble(name = missing_vertices) |>
      dplyr::mutate(
        label = name,
        node_type = dplyr::case_when(
          grepl("^drug\\.|^drug_class\\.", name) ~ "model",
          grepl("^fig\\||^cluster", name) ~ "cluster",
          TRUE ~ "feature"
        ),
        species = NA_character_,
        score = NA_real_,
        breadth = NA_real_,
        node_size = 4
      )

    nodes <- dplyr::bind_rows(nodes, extra_nodes) |>
      dplyr::distinct(name, .keep_all = TRUE)
  }

  graph <- if (nrow(edges) > 0) {
    igraph::graph_from_data_frame(
      d = edges,
      directed = FALSE,
      vertices = nodes
    )
  } else {
    igraph::graph_from_data_frame(
      d = data.frame(from = character(), to = character()),
      directed = FALSE,
      vertices = nodes
    )
  }

  feature_network <- list(
    feature_table = feature_table,
    cluster_table = cluster_table,
    nodes = nodes,
    edges = edges,
    graph = graph
  )

  return(feature_network)
}
#' Plot the feature network with networkD3
#'
#' @param feature_network Output of \code{identifyTopFeatures()}.
#' @param height Widget height in pixels.
#' @param width Widget width.
#'
#' @returns A \code{networkD3} widget.
#' @export
plotFeatureNetworkD3 <- function(feature_network,
                                height = 800,
                                width = "100%"
                              ) {

  stopifnot(is.list(feature_network))
  stopifnot(!is.null(feature_network$nodes))
  stopifnot(!is.null(feature_network$edges))

  nodes <- feature_network$nodes |>
    dplyr::distinct(name, .keep_all = TRUE) |>
    dplyr::mutate(
      id = dplyr::row_number() - 1L,
      group = node_type,
      title = paste0(
        "<b>", label, "</b>",
        ifelse(is.na(species), "", paste0("<br>Species: ", species)),
        ifelse(is.na(score), "", paste0("<br>Score: ", signif(score, 3))),
        ifelse(is.na(breadth), "", paste0("<br>Breadth: ", breadth))
      )
    )

  links <- feature_network$edges |>
    dplyr::filter(!is.na(from), !is.na(to)) |>
    # dplyr::filter(
    #   show_direct_model_cluster | edge_type != "model_cluster"
    # ) |>
    dplyr::left_join(
      nodes |> dplyr::select(name, id),
      by = c("from" = "name")
    ) |>
    dplyr::rename(source = id) |>
    dplyr::left_join(
      nodes |> dplyr::select(name, id),
      by = c("to" = "name")
    ) |>
    dplyr::rename(target = id) |>
    dplyr::filter(!is.na(source), !is.na(target)) |>
    dplyr::mutate(
      value = dplyr::if_else(is.na(weight), 1, weight)
    ) |>
    dplyr::select(source, target, value, edge_type)

  stopifnot(nrow(nodes) > 0)
  stopifnot(nrow(links) > 0)

  colour_scale <- networkD3::JS(
    "d3.scaleOrdinal()
      .domain(['model', 'feature', 'cluster'])
      .range(['#4C78A8', '#F58518', '#54A24B'])"
  )

  networkD3::forceNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "label",
    Group = "group",
    opacity = 0.9,
    zoom = TRUE,
    fontSize = 14,
    # nodeWidth = 24,
    height = height,
    width = width,
    colourScale = colour_scale,
    linkDistance = networkD3::JS(
      "function(d) {
         if (d.edge_type === 'feature_cluster') return 60;
         if (d.edge_type === 'model_feature') return 120;
         return 90;
       }"
    )
  )
}


# final run would be:
# top_features <- topFeaturesPerDrugOrClass(rank_score_quantile = 0.75)
# top_clusters <- summariseClusters(top_features, cluster_feature_parquet = cluster_feature_parquet)
# feature_network <- buildFeatureNetwork(top_features = top_features, top_clusters = top_clusters,
#   cluster_feature_parquet = cluster_feature_parquet, protein_names_parquet = protein_names_parquet)
# plotFeatureNetworkD3(feature_network)
