
#' Filter the optimal model based on MCC threshold and comparison to shuffled data
#'
#' @param all_performance_parquet The path to the all performance parquet file 
#' @param MCC_threshold Numeric The minimum MCC to filter the optimal model (default NULL, no filtering)
#' @param compare_to_shuffled Logical whether to compare the MCC of the model with and without shuffled data (default is TRUE) 
#'
#' @returns a tibble of optimal models
#'
#' @export
#' @examples
filterOptimalModel <- function(all_performance_parquet, 
  MCC_threshold = NULL, 
  compare_to_shuffled = TRUE) 
  {
  stopifnot(file.exists(all_performance_parquet))
  all_perf <- arrow::read_parquet(normalizePath(all_performance_parquet)) 
  
all_perf |>
  dplyr::select(
    species, drug_label, drug_or_class,
    seed, feature_type, feature_subtype,
    model, fit_penalty, fit_mixture,
    shuffled, mcc
  ) |>
  tidyr::pivot_wider(
    names_from = shuffled,
    values_from = mcc,
    names_prefix = "shuffled_"
  ) |>
  dplyr::rename(
    nonshuffled_MCC = shuffled_FALSE,
    shuffled_MCC = shuffled_TRUE
  ) |>
  dplyr::mutate(
    MCC_diff = nonshuffled_MCC - shuffled_MCC
  ) |>
  dplyr::filter(
  if (!is.null(MCC_threshold))(nonshuffled_MCC >= MCC_threshold) else TRUE,
  if (compare_to_shuffled) (MCC_diff > 0 | is.na(shuffled_MCC)) else TRUE
)
}

#' Score features within each seed
#'
#' @param all_top_features_parquet The path to the Parquet file containing all top features with their importance scores.
#' @param core_contribution_threshold The cumulative-contribution cutoff (default 0.75, i.e. 75%) used to flag whether a feature falls within the "core" set of features that jointly account for that share of a seed's total importance.
#' @param exclude_feature_types Feature types to drop before any scoring happens (default \code{"struct"}). struct variables are composite IDs (e.g. \code{polA.group_211.group_2176}, three dot-joined gene/domain identifiers) representing a co-occurrence/structural motif rather than a single molecular entity like the other five scales, and its candidate-variable count (tens of thousands per group) dwarfs the other scales by orders of magnitude — pooling it into this rank_score/contribution machinery would compare a compound signal against five primary ones on an incomparable scale. struct is reserved for post-hoc biological annotation once top clusters are identified, not for scoring/ranking/thresholding here.
#'
#' @returns a tibble of scored top features with their contribution, cumulative contribution, core membership, rank, and rank score within each seed.
#' within a seed
#' contribution is calculated as the importance of a feature divided by the sum of importance scores for all features.
#' cum_contrib is the cumulative contribution when features are sorted in descending order of contribution; features tied on contribution are assigned the same cum_contrib, equal to the cumulative sum through the end of their tied block, so a tie is never split by arbitrary sort order.
#' in_core is TRUE when cum_contrib is at or below core_contribution_threshold; a tied block that pushes the cumulative total past the threshold is excluded in its entirety (conservative: stays at-or-under the threshold rather than overshooting it).
#' Rank is assigned based on the descending order of contribution, with tied contributions receiving the average of the ranks they span, and
#' rank score is calculated as (n_features - rank) / (n_features - 1), where n_features is the total number of features within the same seed.
#' rank score is ranged between 0 and 1, with higher values indicating higher importance.
#'
#' @export
#' @examples
#' scoreFeaturesWithinSeed(all_top_features.parquet)
#'
scoreFeaturesWithinSeed <- function(all_top_features_parquet, core_contribution_threshold = 0.75, exclude_feature_types = "struct") {

  stopifnot(file.exists(all_top_features_parquet))
  scored_top_features <- arrow::read_parquet(normalizePath(all_top_features_parquet)) |>
    dplyr::filter(!shuffled, !feature_type %in% exclude_feature_types) |>
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
      rank = rank(dplyr::desc(contribution), ties.method = "average"),
      n_features = dplyr::n(),
      rank_score = dplyr::if_else(
        n_features > 1,
        (n_features - rank) / (n_features - 1),
        1
      )
    ) |>
    dplyr::arrange(dplyr::desc(contribution), .by_group = TRUE) |>
    dplyr::mutate(running_contrib = cumsum(contribution)) |>
    dplyr::group_by(contribution, .add = TRUE) |>
    dplyr::mutate(
      cum_contrib = max(running_contrib),
      in_core = cum_contrib <= core_contribution_threshold
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-running_contrib)

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
#' mean and median rank score (scaled between 0 and 1 with higher values indicating higher importance),
#' best rank,
#' rank consistency (TRUE if the rank is same across seeds), rank standard deviation (how much the rank_score varies across seeds, computed on the normalized rank_score rather than raw rank so it is comparable across groups with different numbers of features),
#' sign consistency (TRUE if the sign is same across seeds), and sign (the sign of the feature; POS, NEG or MIXED).
#'
#' @export
#' @examples
#' summariseFeaturesAcrossSeeds(scoreFeaturesWithinSeed(all_top_features.parquet))
summariseFeaturesAcrossSeeds <- function(scored_features = scoreFeaturesWithinSeed(all_top_features_parquet)) {
  
  # find max number of seeds
  max_seeds <- scored_features |>
  dplyr::summarise(n_seeds = dplyr::n_distinct(seed)) |>
  dplyr::pull(n_seeds)

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
    seed_ratio = dplyr::n_distinct(seed)/max_seeds,
    # mean_rank = mean(rank, na.rm = TRUE),
    # mean_contribution = mean(contribution, na.rm = TRUE),
    # mean_cum_contrib = mean(cum_contrib, na.rm = TRUE),
    mean_rank_score = mean(rank_score, na.rm = TRUE),
    median_rank = median(rank, na.rm = TRUE),
    median_contribution = median(contribution, na.rm = TRUE),
    median_cum_contrib = median(cum_contrib, na.rm = TRUE),
    median_rank_score = median(rank_score, na.rm = TRUE),
    best_rank = min(rank, na.rm = TRUE),
    rank_consistent = dplyr::n_distinct(rank) == 1,
    rank_score_sd = sd(rank_score, na.rm = TRUE),
    # coefficient of variation (I was trying to see how much closer is SD to mean)
    rank_score_cv = dplyr::if_else(
      mean_rank_score != 0,
      rank_score_sd / mean_rank_score,
      NA_real_
    ),
    in_core_consistent = dplyr::n_distinct(in_core) == 1,
    in_core = if (in_core_consistent) dplyr::first(in_core) else FALSE,
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
#' @param rank_score_quantile The quantile threshold for filtering features based on their median rank scores
#' @param cv_threshold The threshold for filtering features based on the coefficient of variation of rank scores
#' @param cumulative_contribution_threshold The threshold for filtering features based on the median cumulative contribution
#' @param seed_ratio_threshold The threshold for filtering features based on the proportion of seeds in which the feature appears 
#' @param both_subtypes Logical indicating whether to filter features to keep the ones that appears in both binary and count subtypes 
#' @param compare_median_to_sd_rank_score Logical indicating whether to filter features to keep the ones that have higher median rank score than standard deviation
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
                                        cv_threshold = 1,
                                        cumulative_contribution_threshold = 0.75,
                                        #additional filters
                                        seed_ratio_threshold = NULL,
                                        both_subtypes = FALSE,
                                        compare_median_to_sd_rank_score = FALSE 
                                        ) {

top_features <- feature_summary |>
  dplyr::filter(
    if (!is.null(seed_ratio)) seed_ratio == seed_ratio_threshold else TRUE,
    rank_score_cv <= cv_threshold,
    in_core,
    sign_consistent,
    median_cum_contrib <= cumulative_contribution_threshold,
    median_rank_score >= quantile(median_rank_score, rank_score_quantile)
  ) |>
  dplyr::group_by(species, drug_label, drug_or_class, feature_type, variable) |>
    dplyr::mutate( n_subtype = dplyr::n_distinct(feature_subtype),
  subtype_csv = paste(sort(unique(feature_subtype)), collapse = ",") ) |>
  dplyr::ungroup() |>
dplyr::filter(
  if(both_subtypes) (subtype_csv == "binary,counts") else TRUE,
  if(compare_median_to_sd_rank_score) (median_rank_score > rank_score_sd) else TRUE
    )
  
  return(top_features)
}

#' Aggregate mapped features to protein clusters
#'
#' @param top_features The tibble of top features for each drug/class generated from `topFeaturesPerDrugOrClass()`
#' @param cluster_feature_parquet The path to the Parquet file containing the mapping of features to protein clusters
#'
#' @returns a tibble of summarized clusters, including species, drug label, drug or class, cluster, frequency (number of features in the cluster),
#' number of distinct variables, number of distinct feature types, feature types as a comma-separated string,
#' cluster mean rank score (mean of median_rank_score, down-weighted by 1/n_clusters for features that map to multiple clusters so promiscuous domains/COGs don't inflate every cluster they touch),
#' cluster median rank score (median of the same down-weighted score), cluster max rank score (max of the unweighted median_rank_score, i.e. the single strongest feature backing this cluster regardless of its fan-out to other clusters), cluster rank score standard deviation, and cluster best rank.
#' median_rank_score is used throughout (rather than mean_rank_score) for consistency with the seed-noise-robust selection made in `topFeaturesPerDrugOrClass()`.
#'
#' @export
summariseClusters <- function(top_features = topFeaturesPerDrugOrClass(feature_summary),
                               cluster_feature_parquet
                              ) {
  stopifnot(file.exists(cluster_feature_parquet))
  cluster_feature <- arrow::read_parquet(normalizePath(cluster_feature_parquet)) |>
    dplyr::add_count(feature, name = "n_clusters")

  top_clusters <- top_features |>
    dplyr::left_join(
      cluster_feature,
      by = dplyr::join_by(variable == feature),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      n_clusters = dplyr::coalesce(n_clusters, 1L),
      weighted_rank_score = median_rank_score / n_clusters
    ) |>
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
      cluster_mean_rank_score = mean(weighted_rank_score, na.rm = TRUE),
      cluster_median_rank_score = median(weighted_rank_score, na.rm = TRUE),
      cluster_max_rank_score = max(median_rank_score, na.rm = TRUE),
      cluster_rank_score_sd = sd(weighted_rank_score, na.rm = TRUE),
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
#' @returns a tibble of clusters unique to a single drug/class, including drug or class, cluster,
#' cluster name (from protein name annotations), and cluster mean rank score, sorted by cluster mean rank score in descending order.
#'
#' @export
#' @examples
#' findUniqueClusters(summariseClusters(top_features, cluster_feature_parquet), label = "drug", protein_names_parquet)
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
#' Node/edge weights are built from median_rank_score (features) and cluster_mean_rank_score (clusters, itself median-based per `summariseClusters()`), consistent with the seed-noise-robust selection made in `topFeaturesPerDrugOrClass()`.
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
    "median_rank_score"
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
      # mean of median_rank_score across subtype (bin/count) rows for this variable --
      # this is where bin/count reconciliation currently happens (implicitly)
      feature_score = mean(median_rank_score, na.rm = TRUE),
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
      cluster_score = mean(cluster_mean_rank_score, na.rm = TRUE),
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
      score = mean(feature_score, na.rm = TRUE),
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
      score = mean(cluster_score, na.rm = TRUE),
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
      weight = feature_score,
      edge_type = "model_feature"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  cluster_edges <- cluster_table |>
    dplyr::transmute(
      from = model_id,
      to = cluster,
      weight = cluster_score,
      edge_type = "model_cluster"
    ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  feature_cluster_edges <- feature_table |>
    dplyr::left_join(
      cluster_feature |> dplyr::add_count(feature, name = "n_clusters"),
      by = dplyr::join_by(variable == feature),
      relationship = "many-to-many"
    ) |>
    dplyr::filter(!is.na(cluster)) |>
    dplyr::transmute(
      from = variable,
      to = cluster,
      weight = 1 / n_clusters,
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
