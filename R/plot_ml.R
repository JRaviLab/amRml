# This script contains functions to generate ML-related plots for the amR
# project.

#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylim
#' @importFrom graphics barplot
#' @importFrom yardstick pr_curve
NULL

#' Plot a Precision-Recall Curve
#'
#' Generates a precision-recall curve (PRC) for AMR phenotype prediction results.
#' @param test_data_plus_predictions A tibble of test data with added prediction
#' columns (e.g. the output of `predictML()` or `runMLmodels(return_pred=TRUE)`),
#' or a path to a TSV file containing the same.
#'
#' @return A `ggplot2` object showing the precision-recall curve.
#'
#' @details
#' The function uses `yardstick::pr_curve()` to compute the PR curve and then
#' visualizes it using `ggplot2`.
#'
#' @examples
#' data(demo_fit)
#' data(demo_ml_tibble)
#' preds <- predictML(demo_fit, demo_ml_tibble)
#' plotPRC(preds)
#'
#' @export
plotPRC <- function(test_data_plus_predictions) {
  if (is.character(test_data_plus_predictions)) {
    test_data_plus_predictions <- readr::read_tsv(test_data_plus_predictions)
  }
  .checkArgTestDataPlusPredictions(test_data_plus_predictions)
  test_data_plus_predictions <- test_data_plus_predictions |>
    dplyr::mutate(
      genome_drug.resistant_phenotype = factor(
        genome_drug.resistant_phenotype,
        levels = c("Resistant", "Susceptible")
      )
    )

  prc <- yardstick::pr_curve(
    test_data_plus_predictions,
    genome_drug.resistant_phenotype, .pred_Resistant
  ) |>
    ggplot2::ggplot(ggplot2::aes(x = recall, y = precision)) +
    ggplot2::geom_path() +
    ggplot2::ylim(0, 1) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  return(prc)
}

#' Plot a Receiver Operating Characteristic (ROC) Curve
#'
#' Generates a ROC curve for AMR phenotype prediction results.
#'
#' @param test_data_plus_predictions A tibble of test data with prediction
#' columns (e.g. the output of `predictML()` or `runMLmodels(return_pred=TRUE)`),
#' or a path to a TSV file containing the same.
#'
#' @return A ROC curve plotted using `ggplot2::autoplot()`.
#'
#' @examples
#' data(demo_fit)
#' data(demo_ml_tibble)
#' preds <- predictML(demo_fit, demo_ml_tibble)
#' plotROC(preds)
#'
#' @export
plotROC <- function(test_data_plus_predictions) {
  if (is.character(test_data_plus_predictions)) {
    test_data_plus_predictions <- readr::read_tsv(test_data_plus_predictions)
  }
  .checkArgTestDataPlusPredictions(test_data_plus_predictions)
  test_data_plus_predictions <- test_data_plus_predictions |>
    dplyr::mutate(
      genome_drug.resistant_phenotype = factor(
        genome_drug.resistant_phenotype,
        levels = c("Resistant", "Susceptible")
      )
    )

  roc <- yardstick::roc_curve(
    test_data_plus_predictions,
    genome_drug.resistant_phenotype, .pred_Resistant
  ) |>
    ggplot2::autoplot(type = "se") +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  return(roc)
}

#' Plot a Confusion Matrix Heatmap
#'
#' Produces a heatmap visualization of the confusion matrix for AMR predictions.
#'
#' @param test_data_plus_predictions A tibble containing true and predicted
#' phenotype labels (e.g. the output of `predictML()`), or a path to a TSV file
#' containing the same.
#'
#' @return A heatmap (`ggplot2` object) showing the confusion matrix.
#'
#' @examples
#' data(demo_fit)
#' data(demo_ml_tibble)
#' preds <- predictML(demo_fit, demo_ml_tibble)
#' plotCM(preds)
#'
#' @export
plotCM <- function(test_data_plus_predictions) {
  if (is.character(test_data_plus_predictions)) {
    test_data_plus_predictions <- readr::read_tsv(test_data_plus_predictions)
  }
  .checkArgTestDataPlusPredictions(test_data_plus_predictions)
  test_data_plus_predictions <- test_data_plus_predictions |>
    dplyr::mutate(
      genome_drug.resistant_phenotype = factor(
        genome_drug.resistant_phenotype,
        levels = c("Resistant", "Susceptible")
      ),
      .pred_class = factor(
        .pred_class,
        levels = c("Resistant", "Susceptible")
      )
    )
  test_data_plus_predictions |>
    yardstick::conf_mat(
      truth = genome_drug.resistant_phenotype,
      estimate = .pred_class
    ) |>
    ggplot2::autoplot(type = "heatmap")
}

#' Plot Top Feature Importances
#'
#' Creates a bar plot showing the most important features affecting
#' AMR phenotype predictions.
#'
#' @param topfeat A tibble of feature importance scores (e.g. the output of
#' `extractTopFeats()` or `runMLmodels()`), with `Variable`, `Importance`, and
#' `Sign` columns, or a path to a TSV file containing the same.
#' @param n_top_feats Number of top features to display (default: 10).
#'
#' @return A bar plot of variable importance (`ggplot2` object).
#'
#' @examples
#' data(demo_fit)
#' top_feats <- extractTopFeats(demo_fit, n_top_feats = 10)
#' plotTopFeatsVI(top_feats, n_top_feats = 10)
#'
#' @export
plotTopFeatsVI <- function(topfeat, n_top_feats = 10) {
  if (is.character(topfeat)) {
    topfeat <- readr::read_tsv(topfeat)
  }
  .checkArgTibble(topfeat)
  .checkArgNTopFeats(n_top_feats)

  missing_cols <- setdiff(c("Variable", "Importance", "Sign"), colnames(topfeat))
  if (length(missing_cols) > 0) {
    stop(
      "`topfeat` must have `Variable`, `Importance`, and `Sign` columns ",
      "(e.g. the output of `extractTopFeats()`). Missing: ",
      paste(missing_cols, collapse = ", "), "."
    )
  }

  vip <- topfeat |>
    dplyr::slice_max(order_by = Importance, n = n_top_feats) |>
    dplyr::mutate(
      Variable = factor(Variable, levels = rev(Variable)), # preserve order as shown in table
      Sign = factor(Sign, levels = c("POS", "NEG"))
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = Importance, y = Variable, fill = Sign)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c(
        "POS" = "#c6d8d3",
        "NEG" = "#f6c9a1"
      )
    ) +
    ggplot2::labs(
      x = "Importance",
      y = "Features"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  return(vip)
}

#' Compare Baseline Performance With and Without Shuffled Labels
#'
#' Produces a bar plot comparing balanced accuracy for each antibiotic using
#' true AMR labels vs. randomly shuffled labels.
#'
#' @param non_shuffled_label_results A performance tibble with `antibiotic` and
#' `bal_acc` columns (one row per antibiotic), such as an aggregated
#' `performance_tibble` from `runMLPipeline(shuffle_labels = FALSE)`.
#' @param shuffled_label_results The corresponding shuffled-label performance
#' tibble, e.g. from `runMLPipeline(shuffle_labels = TRUE)`.
#'
#' @return A barplot comparing balanced accuracy across antibiotics.
#'
#' @examples
#' non_shuffled <- tibble::tibble(
#'   antibiotic = c("AMP", "CIP", "CRO"),
#'   bal_acc    = c(0.88, 0.81, 0.92)
#' )
#' shuffled <- tibble::tibble(
#'   antibiotic = c("AMP", "CIP", "CRO"),
#'   bal_acc    = c(0.52, 0.49, 0.55)
#' )
#' plotBaselineComparison(non_shuffled, shuffled)
#' @export
plotBaselineComparison <- function(
  non_shuffled_label_results,
  shuffled_label_results
) {
  .checkArgTibble(non_shuffled_label_results)
  .checkArgTibble(shuffled_label_results)

  required_cols <- c("antibiotic", "bal_acc")
  for (results in list(non_shuffled_label_results, shuffled_label_results)) {
    missing_cols <- setdiff(required_cols, colnames(results))
    if (length(missing_cols) > 0) {
      stop(
        "Baseline comparison inputs must have `antibiotic` and `bal_acc` ",
        "columns. Missing: ", paste(missing_cols, collapse = ", "), "."
      )
    }
  }

  drugs <- non_shuffled_label_results |>
    dplyr::pull(antibiotic)

  non_shuffled_bal_acc <- non_shuffled_label_results |>
    dplyr::pull(bal_acc)

  shuffled_bal_acc <- shuffled_label_results |>
    dplyr::pull(bal_acc)

  bal_acc_matrix <- matrix(c(non_shuffled_bal_acc, shuffled_bal_acc),
    nrow = 2, byrow = TRUE
  )

  colnames(bal_acc_matrix) <- drugs
  rownames(bal_acc_matrix) <- c("Non-Shuffled Labels", "Shuffled Labels")

  baseline_comparison_barplot <- barplot(bal_acc_matrix,
    beside = TRUE,
    legend.text = TRUE, col = c("skyblue", "lightpink"),
    ylab = "Balanced accuracy", xlab = "Antibiotic"
  )

  return(baseline_comparison_barplot)
}

#' Plot top features' Fisher's significance
#'
#' This function visualizes all features from \code{runFishers()} ranked by
#' BH-adjusted p-value and explicitly highlights those that pass the
#' significance threshold. Optionally, the top N most significant features
#' can be labeled.
#'
#' @param fisher_df A data frame returned by \code{runFishers()} containing
#'   at minimum the columns:
#'   \itemize{
#'     \item gene
#'     \item adj_p_value
#'     \item sig_after_bh
#'   }
#' @param alpha BH-adjusted significance threshold. Default is 0.05.
#' @param label_top_n Number of top-ranked features to label.
#'   Default is 5. Set to 0 to disable labeling.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' Each point represents a feature.
#' Color explicitly encodes whether a feature passes the BH threshold.
#' Labels are applied only to the top-ranked features to preserve clarity.
#'
#' @examples
#' long <- tibble::tibble(
#'   genome_id = rep(paste0("g", 1:10), each = 2),
#'   feature_id = rep(c("gene_a", "gene_b"), 10),
#'   value = c(
#'     1, 0, 1, 0, 1, 1, 1, 1, 0, 1,
#'     0, 0, 0, 1, 0, 1, 0, 1, 0, 0
#'   ),
#'   genome_drug.resistant_phenotype = rep(
#'     rep(c("Resistant", "Susceptible"), each = 5),
#'     each = 2
#'   )
#' )
#' tmp <- tempfile(fileext = ".parquet")
#' arrow::write_parquet(long, tmp)
#' fisher_results <- runFishers(tmp, Q = 0.05)
#' plotFishers(fisher_results, alpha = 0.05, label_top_n = 2)
#'
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @export
plotFishers <- function(
  fisher_df,
  alpha = 0.05,
  label_top_n = 5
) {
  required_cols <- c("gene", "adj_p_value", "sig_after_bh")
  missing_cols <- setdiff(required_cols, colnames(fisher_df))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  plot_df <- fisher_df |>
    dplyr::arrange(adj_p_value) |>
    dplyr::mutate(
      rank = dplyr::row_number(),
      neg_log10_adj_p = -log10(adj_p_value),
      significance = ifelse(sig_after_bh, "Significant", "Not significant")
    )

  p <- ggplot(plot_df, aes(x = rank, y = neg_log10_adj_p, color = significance)) +
    geom_point(size = 1.6, alpha = 0.8) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      color = "grey40"
    ) +
    scale_color_manual(
      values = c(
        "Significant" = "#0072B2",
        "Not significant" = "#BDBDBD"
      ),
      guide = "none"
    ) +
    labs(
      x = "Feature rank",
      y = expression(-log[10]("BH-adjusted p-value")),
      title = "Ranked feature significance"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold")
    )

  if (label_top_n > 0) {
    label_df <- plot_df |>
      dplyr::slice_head(n = label_top_n)

    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(label = gene),
        size = 3.5,
        color = "black",
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.3
      )
  }

  return(p)
}

# Read `metadata.parquet` from a directory, with a clear error if missing.
.readMetadata <- function(metadata_path) {
  if (!is.character(metadata_path) || length(metadata_path) != 1) {
    stop("`metadata_path` must be a single directory path (character).")
  }
  fp <- file.path(metadata_path, "metadata.parquet")
  if (!file.exists(fp)) {
    stop("Could not find `metadata.parquet` in directory: ", metadata_path)
  }
  arrow::read_parquet(fp)
}

# Accept a data frame directly, or read `filename` from a directory path.
# `arg_name` is used only to produce informative error messages.
.readPerfInput <- function(x, filename, arg_name) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (!is.character(x) || length(x) != 1) {
    stop("`", arg_name, "` must be a data frame or a single directory path.")
  }
  fp <- file.path(x, filename)
  if (!file.exists(fp)) {
    stop("Could not find `", filename, "` in `", arg_name, "` directory: ", x)
  }
  arrow::read_parquet(fp)
}

#' Plot drug phenotype distribution
#'
#' Reads metadata and generates a stacked bar plot showing counts of resistant
#' and susceptible phenotypes per antibiotic.
#'
#' @param metadata_path Character. Path to directory containing `metadata.parquet`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' plotDrugDist(metadata_path = system.file("extdata", "Shigella_flexneri", package = "amRml"))
plotDrugDist <- function(metadata_path = ".") {
  metadata <- .readMetadata(metadata_path)

  drug_dist <- metadata |>
    dplyr::distinct(
      genome.genome_id,
      genome_drug.antibiotic,
      drug_abbr,
      genome_drug.resistant_phenotype
    ) |>
    dplyr::count(
      genome_drug.antibiotic,
      drug_abbr,
      genome_drug.resistant_phenotype
    ) |>
    dplyr::group_by(genome_drug.antibiotic, drug_abbr) |>
    dplyr::mutate(total = sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      label = paste0(genome_drug.antibiotic, " (", drug_abbr, ")"),
      label = forcats::fct_reorder(label, total)
    )

  p <- ggplot2::ggplot(
    drug_dist,
    ggplot2::aes(
      x = label,
      y = n,
      fill = genome_drug.resistant_phenotype
    )
  ) +
    ggplot2::geom_col(color = "black", width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        "Resistant" = "#d4872a",
        "Susceptible" = "#5b8db8"
      ),
      name = "Phenotype"
    ) +
    ggplot2::labs(
      x = "Antibiotic",
      y = "Number of unique genomes"
    ) +
    ggplot2::theme_classic(base_size = 14)

  p
}

#' Plot drug-level model performance
#'
#' Generates heatmaps and ridge plots summarizing model performance (MCC)
#' across drugs and feature types.
#'
#' @param metadata_path Character. Path to a directory containing `metadata.parquet`.
#' @param performance_path A performance tibble (with `drug_label`, `shuffled`,
#' `drug_or_class`, `feature_type`, `feature_subtype`, and `mcc` columns), or a
#' directory path containing `all_perf.parquet`.
#'
#' @return A patchwork ggplot object combining multiple panels.
#' @export
#'
#' @examples
#' # Several models per drug x feature type so the ridge densities can be drawn.
#' performance <- tidyr::expand_grid(
#'   drug_or_class = c("AMP", "CIP", "NAL"),
#'   feature_type = c("genes", "proteins"),
#'   feature_subtype = c("binary", "counts"),
#'   replicate = 1:8
#' )
#' performance$drug_label <- "drug"
#' performance$shuffled <- FALSE
#' performance$mcc <- 0.6 + 0.3 * sin(seq_len(nrow(performance)))
#' plotDrugPerf(
#'   metadata_path = system.file("extdata", "Shigella_flexneri", package = "amRml"),
#'   performance_path = performance
#' )
plotDrugPerf <- function(metadata_path = ".", performance_path = ".") {
  metadata <- .readMetadata(metadata_path)
  performance <- .readPerfInput(
    performance_path, "all_perf.parquet", "performance_path"
  )

  plot_df <- metadata |>
    dplyr::distinct(genome.genome_id, genome_drug.antibiotic, drug_abbr) |>
    dplyr::count(genome_drug.antibiotic, drug_abbr, name = "total")

  ######################## drug performances #################################
  median_drug <- performance |>
    dplyr::filter(
      drug_label == "drug",
      !shuffled # keep real models; remove if you want both
    ) |>
    dplyr::group_by(drug_or_class, feature_type, feature_subtype) |>
    dplyr::summarise(median_mcc = median(mcc, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(plot_df, by = c("drug_or_class" = "drug_abbr")) |>
    dplyr::mutate(drug_or_class = reorder(drug_or_class, total))

  drug_p1 <- ggplot2::ggplot(
    median_drug,
    ggplot2::aes(
      x = feature_type,
      y = drug_or_class,
      fill = median_mcc
    )
  ) +
    ggplot2::geom_tile(color = "grey90", width = 0.9) +
    ggplot2::scale_fill_gradientn(
      colors = c(
        "#EAF2FF", # very light blue
        "#5F84C9", # medium blue
        "#0F2A5A" # dark blue
      ),
      values = scales::rescale(c(-1, 0, 1)),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      labels = scales::label_number(accuracy = 0.01),
      name = "Median MCC"
    ) +
    ggplot2::labs(x = "Feature type") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    ggplot2::coord_fixed()

  median_feature <- performance |>
    dplyr::filter(
      drug_label == "drug",
      !shuffled
    ) |>
    dplyr::group_by(drug_or_class, feature_type) |>
    dplyr::summarise(median_mcc = median(mcc, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(plot_df, by = c("drug_or_class" = "drug_abbr")) |>
    dplyr::mutate(drug_or_class = reorder(drug_or_class, total))


  feat_pal <- c(
    "args"     = "#56B4E9", # sky blue
    "cogs"     = "#E69F00", # orange
    "genes"    = "#009E73", # bluish green
    "domains"  = "#F0E442", # yellow
    "proteins" = "#CC79A7", # reddish purple
    "struct"   = "#D55E00" # vermillion
  )

  rc_perf <- ggplot2::ggplot(
    median_feature |>
      dplyr::distinct(
        drug_or_class,
        feature_type, median_mcc
      ),
    ggplot2::aes(x = median_mcc, y = drug_or_class)
  ) +
    ggridges::geom_density_ridges(
      scale = 0.75,
      rel_min_height = 0.01,
      alpha = 0.4,
      fill = "grey90",
      colour = "grey70"
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(height = 0.1),
      size = 1.5,
      alpha = 0.8,
      ggplot2::aes(color = feature_type)
    ) +
    ggplot2::scale_color_manual(values = feat_pal, name = "Feature type") +
    ggplot2::stat_summary(
      fun = median,
      geom = "point",
      size = 2,
      color = "black"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black")
    )

  final_plot <- patchwork::wrap_plots(
    drug_p1, rc_perf,
    widths = c(2, 2), # adjust proportions
    guides = "collect"
  ) &
    ggplot2::theme(
      legend.position = "bottom"
    )

  final_plot
}

#' Plot cross-drug prediction as a heatmap
#'
#' Creates a heatmap showing cross-drug model performance (MCC), where models
#' trained on one drug are evaluated on another.
#'
#' @param cross_test_performance_path A cross-drug performance tibble (with
#' `drug`, `test_drug`, and `mcc` columns), or a directory path
#' containing `cross_drug_perf.parquet`.
#' @param drug_performance_path A performance tibble (with `drug_label`,
#' `drug_or_class`, and `mcc` columns), or a directory path containing
#' `all_perf.parquet`.
#' @param metadata_path Character. Path to a directory containing `metadata.parquet`.
#'
#' @return A ComplexHeatmap object.
#' @export
#'
#' @examples
#' cross_drug <- tibble::tibble(
#'   drug = c("AMP", "AMP", "CIP", "CIP", "NAL", "NAL"),
#'   test_drug = c("CIP", "NAL", "AMP", "NAL", "AMP", "CIP"),
#'   mcc = c(0.3, 0.2, 0.4, 0.25, 0.15, 0.35)
#' )
#' performance <- tibble::tibble(
#'   drug_label = "drug",
#'   drug_or_class = c("AMP", "CIP", "NAL"),
#'   mcc = c(0.8, 0.7, 0.6)
#' )
#' plotCrossDrug(
#'   cross_test_performance_path = cross_drug,
#'   drug_performance_path = performance,
#'   metadata_path = system.file("extdata", "Shigella_flexneri", package = "amRml")
#' )
plotCrossDrug <- function(
  cross_test_performance_path = ".",
  drug_performance_path = ".",
  metadata_path = "."
) {
  cross_drug <- .readPerfInput(
    cross_test_performance_path, "cross_drug_perf.parquet",
    "cross_test_performance_path"
  )
  performance <- .readPerfInput(
    drug_performance_path, "all_perf.parquet", "drug_performance_path"
  )
  metadata <- .readMetadata(metadata_path)

  heatmap_df <- cross_drug |>
    dplyr::filter(!is.na(drug), !is.na(test_drug)) |>
    dplyr::group_by(drug, test_drug) |>
    dplyr::summarise(median_mcc = median(mcc, na.rm = TRUE), .groups = "drop")

  same_drugs <- performance |>
    dplyr::filter(
      drug_label == "drug",
      drug_or_class %in% (cross_drug |>
        dplyr::distinct(drug) |>
        dplyr::pull())
    ) |>
    dplyr::group_by(drug_or_class) |>
    dplyr::summarise(median_mcc = median(mcc, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(test_drug = drug_or_class) |>
    dplyr::distinct(drug_or_class, test_drug, median_mcc) |>
    dplyr::rename(drug = drug_or_class)

  heatmap_df <- heatmap_df |>
    dplyr::add_row(same_drugs) |>
    dplyr::left_join(
      metadata |>
        dplyr::distinct(drug_abbr, class_abbr),
      by = c("drug" = "drug_abbr")
    ) |>
    dplyr::rename("drug_class" = "class_abbr") |>
    dplyr::left_join(
      metadata |>
        dplyr::distinct(drug_abbr, class_abbr),
      by = c("test_drug" = "drug_abbr")
    )

  # Row annotation (already similar to what you did)
  annotation_row <- heatmap_df |>
    dplyr::distinct(drug, drug_class) |>
    tibble::column_to_rownames("drug")

  # Column annotation
  annotation_col <- heatmap_df |>
    dplyr::distinct(test_drug, class_abbr) |>
    tibble::column_to_rownames("test_drug")

  mat <- heatmap_df |>
    dplyr::select(drug, test_drug, median_mcc) |>
    tidyr::pivot_wider(names_from = test_drug, values_from = median_mcc) |>
    tibble::column_to_rownames("drug") |>
    as.matrix()

  row_order <- heatmap_df |>
    dplyr::distinct(drug, drug_class) |>
    dplyr::arrange(drug_class, drug) |>
    dplyr::pull(drug)

  col_order <- heatmap_df |>
    dplyr::distinct(test_drug, class_abbr) |>
    dplyr::arrange(class_abbr, test_drug) |>
    dplyr::pull(test_drug)

  mat <- mat[row_order, col_order]

  # Align annotations
  annotation_row <- annotation_row[row_order, , drop = FALSE]
  annotation_col <- annotation_col[col_order, , drop = FALSE]


  # Collect all classes from both row and column
  classes <- base::union(
    annotation_row$drug_class,
    annotation_col$class_abbr
  )

  # One named color vector shared by row and column annotations.
  class_colors <- stats::setNames(
    colorRampPalette(RColorBrewer::brewer.pal(9, "Pastel1"))(length(classes)),
    classes
  )


  heat_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100)
  max_val <- max(abs(mat), na.rm = TRUE)

  # ---- Convert annotations ----
  ha_row <- ComplexHeatmap::rowAnnotation(
    drug_class = annotation_row$drug_class,
    col = list(drug_class = class_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  ha_col <- ComplexHeatmap::HeatmapAnnotation(
    class_abbr = annotation_col$class_abbr,
    col = list(class_abbr = class_colors),
    show_annotation_name = FALSE, na_col = "grey3",
    annotation_legend_param = list(labels_gp = grid::gpar(fontsize = 10))
  )

  # ---- Color function (instead of breaks + palette) ----
  col_fun <- circlize::colorRamp2(
    seq(-max_val, max_val, length.out = length(heat_colors)),
    heat_colors
  )
  # ---- Heatmap ----
  cross_drug_hm <- ComplexHeatmap::Heatmap(
    mat,
    name = "median_mcc",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = row_order,
    column_order = col_order,
    left_annotation = ha_row,
    top_annotation = ha_col,
    width = grid::unit(ncol(mat) * 1.5, "cm"),
    height = grid::unit(nrow(mat) * 1.5, "cm"),
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = "tested on",
    row_title = "trained on",
    column_title_side = "bottom",
    row_title_side = "right",
    row_names_gp = grid::gpar(fontsize = 14),
    column_names_gp = grid::gpar(fontsize = 14),
    column_names_rot = 45,

    # remove borders like pheatmap
    rect_gp = grid::gpar(col = NA),

    # legends
    show_heatmap_legend = TRUE,
    use_raster = FALSE
  )

  cross_drug_hm
}

#' Plot stratified model performance
#'
#' Visualizes model performance (MCC) stratified by year or country,
#' comparing within-group vs cross-group evaluation.
#'
#' @param year_or_country Character. Either "year" or "country".
#' @param stratified_performance_path Character. Path to stratified performance files.
#' @param stratified_cross_performance_path Character. Path to cross-stratified performance files.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' plotStratifiedPerf("year",
#'   stratified_performance_path = "data/Campylobacter/ML_year_performance",
#'   stratified_cross_performance_path = "data/Campylobacter/cross_test_ML_year_performance"
#' )
#' }
plotStratifiedPerf <- function(year_or_country = "year",
                               stratified_performance_path = ".",
                               stratified_cross_performance_path = ".") {
  if (!year_or_country %in% c("year", "country")) {
    stop("`year_or_country` must be either \"year\" or \"country\".")
  }

  perf_file <- file.path(
    stratified_performance_path,
    paste0(year_or_country, "_perf.parquet")
  )
  cross_file <- file.path(
    stratified_cross_performance_path,
    paste0("cross_", year_or_country, "_perf.parquet")
  )
  if (!file.exists(perf_file)) {
    stop("Could not find stratified performance file: ", perf_file)
  }
  if (!file.exists(cross_file)) {
    stop("Could not find cross-stratified performance file: ", cross_file)
  }

  perf <- arrow::read_parquet(perf_file)
  cross_test <- arrow::read_parquet(cross_file)

  all <- perf |>
    dplyr::mutate(strat_value_test = strat_value) |>
    dplyr::select(
      drug_label, drug_or_class,
      strat_value, strat_value_test, feature_type, feature_subtype, mcc
    ) |>
    dplyr::bind_rows(cross_test |>
      dplyr::select(
        drug_label, drug_or_class,
        strat_value, strat_value_test, feature_type,
        feature_subtype, mcc
      )) |>
    dplyr::mutate(category = dplyr::if_else(
      strat_value == strat_value_test, "same", "different"
    ))

  fill_vals <- c(
    "same" = "#b3cde3",
    "different" = "#fbb4ae"
  )

  plot <- ggplot2::ggplot(
    all |>
      dplyr::filter(drug_label == "drug", !is.na(mcc)),
    ggplot2::aes(x = mcc, y = drug_or_class, fill = category)
  ) +
    ggridges::geom_density_ridges(
      alpha = 0.5,
      scale = 1,
      rel_min_height = 0.01,
      position = "identity"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::scale_fill_manual(values = fill_vals) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = if (year_or_country == "year") {
        "Performance of AMR drug models with temporal holdouts"
      } else {
        "Performance of AMR drug models with geographical holdouts"
      },
      x = "MCC",
      y = "Drug",
      fill = "Tested on"
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black", size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 14, colour = "black"),
      axis.title.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
  plot
}

#' Plot multi-drug resistance (MDR) model performance
#'
#' Generates a violin plot of performance across feature types and a
#' prediction confusion-style tile plot for MDR models.
#'
#' @param MDR_performance_path Character. Path to a directory containing
#' `MDR_perf.parquet`.
#' @param MDR_top_feature_path Character. Path to a directory containing
#' `MDR_top_features.parquet`. Currently unused (reserved for an upcoming
#' feature-cluster panel).
#' @param MDR_pred_path Character. Path to a directory containing
#' `MDR_pred.parquet`.
#'
#' @return A named list of two `ggplot` objects: `performance` (violin plot)
#' and `prediction` (predicted-vs-true tile plot).
#' @export
#'
#' @examples
#' \dontrun{
#' plotMDR(
#'   MDR_performance_path = "data/Campylobacter/MDR_ML_performance",
#'   MDR_top_feature_path = "data/Campylobacter/MDR_ML_top_features",
#'   MDR_pred_path = "data/Campylobacter/MDR_ML_pred"
#' )
#' }
plotMDR <- function(MDR_performance_path = ".", MDR_top_feature_path = ".",
                    MDR_pred_path = ".") {
  perf_file <- file.path(MDR_performance_path, "MDR_perf.parquet")
  pred_file <- file.path(MDR_pred_path, "MDR_pred.parquet")
  if (!file.exists(perf_file)) {
    stop("Could not find `MDR_perf.parquet` in: ", MDR_performance_path)
  }
  if (!file.exists(pred_file)) {
    stop("Could not find `MDR_pred.parquet` in: ", MDR_pred_path)
  }

  MDR_perf <- arrow::read_parquet(perf_file)

  # ---- Violin plot ----
  perf_plot <- ggplot2::ggplot(
    MDR_perf,
    ggplot2::aes(x = feature_type, y = mcc)
  ) +

    # violins (overall distribution per feature type)
    ggplot2::geom_violin(fill = "grey85", color = NA, alpha = 0.8) +

    # points (colored by binary vs counts)
    ggplot2::geom_jitter(
      ggplot2::aes(color = feature_subtype),
      width = 0.12, size = 2, alpha = 0.8
    ) +
    ggplot2::scale_color_manual(values = c(
      "binary" = "#7B9CB5",
      "counts" = "#CC8644"
    )) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      # title = "MDR model performances",
      # subtitle = "Violin = distribution per feature type; points = binary vs counts",
      x = "Feature type",
      y = "MCC",
      color = "Feature\nsubtype"
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black", size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 14, colour = "black"),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "none",
      title = ggplot2::element_text(face = "bold"),
      panel.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(), # remove vertical lines
      panel.grid.major.y = ggplot2::element_line(color = "grey80"), # keep horizontal lines

      axis.line = ggplot2::element_line(color = "black")
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1))

  MDR_pred <- arrow::read_parquet(pred_file) |>
    dplyr::mutate(
      diff_top2 = purrr::pmap_dbl(dplyr::across(dplyr::contains(".pred") & dplyr::where(is.numeric)), function(...) {
        x <- c(...)
        sx <- sort(x, decreasing = TRUE)
        sx[1] - sx[2]
      }) # Difference between prediction probabilities of top two classes
    ) |>
    dplyr::select(
      genome_id, resistant_classes, .pred_class, diff_top2,
      feature_type, feature_subtype, seed
    ) |>
    dplyr::group_by(resistant_classes, .pred_class, feature_type) |>
    dplyr::summarise(mean_margin = mean(diff_top2), n = n(), .groups = "drop") |>
    dplyr::group_by(resistant_classes, feature_type) |> # normalize within true class
    dplyr::mutate(sum = sum(n), prop = n / sum(n)) |>
    dplyr::ungroup()

  MDR_pred_plot <- ggplot2::ggplot(
    MDR_pred,
    ggplot2::aes(
      x = resistant_classes,
      y = .pred_class
    )
  ) +
    ggplot2::geom_tile(ggplot2::aes(fill = prop)) +
    ggplot2::geom_point(ggplot2::aes(size = mean_margin), color = "black") +
    ggplot2::facet_wrap(~feature_type) +
    ggplot2::scale_fill_distiller(
      palette = "RdBu",
      direction = 1, # flip with -1 if needed
      name = "Prediction proportion"
    ) +
    ggplot2::labs(x = "true class", y = "predicted class") +
    ggplot2::scale_size(range = c(1, 6), name = "Mean margin") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black", size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "right",
      title = ggplot2::element_text(face = "bold")
    )

  list(
    performance = perf_plot,
    prediction = MDR_pred_plot
  )
}

#' Compare shuffled vs real model performance
#'
#' Creates boxplots comparing performance (MCC) between real and shuffled labels
#' across feature types.
#'
#' @param metadata_path Character. Unused; retained for backward compatibility.
#' @param performance_path A performance tibble (with `feature_type`,
#' `feature_subtype`, `mcc`, and `shuffled` columns), or a directory path
#' containing `all_perf.parquet`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' performance <- tibble::tibble(
#'   feature_type = rep(c("genes", "proteins"), each = 4),
#'   feature_subtype = rep(c("binary", "counts"), times = 4),
#'   mcc = c(0.7, 0.6, 0.65, 0.55, 0.1, 0.05, 0.08, 0.02),
#'   shuffled = rep(c(FALSE, TRUE), each = 4)
#' )
#' plotShuffleVsReal(performance_path = performance)
plotShuffleVsReal <- function(metadata_path = ".", performance_path = ".") {
  performance <- .readPerfInput(
    performance_path, "all_perf.parquet", "performance_path"
  )

  performance |>
    dplyr::mutate(
      shuffled_label = dplyr::if_else(shuffled, "shuffled", "real")
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = feature_subtype, y = mcc, fill = shuffled_label)) +
    ggplot2::geom_boxplot(
      width = 0.55, outlier.size = 0.8, outlier.alpha = 0.4,
      outlier.color = "grey50", linewidth = 0.4
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    ggplot2::scale_fill_manual(
      values = c("real" = "#7B9CB5", "shuffled" = "#C4B8A8"),
      name = NULL
    ) +
    ggplot2::scale_y_continuous(limits = c(-0.2, 1), breaks = seq(-0.2, 1, 0.2)) +
    ggplot2::facet_wrap(~feature_type, nrow = 1) +
    ggplot2::theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#E5E2D9", linewidth = 0.4),
      strip.text = element_text(color = "grey30", face = "bold", size = 10),
      # strip.background = element_rect(fill = "#EEEAE0", color = NA),
      axis.text = element_text(color = "grey45"),
      legend.position = "top",
      legend.text = element_text(color = "grey40", size = 10)
    ) +
    labs(
      x = NULL, y = "MCC"
    )
}

#' Plot top contributing feature clusters
#'
#' Identifies top contributing clusters across feature types and drugs,
#' and visualizes their relative contributions.
#'
#' @param top_feat_path Character. Path to `all_top_features.parquet`.
#' @param cluster_feature_path Character. Path to `cluster_feature.parquet`.
#' @param protein_names_path Character. Path to `protein_names.parquet`.
#' @param top_n Integer. Number of top features to retain per model.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' plotTopClusters(
#'   top_feat_path = "data/Campylobacter/ML_top_features",
#'   cluster_feature_path = "data/Campylobacter/",
#'   protein_names_path = "data/Campylobacter/",
#'   top_n = 10
#' )
#' }
plotTopClusters <- function(top_feat_path = ".", cluster_feature_path = ".",
                            protein_names_path = ".", top_n = 10) {
  top_feat <- .readPerfInput(
    top_feat_path, "all_top_features.parquet", "top_feat_path"
  )
  cluster_feature <- .readPerfInput(
    cluster_feature_path, "cluster_feature.parquet", "cluster_feature_path"
  )
  protein_names <- .readPerfInput(
    protein_names_path, "protein_names.parquet", "protein_names_path"
  )

  # which clusters appear in top n across feature types per drug
  # join top features with cluster mapping, filter out struct and shuffled
  top_feat_clean <- top_feat |>
    dplyr::filter(!shuffled, feature_type != "struct", drug_label == "drug") |>
    dplyr::group_by(drug_or_class, feature_type, feature_subtype, seed) |>
    dplyr::slice_max(Importance, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(Variable = gsub(".NCBIFAM", "", Variable)) |>
    dplyr::mutate(Variable = gsub("^X", "", Variable)) |>
    dplyr::mutate(Variable = dplyr::if_else(
      feature_type == "domains", gsub("_.*", "", Variable), Variable
    )) |>
    dplyr::mutate(Variable = dplyr::if_else(
      feature_type == "proteins", gsub("fig.", "fig|", Variable), Variable
    )) |>
    dplyr::left_join(cluster_feature, by = c("Variable" = "feature")) |>
    dplyr::mutate(
      cluster = dplyr::coalesce(cluster, Variable), # fallback to Variable if no match
      Importance_signed = dplyr::if_else(Sign == "NEG", -Importance, Importance)
    )

  shared_mat <- top_feat_clean |>
    dplyr::group_by(drug_or_class, feature_type, cluster) |>
    dplyr::summarise(abs_imp = median(Importance, na.rm = TRUE), .groups = "drop") |>
    # convert to contribution within each feature_type
    dplyr::group_by(drug_or_class, feature_type) |>
    dplyr::mutate(contribution = abs_imp / sum(abs_imp, na.rm = TRUE)) |>
    # pick top n contributors
    dplyr::slice_max(contribution, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::add_count(drug_or_class, cluster, name = "n_feat_types") |>
    dplyr::left_join(protein_names, by = c("cluster" = "proteinID")) |>
    dplyr::mutate(
      proteinName = stringr::str_trunc(proteinName, 50),
      proteinName = forcats::fct_reorder(proteinName, n_feat_types)
    )

  feat_plot <- ggplot2::ggplot(
    shared_mat,
    ggplot2::aes(
      x = feature_type,
      y = proteinName,
      fill = contribution
    )
  ) +
    ggplot2::geom_tile(color = "#FAFAF7", linewidth = 0.5, width = 0.9, height = 0.9) +
    # coord_fixed() +
    ggplot2::scale_fill_distiller(
      palette = "RdPu",
      direction = 1,
      name = "contribution",
      na.value = "#EEEAE0"
    ) +
    ggplot2::facet_wrap(~drug_or_class, scales = "free_y") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E5E2D9", linewidth = 0.4),
      strip.text = ggplot2::element_text(color = "grey30", face = "bold", size = 10),
      strip.background = ggplot2::element_rect(fill = "#EEEAE0", color = NA),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(color = "black", angle = 30, hjust = 1, size = 6),
      axis.text.y = ggplot2::element_text(color = "black", size = 6),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(color = "grey40", size = 10),
      legend.title = ggplot2::element_text(color = "grey40", size = 10)
    )

  feat_plot
}
