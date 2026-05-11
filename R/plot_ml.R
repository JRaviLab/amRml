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
#' @importFrom tune extract_fit_parsnip
#' @importFrom vip vip
#' @importFrom yardstick pr_curve
#' @importFrom graphics barplot
NULL

#' plotPRC()
#'
#' Plots the precision-recall curve given a set of test data plus predicted AMR
#' phenotypes.
#'
#' @param test_data_plus_predictions Test data (tibble) with an added column for
#' predicted phenotype labels, such as the output of `predict()`.
#' @return A precision-recall curve as a `ggplot2` object
#' @export
plotPRC <- function(test_data_plus_predictions) {
  .checkArgTestDataPlusPredictions(test_data_plus_predictions)

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

#' plotTopFeatsVI()
#'
#' Generates a plot showing the top features and their variable importance
#' scores.
#'
#' @param fit Best model fit, such as the output of `fitBestModel()`
#' @param n_top_feats [num] Number of top features to plot
#' @return Variable importance plot (a `ggplot2` object)
#' @export
plotTopFeatsVI <- function(fit, n_top_feats = 10) {
  .checkArgWflow(fit)
  .checkArgNTopFeats(n_top_feats)

  vip <- fit |>
    tune::extract_fit_parsnip() |>
    vip::vip(num_features = n_top_feats) +
    ggplot2::xlab("Top features") +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  return(vip)
}

#' plotDefaultEval()
#'
#' Plots performance metric or runtime vs. training data proportion or number
#' of cross-validation folds, colored by model.
#'
#' @param default_eval_tibble Output of `findOptimalMLDefaults()`
#' @param x_default_eval [chr] x value of default evaluation plot: "train_prop"
#' or "n_fold"
#' @param y_default_eval [chr] y value of default evaluation plot. It can be
#' "avg_runtime_sec" or one of the following performance metrics:
#' "avg_f1_score", "avg_log2_apop", "avg_bal_acc", "avg_mcc", or "avg_nmcc"
#' @param xlab [chr] Label for x axis
#' @param ylab [chr] Label for y axis
#' @return A `ggplot2` scatterplot (performance metric or runtime vs.
#' `train_prop` or `n_fold`), colored by model
#' @export
plotDefaultEval <- function(
  default_eval_tibble, x_default_eval = "train_prop",
  y_default_eval = "avg_f1_score", xlab = "Train data proportion",
  ylab = "Average F1 score"
) {
  .checkArgTibble(default_eval_tibble)
  .checkArgXDefaultEval(x_default_eval)
  .checkArgYDefaultEval(y_default_eval)
  .checkArgXYLabs(xlab = xlab, ylab = ylab)

  if (x_default_eval == "n_fold") {
    default_eval_tibble <- default_eval_tibble |>
      dplyr::filter(train_prop == 0.8)
  } else {
    default_eval_tibble <- default_eval_tibble |>
      dplyr::filter(train_prop != 0.8)
  }

  default_eval_plot <- ggplot2::ggplot(
    default_eval_tibble,
    ggplot2::aes(
      x = unlist(default_eval_tibble[x_default_eval]),
      y = unlist(default_eval_tibble[y_default_eval]), color = model
    )
  ) +
    ggplot2::geom_line(size = 1.5) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 1.5),
      axis.ticks = ggplot2::element_line(linewidth = 1.5, colour = "black"),
      axis.text = ggplot2::element_text(size = 16, colour = "black"),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 16),
      legend.title = ggplot2::element_text(size = 16, face = "bold")
    ) +
    ggplot2::labs(x = xlab, y = ylab, color = "Model")

  return(default_eval_plot)
}

#' plotBaselineComparison()
#'
#' Generates a bar plot that compares model performance with and without
#' randomly shuffled AMR phenotype labels.
#'
#' @importFrom graphics barplot
#'
#' @param non_shuffled_label_results Output of `runMLPipeline()`
#' (`shuffle_labels = FALSE`)
#' @param shuffled_label_results Output of `runMLPipeline()`
#' (`shuffle_labels = TRUE`)
#' @return A bar plot with balanced accuracy comparisons per antibiotic
#' @export
plotBaselineComparison <- function(
  non_shuffled_label_results,
  shuffled_label_results
) {
  .checkArgTibble(non_shuffled_label_results)
  .checkArgTibble(shuffled_label_results)

  drugs <- non_shuffled_label_results |>
    dplyr::select(antibiotic) |>
    dplyr::pull()

  non_shuffled_bal_acc <- non_shuffled_label_results |>
    dplyr::select(bal_acc) |>
    dplyr::pull()

  shuffled_bal_acc <- shuffled_label_results |>
    dplyr::select(bal_acc) |>
    dplyr::pull()

  bal_acc_matrix <- matrix(c(non_shuffled_bal_acc, shuffled_bal_acc),
    nrow = 2, byrow = TRUE
  )

  colnames(bal_acc_matrix) <- drugs
  rownames(bal_acc_matrix) <- c("Non-shuffled labels", "Shuffled labels")

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
#' \dontrun{
#' plotFishers(fisher_results)
#' plotFishers(fisher_results, label_top_n = 0)
#' }
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

  plot_df <- fisher_df %>%
    dplyr::arrange(adj_p_value) %>%
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
    label_df <- plot_df %>%
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
