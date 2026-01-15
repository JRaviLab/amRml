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
  
  prc <- yardstick::pr_curve(test_data_plus_predictions, 
    genome_drug.resistant_phenotype, .pred_Resistant) |> 
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
  .checkArgWflow(fit); .checkArgNTopFeats(n_top_feats)
  
  vip <- fit |> tune::extract_fit_parsnip() |> 
    vip::vip(num_features = n_top_feats) + 
    ggplot2::xlab("Top Features") + 
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
#' "avg_f1_score", "avg_log2_apop", "avg_bal_acc", or "avg_nmcc"
#' @param xlab [chr] Label for x axis
#' @param ylab [chr] Label for y axis
#' @return A `ggplot2` scatterplot (performance metric or runtime vs. 
#' `train_prop` or `n_fold`), colored by model
#' @export
plotDefaultEval <- function(default_eval_tibble, x_default_eval = "train_prop", 
  y_default_eval = "avg_f1_score", xlab = "Train Data Proportion", 
  ylab = "Average F1 Score") {
  .checkArgTibble(default_eval_tibble); .checkArgXDefaultEval(x_default_eval)
  .checkArgYDefaultEval(y_default_eval)
  .checkArgXYLabs(xlab = xlab, ylab = ylab)
  
  if(x_default_eval == "n_fold") {
    default_eval_tibble <- default_eval_tibble |> 
      dplyr::filter(train_prop == 0.8)
  } else {
    default_eval_tibble <- default_eval_tibble |> 
      dplyr::filter(train_prop != 0.8)
  }
  
  default_eval_plot <- ggplot2::ggplot(default_eval_tibble, 
    ggplot2::aes(x = unlist(default_eval_tibble[x_default_eval]), 
      y = unlist(default_eval_tibble[y_default_eval]), color = model)) + 
    ggplot2::geom_line(size = 1.5) + 
    ggplot2::geom_point(size = 3) + 
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 1.5),
      axis.ticks = ggplot2::element_line(linewidth = 1.5, colour = "black"), 
      axis.text = ggplot2::element_text(size = 16, colour = "black"), 
      axis.title = ggplot2::element_text(size = 16, face = "bold"), 
      panel.grid = ggplot2::element_blank(), 
      panel.background = ggplot2::element_blank(), 
      legend.text = ggplot2::element_text(size = 16), 
      legend.title = ggplot2::element_text(size = 16, face = "bold")) +
    ggplot2::labs(x = xlab, y = ylab, color = "Model")
  
  return(default_eval_plot)
}

#' getBaselineComparisonBarplot()
#' 
#' Generates a bar plot that compares model performance with and without 
#' randomly shuffled AMR phenotype labels.
#' 
#' @param non_shuffled_label_results Output of `runMLPipeline()`
#' (`shuffle_labels = FALSE`)
#' @param shuffled_label_results Output of `runMLPipeline()` 
#' (`shuffle_labels = TRUE`)
#' @return A bar plot with balanced accuracy comparisons per antibiotic
#' @export
getBaselineComparisonBarplot <- function(non_shuffled_label_results, 
  shuffled_label_results) {
  .checkArgTibble(non_shuffled_label_results)
  .checkArgTibble(shuffled_label_results)
  
  drugs <- non_shuffled_label_results |> dplyr::select(antibiotic) |> 
    dplyr::pull()
  
  non_shuffled_bal_acc <- non_shuffled_label_results |> 
    dplyr::select(bal_acc) |> 
    dplyr::pull()
  
  shuffled_bal_acc <- shuffled_label_results |> dplyr::select(bal_acc) |> 
    dplyr::pull()
  
  bal_acc_matrix <- matrix(c(non_shuffled_bal_acc, shuffled_bal_acc), 
    nrow = 2, byrow = TRUE)
  
  colnames(bal_acc_matrix) <- drugs
  rownames(bal_acc_matrix) <- c("Non-Shuffled Labels", "Shuffled Labels")
  
  baseline_comparison_barplot <- barplot(bal_acc_matrix, beside = TRUE, 
    legend.text = TRUE, col = c("skyblue", "lightpink"), 
    ylab = "Balanced Accuracy", xlab = "Antibiotic")
  
  return(baseline_comparison_barplot)
}
