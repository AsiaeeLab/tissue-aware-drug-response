#' =============================================================================
#' Drug Response Prediction - Plotting Utilities
#' =============================================================================
#'
#' This script provides functions for creating publication-quality figures.
#'
#' Usage:
#'   library(here)
#'   source(here("methods", "04_plot_utils.R"))
#'
#' =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
  library(scales)
})

#' -----------------------------------------------------------------------------
#' Color Palettes
#' -----------------------------------------------------------------------------

# Define consistent color schemes
COLORS <- list(
  wrong = "#E41A1C",      # Red for wrong/leakage

  right = "#377EB8",      # Blue for right/correct
  dsen = "#FF7F00",       # Orange for DSEN
  en = "#4DAF4A",         # Green for elastic net
  neutral = "#666666",    # Gray for neutral
  highlight = "#984EA3"   # Purple for highlights
)

# Color palette for tissues
get_tissue_colors <- function(n_tissues) {
  if (n_tissues <= 8) {
    return(brewer.pal(max(3, n_tissues), "Set2"))
  } else if (n_tissues <= 12) {
    return(brewer.pal(n_tissues, "Set3"))
  } else {
    return(colorRampPalette(brewer.pal(12, "Set3"))(n_tissues))
  }
}

#' -----------------------------------------------------------------------------
#' Publication Theme
#' -----------------------------------------------------------------------------

#' ggplot2 theme for publication-quality figures
theme_publication <- function(base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Panel
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),

      # Axis
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),

      # Legend
      legend.background = element_blank(),
      legend.key = element_blank(),

      # Strip (for facets)
      strip.background = element_rect(fill = "gray95", color = "black"),
      strip.text = element_text(face = "bold"),

      # Plot
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

#' -----------------------------------------------------------------------------
#' CV Curve Plots
#' -----------------------------------------------------------------------------

#' Plot CV error curve with standard error bars
#'
#' @param cv_result CV result object with lambda, cvm, cvse
#' @param title Plot title
#' @param show_lambda_min Whether to mark lambda.min
#' @return ggplot object
plot_cv_curve <- function(cv_result, title = "CV Error", show_lambda_min = TRUE) {

  df <- data.frame(
    lambda = cv_result$lambda,
    cvm = cv_result$cvm,
    cvse = cv_result$cvse,
    upper = cv_result$cvm + cv_result$cvse,
    lower = cv_result$cvm - cv_result$cvse
  )

  p <- ggplot(df, aes(x = log(lambda), y = cvm)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = COLORS$right) +
    geom_line(color = COLORS$right, linewidth = 1) +
    geom_point(color = COLORS$right, size = 1) +
    labs(x = expression(log(lambda)), y = "Cross-Validation Error (MSE)",
         title = title) +
    theme_publication()

  if (show_lambda_min && !is.null(cv_result$lambda.min)) {
    lambda_min_idx <- which.min(cv_result$cvm)
    p <- p +
      geom_vline(xintercept = log(cv_result$lambda.min),
                 linetype = "dashed", color = COLORS$highlight) +
      annotate("text", x = log(cv_result$lambda.min), y = max(df$cvm),
               label = expression(lambda[min]), hjust = -0.2, vjust = 1,
               color = COLORS$highlight)
  }

  return(p)
}

#' Plot comparison of wrong vs right CV curves
#'
#' @param wrong_result Wrong CV result
#' @param right_result Right CV result
#' @param drug Drug name for title
#' @return ggplot object
plot_cv_comparison <- function(wrong_result, right_result, drug = "") {

  # Prepare data for wrong
  df_wrong <- data.frame(
    lambda = wrong_result$lambda,
    cvm = wrong_result$cvm,
    cvse = wrong_result$cvse,
    method = "Wrong CV (data leakage)"
  )

  # Prepare data for right
  df_right <- data.frame(
    lambda = right_result$lambda,
    cvm = right_result$cvm,
    cvse = right_result$cvse,
    method = "Right CV (per-fold selection)"
  )

  df <- rbind(df_wrong, df_right)

  p <- ggplot(df, aes(x = log(lambda), y = cvm, color = method, fill = method)) +
    geom_ribbon(aes(ymin = cvm - cvse, ymax = cvm + cvse), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Wrong CV (data leakage)" = COLORS$wrong,
                                   "Right CV (per-fold selection)" = COLORS$right)) +
    scale_fill_manual(values = c("Wrong CV (data leakage)" = COLORS$wrong,
                                  "Right CV (per-fold selection)" = COLORS$right)) +
    labs(x = expression(log(lambda)), y = "Cross-Validation Error (MSE)",
         title = ifelse(drug != "", paste("CV Comparison:", drug), "CV Comparison"),
         color = "", fill = "") +
    theme_publication() +
    theme(legend.position = "bottom")

  return(p)
}

#' -----------------------------------------------------------------------------
#' MSE Comparison Plots
#' -----------------------------------------------------------------------------

#' Scatter plot comparing MSE from two methods
#'
#' @param df Data frame with mse_wrong and mse_right columns
#' @param title Plot title
#' @return ggplot object
plot_mse_scatter <- function(df, title = "MSE Comparison: Wrong vs Right CV") {

  # Determine which method is better for each drug
  df$better <- ifelse(df$mse_wrong < df$mse_right, "Wrong", "Right")

  p <- ggplot(df, aes(x = mse_wrong, y = mse_right, color = better)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Wrong" = COLORS$wrong, "Right" = COLORS$right),
                       labels = c("Wrong" = "Wrong CV better (leakage)",
                                   "Right" = "Right CV better")) +
    labs(x = "MSE (Wrong CV)", y = "MSE (Right CV)",
         title = title,
         subtitle = sprintf("n = %d drugs", nrow(df)),
         color = "") +
    coord_fixed() +
    theme_publication() +
    theme(legend.position = "bottom")

  return(p)
}

#' Histogram of MSE percent increase
#'
#' @param pct_increase Vector of percent increases
#' @param title Plot title
#' @return ggplot object
plot_inflation_histogram <- function(pct_increase, title = "MSE Inflation Distribution") {

  df <- data.frame(pct = pct_increase)
  mean_pct <- mean(pct_increase, na.rm = TRUE)
  median_pct <- median(pct_increase, na.rm = TRUE)

  p <- ggplot(df, aes(x = pct)) +
    geom_histogram(bins = 30, fill = COLORS$right, color = "white", alpha = 0.8) +
    geom_vline(xintercept = mean_pct, color = COLORS$wrong,
               linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = median_pct, color = COLORS$highlight,
               linetype = "dotted", linewidth = 1) +
    annotate("text", x = mean_pct, y = Inf,
             label = sprintf("Mean: %.1f%%", mean_pct),
             hjust = -0.1, vjust = 2, color = COLORS$wrong) +
    annotate("text", x = median_pct, y = Inf,
             label = sprintf("Median: %.1f%%", median_pct),
             hjust = -0.1, vjust = 4, color = COLORS$highlight) +
    labs(x = "MSE Increase (%)", y = "Number of Drugs",
         title = title) +
    theme_publication()

  return(p)
}

#' -----------------------------------------------------------------------------
#' Feature Stability Plots
#' -----------------------------------------------------------------------------

#' Histogram of Jaccard similarity
#'
#' @param jaccard Vector of Jaccard values
#' @param threshold Bootstrap threshold used
#' @param title Plot title
#' @return ggplot object
plot_jaccard_histogram <- function(jaccard, threshold = 0.8,
                                    title = "Feature Set Overlap") {

  df <- data.frame(jaccard = jaccard)
  zero_pct <- mean(jaccard == 0, na.rm = TRUE) * 100

  p <- ggplot(df, aes(x = jaccard)) +
    geom_histogram(bins = 20, fill = COLORS$right, color = "white", alpha = 0.8) +
    geom_vline(xintercept = mean(jaccard, na.rm = TRUE),
               color = COLORS$wrong, linetype = "dashed", linewidth = 1) +
    annotate("text", x = 0.5, y = Inf,
             label = sprintf("%.0f%% with zero overlap", zero_pct),
             hjust = 0.5, vjust = 2, color = COLORS$wrong, size = 4) +
    labs(x = sprintf("Jaccard Similarity (threshold = %.1f)", threshold),
         y = "Number of Drugs",
         title = title) +
    scale_x_continuous(limits = c(0, 1)) +
    theme_publication()

  return(p)
}

#' Bar plot comparing feature counts
#'
#' @param n_wrong Number of features from wrong CV
#' @param n_right Number of features from right CV
#' @param title Plot title
#' @return ggplot object
plot_feature_counts <- function(n_wrong, n_right,
                                 title = "Feature Count Comparison") {

  df <- data.frame(
    method = rep(c("Wrong CV", "Right CV"), each = length(n_wrong)),
    count = c(n_wrong, n_right)
  )

  df_summary <- data.frame(
    method = c("Wrong CV", "Right CV"),
    mean_count = c(mean(n_wrong, na.rm = TRUE), mean(n_right, na.rm = TRUE)),
    se_count = c(sd(n_wrong, na.rm = TRUE) / sqrt(sum(!is.na(n_wrong))),
                  sd(n_right, na.rm = TRUE) / sqrt(sum(!is.na(n_right))))
  )

  p <- ggplot(df_summary, aes(x = method, y = mean_count, fill = method)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count),
                  width = 0.2) +
    geom_text(aes(label = sprintf("%.1f", mean_count)),
              vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("Wrong CV" = COLORS$wrong,
                                  "Right CV" = COLORS$right)) +
    labs(x = "", y = "Mean Number of Selected Features",
         title = title) +
    theme_publication() +
    theme(legend.position = "none")

  return(p)
}

#' -----------------------------------------------------------------------------
#' DSEN Plots
#' -----------------------------------------------------------------------------

#' Heatmap of DSEN coefficients by tissue
#'
#' @param coef_matrix Matrix of coefficients (features x tissues)
#' @param drug Drug name for title
#' @param top_n Number of top features to show
#' @return ggplot object
plot_dsen_heatmap <- function(coef_matrix, drug = "", top_n = 15) {

  # Select top features by overall magnitude
  row_maxabs <- apply(abs(coef_matrix), 1, max, na.rm = TRUE)
  top_idx <- order(row_maxabs, decreasing = TRUE)[1:min(top_n, nrow(coef_matrix))]

  mat <- coef_matrix[top_idx, , drop = FALSE]

  # Convert to long format
  df <- expand.grid(
    Feature = rownames(mat),
    Tissue = colnames(mat)
  )
  df$Coefficient <- as.vector(mat)
  df$Feature <- factor(df$Feature, levels = rev(rownames(mat)))

  # Set color limits symmetrically
  max_abs <- max(abs(df$Coefficient), na.rm = TRUE)

  p <- ggplot(df, aes(x = Tissue, y = Feature, fill = Coefficient)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = COLORS$right, mid = "white", high = COLORS$wrong,
                         midpoint = 0, limits = c(-max_abs, max_abs)) +
    labs(x = "", y = "",
         title = ifelse(drug != "", paste("DSEN Coefficients:", drug),
                        "DSEN Coefficients"),
         fill = "Coefficient") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))

  return(p)
}

#' DSEN vs EN comparison scatter plot
#'
#' @param rmse_dsen RMSE values for DSEN
#' @param rmse_en RMSE values for EN
#' @param title Plot title
#' @return ggplot object
plot_dsen_vs_en <- function(rmse_dsen, rmse_en, title = "DSEN vs Elastic Net") {

  df <- data.frame(
    dsen = rmse_dsen,
    en = rmse_en,
    better = ifelse(rmse_dsen < rmse_en, "DSEN", "EN")
  )

  n_dsen_better <- sum(df$better == "DSEN")
  n_en_better <- sum(df$better == "EN")

  p <- ggplot(df, aes(x = dsen, y = en, color = better)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("DSEN" = COLORS$dsen, "EN" = COLORS$en)) +
    labs(x = "RMSE (DSEN)", y = "RMSE (Elastic Net)",
         title = title,
         subtitle = sprintf("DSEN better: %d | EN better: %d",
                           n_dsen_better, n_en_better),
         color = "Better method") +
    coord_fixed() +
    theme_publication() +
    theme(legend.position = "bottom")

  return(p)
}

#' -----------------------------------------------------------------------------
#' Cross-Tissue Plots
#' -----------------------------------------------------------------------------
#' Heatmap for cross-tissue prediction performance
#'
#' @param perf_matrix Matrix of performance (train tissue x test tissue)
#' @param metric Name of metric for label
#' @return ggplot object
plot_cross_tissue_heatmap <- function(perf_matrix, metric = "RMSE") {

  df <- expand.grid(
    Train = rownames(perf_matrix),
    Test = colnames(perf_matrix)
  )
  df$Performance <- as.vector(perf_matrix)

  p <- ggplot(df, aes(x = Test, y = Train, fill = Performance)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "plasma") +
    labs(x = "Test Tissue", y = "Training Tissue",
         title = sprintf("Cross-Tissue Prediction (%s)", metric),
         fill = metric) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}

#' -----------------------------------------------------------------------------
#' Save Functions
#' -----------------------------------------------------------------------------

#' Save plot to PDF and PNG
#'
#' @param plot ggplot object
#' @param filename Base filename (without extension)
#' @param dir Output directory
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi DPI for PNG
save_publication_plot <- function(plot, filename, dir = results_dirs$figures,
                                   width = 8, height = 6, dpi = 300) {

  # Save PDF
  pdf_path <- file.path(dir, paste0(filename, ".pdf"))
  ggsave(pdf_path, plot, width = width, height = height)

  # Save PNG
  png_path <- file.path(dir, paste0(filename, ".png"))
  ggsave(png_path, plot, width = width, height = height, dpi = dpi)

  message(sprintf("Saved: %s (.pdf and .png)", filename))
  invisible(list(pdf = pdf_path, png = png_path))
}

message("Plotting utilities loaded.")
