#!/usr/bin/env Rscript
# ==============================================================================
# Generate publication figures from threshold_25_v2 DSEN results
#
# Produces:
#   Fig 2: DSEN performance (waterfall, difficulty scatter, lamrat distribution)
#   Supp Fig 1: Feature sparsity (EN vs DSEN)
#   Supp Fig 2: Shared vs tissue-specific features
#   Supp Fig 3: Lambda ratio distribution (expanded)
#   Supp Fig 4: Target recovery (bar chart + case study heatmap)
#
# Fig 1 (leakage) and Table 1 (audit) are NOT regenerated here.
#
# Prerequisites:
#   - Run methods/run_dsen_production.R with MIN_SAMPLES_PER_TISSUE=25
#   - Run methods/run_dsen_bootstrap.R
#   - Run experiments/11_tissue_target_recovery.Rmd
#
# Output directories:
#   doc/nature_methods/figures/
#   doc/nature_methods/latex/
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(patchwork)
  library(scales)
})

source(here("methods", "00_paths.R"))
source(here("experiments", "figure_style_nature_methods.R"))

# ==============================================================================
# Output directories
# ==============================================================================
fig_dir <- here("doc", "nature_methods", "figures")
latex_dir <- here("doc", "nature_methods", "latex")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(latex_dir, showWarnings = FALSE, recursive = TRUE)

colors <- nm_colors

theme_pub <- function(base_size = 8) {
  nm_theme_compact(base_size = base_size)
}

# ==============================================================================
# Load data
# ==============================================================================
cat("Loading v2 results...\n")

dsen_df <- read_csv(here("results", "tables", "threshold_25_v2", "dsen_full_results.csv"),
                    show_col_types = FALSE) %>%
  filter(status == "success")
cat("  DSEN full results:", nrow(dsen_df), "drugs\n")

quartile_df <- read_csv(here("results", "tables", "threshold_25_v2", "difficulty_quartile.csv"),
                        show_col_types = FALSE)

cat("Loading v2 bootstrap results...\n")
load(here("results", "dsen_results", "threshold_25_v2", "dsen_bootstrap_results.Rda"))
bootstrap_results <- all_results
rm(all_results)
cat("  Bootstrap results:", length(bootstrap_results), "drugs\n")

# Target recovery data
load(here("results", "dsen_results", "threshold_25_v2", "tissue_target_recovery.Rda"))
recovery_df <- tissue_target_recovery$per_drug
summary_recovery <- tissue_target_recovery$summary
cat("  Target recovery data loaded\n")

recovery_summary_csv <- read_csv(
  here("results", "tables", "threshold_25_v2", "target_recovery_extended_summary.csv"),
  show_col_types = FALSE
)

# ==============================================================================
# Figure 2: DSEN Performance (3-panel) -- 183mm width
# ==============================================================================
cat("\nGenerating Figure 2: DSEN Performance...\n")

plot_df <- dsen_df %>%
  arrange(pct_improvement) %>%
  mutate(
    rank = row_number(),
    direction = ifelse(pct_improvement > 0, "DSEN Better", "EN Better")
  )

n_total <- nrow(dsen_df)
n_dsen_better <- sum(dsen_df$pct_improvement > 0)
n_en_better <- n_total - n_dsen_better
pct_dsen_wins <- 100 * n_dsen_better / n_total
mean_improvement <- mean(dsen_df$pct_improvement)
median_improvement <- median(dsen_df$pct_improvement)

# Panel A: Waterfall (no text annotations; color-only)
fig2a <- ggplot(plot_df, aes(x = rank, y = pct_improvement, fill = direction)) +
  geom_col(width = 1, color = NA) +
  scale_fill_manual(
    values = c("DSEN Better" = colors$win, "EN Better" = colors$loss),
    name = NULL
  ) +
  geom_hline(yintercept = 0, color = colors$dark, linewidth = 0.5) +
  geom_hline(yintercept = mean_improvement,
             linetype = "dashed", color = colors$neutral, linewidth = 0.5) +
  labs(
    x = "Drug (ranked by improvement)",
    y = "MSE improvement (%)",
    tag = "a"
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  theme_pub() +
  theme(legend.position = "none")

# Panel B: Difficulty scatter (no text annotations on quartile bands)
quartile_bounds <- dsen_df %>%
  mutate(quartile = ntile(mse_en, 4)) %>%
  group_by(quartile) %>%
  summarize(xmin = min(mse_en), xmax = max(mse_en), .groups = "drop")

fig2b <- dsen_df %>%
  mutate(quartile = ntile(mse_en, 4), wins = pct_improvement > 0) %>%
  ggplot(aes(x = mse_en, y = pct_improvement)) +
  geom_rect(
    data = quartile_bounds,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = colors$light_gray,
    alpha = 0.55
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = colors$neutral, linewidth = 0.45) +
  geom_point(aes(color = wins), alpha = 0.72, size = 1.55) +
  scale_color_manual(values = c("TRUE" = colors$win, "FALSE" = colors$loss), guide = "none") +
  scale_x_continuous(
    trans = "sqrt",
    breaks = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1),
    labels = number_format(accuracy = 0.001)
  ) +
  labs(
    x = "Baseline EN MSE (sqrt scale)",
    y = "DSEN improvement (%)",
    tag = "b"
  ) +
  theme_pub()

# Panel C (NEW): Lambda ratio (rho) distribution -- colored by shared- vs tissue-favoring
lamrat_counts <- dsen_df %>%
  count(lamrat_dsen) %>%
  mutate(
    pct = n / sum(n) * 100,
    favors = ifelse(lamrat_dsen < 1, "Shared-favoring", "Tissue-favoring")
  )

fig2c <- ggplot(lamrat_counts, aes(x = factor(lamrat_dsen), y = n, fill = favors)) +
  geom_col(width = 0.7, alpha = 0.90) +
  scale_fill_manual(
    values = c("Shared-favoring" = colors$dsen_shared, "Tissue-favoring" = colors$dsen_tissue),
    name = NULL
  ) +
  labs(
    x = expression("Penalty multiplier " * rho),
    y = "Drugs",
    tag = "c"
  ) +
  theme_pub() +
  theme(
    legend.position = c(0.30, 0.88),
    legend.background = element_rect(fill = alpha("white", 0.90), color = NA),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    panel.grid.major.x = element_blank()
  )

fig2 <- fig2a + fig2b + fig2c +
  plot_layout(widths = c(1.25, 1, 0.85))

nm_save(fig2, file.path(fig_dir, "fig2_dsen_overview.pdf"), width_mm = 183, height_mm = 65)
file.copy(file.path(fig_dir, "fig2_dsen_overview.pdf"),
          file.path(latex_dir, "fig2_dsen_overview.pdf"), overwrite = TRUE)
cat("  Saved fig2_dsen_overview.pdf\n")

# ==============================================================================
# Supplementary Fig 1: Feature Sparsity (EN vs DSEN) -- 183mm
# ==============================================================================
cat("\nGenerating Supp Fig 1: Feature Sparsity...\n")

THRESH <- 0.50
feature_df <- map_dfr(bootstrap_results, function(r) {
  if (r$status != "success") return(NULL)
  nboot <- r$en_bootstrap$nboot
  en_stable <- sum(r$en_bootstrap$cnt / nboot >= THRESH)
  dsen_shared_stable <- sum(r$dsen_bootstrap$cnt / r$dsen_bootstrap$nboot >= THRESH)
  tibble(
    drug = r$drug,
    n_en = en_stable,
    n_dsen_shared = dsen_shared_stable
  )
})

# Panel A: Scatter of EN vs DSEN shared feature counts
figS1a <- ggplot(feature_df, aes(x = n_en, y = n_dsen_shared)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = colors$neutral) +
  geom_point(alpha = 0.5, size = 1.2, color = colors$en) +
  labs(
    x = "Stable features (EN)",
    y = "Stable features (DSEN shared)",
    tag = "a"
  ) +
  theme_pub()

# Panel B: Distribution of EN feature counts (no text annotations)
figS1b <- ggplot(feature_df, aes(x = n_en)) +
  geom_histogram(bins = 30, fill = colors$en, color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(feature_df$n_en), color = colors$dark, linewidth = 0.5) +
  labs(x = "Stable EN features", y = "Drugs", tag = "b") +
  theme_pub()

# Panel C: Distribution of DSEN shared feature counts (no text annotations)
figS1c <- ggplot(feature_df, aes(x = n_dsen_shared)) +
  geom_histogram(bins = 30, fill = colors$dsen_shared, color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(feature_df$n_dsen_shared), color = colors$dark, linewidth = 0.5) +
  labs(x = "Stable DSEN shared features", y = "Drugs", tag = "c") +
  theme_pub()

figS1 <- figS1a + figS1b + figS1c +
  plot_layout(widths = c(1, 1, 1))

ggsave(file.path(fig_dir, "figS1_feature_sparsity.pdf"), figS1,
       width = 183/25.4, height = 65/25.4, device = cairo_pdf)
file.copy(file.path(fig_dir, "figS1_feature_sparsity.pdf"),
          file.path(latex_dir, "figS1_feature_sparsity.pdf"), overwrite = TRUE)
cat("  Saved figS1_feature_sparsity.pdf\n")

# ==============================================================================
# Supplementary Fig 2: Shared vs Tissue-Specific Features -- 183mm
# ==============================================================================
cat("\nGenerating Supp Fig 2: Shared vs Tissue-Specific Features...\n")

coef_df <- map_dfr(bootstrap_results, function(r) {
  if (r$status != "success") return(NULL)
  nboot <- r$dsen_bootstrap$nboot

  shared_freq <- r$dsen_bootstrap$cnt / nboot
  n_shared <- sum(shared_freq >= THRESH)

  tissue_freq <- r$dsen_bootstrap$tissue_cnt / nboot
  n_tissue <- sum(apply(tissue_freq >= THRESH, 1, any))

  tibble(drug = r$drug, n_shared = n_shared, n_tissue = n_tissue)
})

n_no_shared <- sum(coef_df$n_shared == 0)

figS2a <- ggplot(coef_df, aes(x = n_shared, y = n_tissue)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = colors$neutral) +
  geom_point(alpha = 0.5, size = 1.2, color = colors$primary_orange) +
  labs(
    x = "Shared stable features",
    y = "Tissue-specific stable features",
    tag = "a"
  ) +
  theme_pub()

figS2b <- ggplot(coef_df, aes(x = n_shared)) +
  geom_histogram(bins = 30, fill = colors$dsen_shared, color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(coef_df$n_shared), color = colors$dark, linewidth = 0.5) +
  labs(
    x = "Number of shared stable features",
    y = "Drugs",
    tag = "b"
  ) +
  theme_pub()

figS2 <- figS2a + figS2b +
  plot_layout(widths = c(1, 1))

ggsave(file.path(fig_dir, "figS2_shared_vs_tissue.pdf"), figS2,
       width = 183/25.4, height = 75/25.4, device = cairo_pdf)
file.copy(file.path(fig_dir, "figS2_shared_vs_tissue.pdf"),
          file.path(latex_dir, "figS2_shared_vs_tissue.pdf"), overwrite = TRUE)
cat("  Saved figS2_shared_vs_tissue.pdf\n")

# ==============================================================================
# Supplementary Fig 3: Lambda Ratio Distribution (expanded) -- 183mm
# ==============================================================================
cat("\nGenerating Supp Fig 3: Lambda Ratio Distribution...\n")

figS3 <- ggplot(lamrat_counts, aes(x = factor(lamrat_dsen), y = n, fill = favors)) +
  geom_col(width = 0.7, alpha = 0.85) +
  scale_fill_manual(
    values = c("Shared-favoring" = colors$dsen_shared, "Tissue-favoring" = colors$dsen_tissue),
    name = NULL
  ) +
  labs(
    x = expression("Shared-block penalty multiplier " * rho),
    y = "Number of drugs"
  ) +
  theme_pub() +
  theme(
    legend.position = c(0.35, 0.90),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "figS3_lambda_ratio_distribution.pdf"), figS3,
       width = 183/25.4, height = 100/25.4, device = cairo_pdf)
file.copy(file.path(fig_dir, "figS3_lambda_ratio_distribution.pdf"),
          file.path(latex_dir, "figS3_lambda_ratio_distribution.pdf"), overwrite = TRUE)
cat("  Saved figS3_lambda_ratio_distribution.pdf\n")

# ==============================================================================
# Supplementary Fig 4: Target Recovery -- 183mm
# ==============================================================================
cat("\nGenerating Supp Fig 4: Target Recovery...\n")

# Panel A: Grouped bar chart of recovery rates across thresholds
bar_data <- recovery_summary_csv %>%
  select(threshold, en_any_pct, dsen_shared_any_pct, dsen_combined_any_pct) %>%
  pivot_longer(
    cols = c(en_any_pct, dsen_shared_any_pct, dsen_combined_any_pct),
    names_to = "method",
    values_to = "recovery_pct"
  ) %>%
  mutate(
    method = factor(method,
                    levels = c("en_any_pct", "dsen_shared_any_pct", "dsen_combined_any_pct"),
                    labels = c("EN", "DSEN shared", "DSEN combined")),
    threshold = factor(sprintf("%.0f%%", 100 * threshold))
  )

panel_a <- ggplot(bar_data, aes(x = threshold, y = recovery_pct, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
  scale_fill_manual(
    values = c("EN" = colors$en,
               "DSEN shared" = colors$dsen_shared,
               "DSEN combined" = colors$dsen_combined),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, max(bar_data$recovery_pct) * 1.12),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Bootstrap frequency threshold",
    y = "Drugs with any target recovered (%)",
    tag = "a"
  ) +
  theme_pub() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    panel.grid.major.x = element_blank()
  )

# Panel B: Case study heatmap at threshold 0.50
drug_targets_raw <- read_csv(here("data", "drug_targets_curated.csv"), show_col_types = FALSE)
gene_mapping_raw <- read_csv(here("data", "gene_name_mapping.csv"), show_col_types = FALSE)

standardize_gene <- function(gene, mapping_df) {
  match_idx <- match(toupper(gene), toupper(mapping_df$alias))
  if (!is.na(match_idx)) return(mapping_df$hgnc_symbol[match_idx])
  return(gene)
}

targets_std <- drug_targets_raw %>%
  filter(!is.na(target_gene), source != "not_found") %>%
  mutate(target_gene = sapply(target_gene, standardize_gene, mapping_df = gene_mapping_raw))

recovery_050 <- recovery_df %>%
  filter(n_targets > 0, threshold == 0.50)

case_drugs_combined <- recovery_050 %>%
  filter(dsen_combined_any == TRUE) %>%
  pull(drug)

key_drugs <- c("PLX-4720", "Dabrafenib", "Nutlin-3a(-)", "Nutlin-3a",
               "NVP-TAE684", "Quizartinib", "Sunitinib",
               "Navitoclax", "Ruxolitinib", "Fedratinib",
               "Afatinib", "Cetuximab", "Gefitinib")

case_data <- recovery_050 %>%
  filter(drug %in% union(case_drugs_combined, key_drugs)) %>%
  left_join(
    targets_std %>%
      group_by(drug_name) %>%
      summarize(primary_target = first(target_gene), .groups = "drop"),
    by = c("drug" = "drug_name")
  ) %>%
  filter(!is.na(primary_target)) %>%
  mutate(drug_label = paste0(drug, " (", primary_target, ")")) %>%
  select(drug, drug_label, primary_target,
         en_any, dsen_shared_any, dsen_tissue_any) %>%
  distinct()

if (nrow(case_data) > 12) {
  case_data <- case_data %>%
    mutate(
      priority = case_when(
        dsen_tissue_any & !dsen_shared_any ~ 1,
        dsen_shared_any & dsen_tissue_any ~ 2,
        dsen_shared_any & !dsen_tissue_any ~ 3,
        en_any ~ 4,
        TRUE ~ 5
      )
    ) %>%
    arrange(priority) %>%
    head(12) %>%
    select(-priority)
}

case_long <- case_data %>%
  pivot_longer(
    cols = c(en_any, dsen_shared_any, dsen_tissue_any),
    names_to = "method",
    values_to = "recovered"
  ) %>%
  mutate(
    method = factor(method,
                    levels = c("en_any", "dsen_shared_any", "dsen_tissue_any"),
                    labels = c("EN", "DSEN\nShared", "DSEN\nTissue")),
    drug_label = fct_reorder(drug_label, as.numeric(recovered), .fun = sum),
    status = ifelse(recovered, "Recovered", "Not found")
  )

panel_b <- ggplot(case_long, aes(x = method, y = drug_label, fill = status)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_manual(
    values = c("Recovered" = colors$dsen_combined, "Not found" = colors$not_found),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = NULL,
    tag = "b"
  ) +
  theme_pub() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 6),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6)
  )

figS4 <- panel_a + panel_b +
  plot_layout(widths = c(1.2, 1))

ggsave(file.path(fig_dir, "fig4_target_recovery.pdf"), figS4,
       width = 183/25.4, height = 100/25.4, device = cairo_pdf)
file.copy(file.path(fig_dir, "fig4_target_recovery.pdf"),
          file.path(latex_dir, "fig4_target_recovery.pdf"), overwrite = TRUE)
cat("  Saved fig4_target_recovery.pdf\n")

cat("\n=== All figures generated ===\n")
cat("Output:", fig_dir, "\n")
cat("Copied to:", latex_dir, "\n")
