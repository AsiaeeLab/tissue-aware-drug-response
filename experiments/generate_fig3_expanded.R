#!/usr/bin/env Rscript
# ==============================================================================
# Figure 3: Target recovery + feature structure (3-panel)
#
# Panel a: Recovery curves (line plot across thresholds)
# Panel b: Shared vs tissue scatter
# Panel c (NEW): Feature modality breakdown (stacked bar)
#
# Output: doc/nature_methods/figures/fig3_dsen_characterization.pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(patchwork)
  library(purrr)
  library(readr)
  library(scales)
  library(tidyr)
})

source(here("methods", "00_paths.R"))
source(here("experiments", "figure_style_nature_methods.R"))

# ==============================================================================
# Load data
# ==============================================================================

recovery_summary <- read_csv(
  here("results", "tables", "threshold_25_v2", "target_recovery_extended_summary.csv"),
  show_col_types = FALSE
)

feature_summary <- read_csv(
  here("results", "tables", "threshold_25_v2", "feature_stability_summary.csv"),
  show_col_types = FALSE
)

# Load bootstrap results for modality breakdown
cat("Loading bootstrap results for modality analysis...\n")
load(here("results", "dsen_results", "threshold_25_v2", "dsen_bootstrap_results.Rda"))
bootstrap_results <- all_results
rm(all_results)
cat("  Bootstrap results:", length(bootstrap_results), "drugs\n")

# ==============================================================================
# Helper: classify feature modality from name
# ==============================================================================
rppa_mapping <- read_csv(here("data", "rppa_gene_mapping.csv"), show_col_types = FALSE)

get_feature_modality <- function(feature_name) {
  if (grepl("_mut\\.", feature_name)) return("Mutation")
  if (grepl("_Expr$", feature_name)) return("Expression")
  if (grepl("_CN$", feature_name)) return("CN")
  if (feature_name %in% rppa_mapping$rppa_feature) return("RPPA")
  "Other"
}

# ==============================================================================
# Panel a: Recovery curves (direct line-end labels, no legend box)
# ==============================================================================

recovery_long <- recovery_summary %>%
  transmute(
    threshold,
    EN = en_any_pct,
    `DSEN shared` = dsen_shared_any_pct,
    `DSEN tissue` = dsen_tissue_any_pct,
    `DSEN combined` = dsen_combined_any_pct
  ) %>%
  pivot_longer(cols = -threshold, names_to = "method", values_to = "recovery_pct") %>%
  mutate(
    method = factor(
      method,
      levels = c("DSEN combined", "DSEN tissue", "DSEN shared", "EN")
    )
  )

# Get rightmost points for direct labels
label_data <- recovery_long %>%
  group_by(method) %>%
  filter(threshold == max(threshold)) %>%
  ungroup()

panel_a <- ggplot(
  recovery_long,
  aes(x = threshold, y = recovery_pct, color = method, group = method)
) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = method), size = 2.0) +
  # Direct line-end labels instead of legend
  geom_text(
    data = label_data,
    aes(x = threshold + 0.04, label = method),
    hjust = 0,
    size = 2.2,
    family = nm_base_family,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    trans = "sqrt",
    breaks = c(0.005, 0.05, 0.1, 0.25, 0.5, 0.75),
    labels = c("0.005", "0.05", "0.1", "0.25", "0.5", "0.75"),
    expand = expansion(mult = c(0.02, 0.30))
  ) +
  scale_color_manual(
    values = c(
      "EN" = nm_colors$en,
      "DSEN shared" = nm_colors$dsen_shared,
      "DSEN tissue" = nm_colors$dsen_tissue,
      "DSEN combined" = nm_colors$dsen_combined
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("EN" = 16, "DSEN shared" = 17, "DSEN tissue" = 15, "DSEN combined" = 18),
    guide = "none"
  ) +
  labs(
    x = "Bootstrap frequency threshold (sqrt scale)",
    y = "Target recovery (%)",
    tag = "a"
  ) +
  nm_theme_compact()

# ==============================================================================
# Panel b: Shared vs tissue scatter (no text annotations)
# ==============================================================================

zero_shared_pct <- 100 * mean(feature_summary$dsen_shared_zero, na.rm = TRUE)
sparser_pct <- 100 * mean(feature_summary$dsen_shared_sparser_than_en, na.rm = TRUE)

feature_plot_df <- feature_summary %>%
  mutate(shared_zero_group = if_else(dsen_shared_zero, "Zero shared", "Non-zero shared"))

panel_b <- ggplot(
  feature_plot_df,
  aes(x = dsen_shared_stable_features, y = dsen_tissue_stable_features, color = shared_zero_group)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = nm_colors$neutral, linewidth = 0.45) +
  geom_point(alpha = 0.80, size = 1.65) +
  scale_color_manual(
    values = c("Non-zero shared" = nm_colors$primary_orange, "Zero shared" = nm_colors$leakage),
    guide = "none"
  ) +
  labs(
    x = "Stable shared features",
    y = "Stable tissue-specific features",
    tag = "b"
  ) +
  nm_theme_compact()

# ==============================================================================
# Panel c (NEW): Feature modality breakdown (stacked bar)
# ==============================================================================

THRESH <- 0.50

modality_data <- map_dfr(bootstrap_results, function(r) {
  if (r$status != "success") return(NULL)

  nboot_en <- r$en_bootstrap$nboot
  nboot_dsen <- r$dsen_bootstrap$nboot

  # EN stable features
  en_stable <- names(r$en_bootstrap$cnt)[r$en_bootstrap$cnt / nboot_en >= THRESH]

  # DSEN shared stable features
  shared_stable <- names(r$dsen_bootstrap$cnt)[r$dsen_bootstrap$cnt / nboot_dsen >= THRESH]

  # DSEN tissue stable features (any tissue)
  tissue_stable <- rownames(r$dsen_bootstrap$tissue_cnt)[
    apply(r$dsen_bootstrap$tissue_cnt / nboot_dsen >= THRESH, 1, any)
  ]

  bind_rows(
    if (length(en_stable) > 0) {
      tibble(method = "EN", feature = en_stable)
    },
    if (length(shared_stable) > 0) {
      tibble(method = "DSEN shared", feature = shared_stable)
    },
    if (length(tissue_stable) > 0) {
      tibble(method = "DSEN tissue", feature = tissue_stable)
    }
  ) %>%
    mutate(drug = r$drug)
})

# Classify modality for each feature
modality_data <- modality_data %>%
  mutate(modality = vapply(feature, get_feature_modality, character(1)))

# Aggregate: count features per method per modality
modality_agg <- modality_data %>%
  count(method, modality) %>%
  group_by(method) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    method = factor(method, levels = c("EN", "DSEN shared", "DSEN tissue")),
    modality = factor(modality, levels = c("Expression", "Mutation", "CN", "RPPA", "Other"))
  )

panel_c <- ggplot(modality_agg, aes(x = method, y = pct, fill = modality)) +
  geom_col(width = 0.65, alpha = 0.90) +
  scale_fill_manual(values = nm_modality_colors, name = "Feature type") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = NULL,
    y = "Feature composition (%)",
    tag = "c"
  ) +
  nm_theme_compact() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.30, "cm"),
    legend.text = element_text(size = 6),
    panel.grid.major.x = element_blank()
  )

# ==============================================================================
# Combine: top row (a, b), bottom row (c spanning full width) or 1x3
# ==============================================================================

fig3 <- (panel_a + panel_b + panel_c) +
  plot_layout(widths = c(1.15, 1, 0.85))

out_path <- here("doc", "nature_methods", "figures", "fig3_dsen_characterization.pdf")
nm_save(fig3, out_path, width_mm = 183, height_mm = 84)
file.copy(
  out_path,
  here("doc", "nature_methods", "latex", "fig3_dsen_characterization.pdf"),
  overwrite = TRUE
)

cat("Figure 3 saved:", out_path, "\n")
cat(sprintf("  Zero shared features: %.1f%%\n", zero_shared_pct))
cat(sprintf("  Shared block sparser than EN: %.1f%%\n", sparser_pct))
