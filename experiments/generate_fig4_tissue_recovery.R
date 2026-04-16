#!/usr/bin/env Rscript
# Generate Figure 4: Tissue-specific target recovery
# Panel a: Grouped bar chart (EN vs DSEN shared vs DSEN combined)
# Panel b: Case study heatmap for key drugs

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(patchwork)
  library(scales)
})

source(here("methods", "00_paths.R"))

# ==============================================================================
# Load precomputed results
# ==============================================================================
results_file <- here("results", "dsen_results", "tissue_target_recovery.Rda")
if (!file.exists(results_file)) {
  stop("Run experiments/11_tissue_target_recovery.Rmd first to generate results")
}
load(results_file)

# ==============================================================================
# Panel A: Grouped bar chart of recovery rates
# ==============================================================================

# Reshape summary for plotting
bar_data <- summary_df %>%
  select(threshold, en_any_pct, dsen_shared_any_pct, dsen_combined_any_pct) %>%
  pivot_longer(
    cols = c(en_any_pct, dsen_shared_any_pct, dsen_combined_any_pct),
    names_to = "method",
    values_to = "recovery_pct"
  ) %>%
  mutate(
    method = factor(method,
                    levels = c("en_any_pct", "dsen_shared_any_pct", "dsen_combined_any_pct"),
                    labels = c("Elastic Net", "DSEN (shared)", "DSEN (shared + tissue)")),
    threshold = factor(sprintf("%.0f%%", 100 * threshold))
  )

panel_a <- ggplot(bar_data, aes(x = threshold, y = recovery_pct, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", recovery_pct)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(
    values = c("Elastic Net" = "gray60",
               "DSEN (shared)" = "#6BAED6",
               "DSEN (shared + tissue)" = "#08519C"),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, max(bar_data$recovery_pct) * 1.15),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Bootstrap frequency threshold",
    y = "Drugs with any target recovered (%)",
    tag = "a"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.tag = element_text(face = "bold", size = 14)
  )

# ==============================================================================
# Panel B: Case study heatmap
# ==============================================================================

# Get threshold 0.50 results
recovery_050 <- recovery_df %>%
  filter(n_targets > 0, threshold == 0.50)

# Load targets for annotation
drug_targets_raw <- read_csv(here("data", "drug_targets_curated.csv"), show_col_types = FALSE)
gene_mapping_raw <- read_csv(here("data", "gene_name_mapping.csv"), show_col_types = FALSE)

# Standardize
standardize_gene <- function(gene, mapping_df) {
  match_idx <- match(toupper(gene), toupper(mapping_df$alias))
  if (!is.na(match_idx)) return(mapping_df$hgnc_symbol[match_idx])
  return(gene)
}

targets_std <- drug_targets_raw %>%
  filter(!is.na(target_gene), source != "not_found") %>%
  mutate(target_gene = sapply(target_gene, standardize_gene, mapping_df = gene_mapping_raw))

# Select case study drugs: those with combined recovery + key drugs
case_drugs_combined <- recovery_050 %>%
  filter(dsen_combined_n_found > 0) %>%
  pull(drug)

# Also include key drugs from the prompt if they have targets
key_drugs <- c("PLX-4720", "Dabrafenib", "Nutlin-3a(-)", "Nutlin-3a",
               "NVP-TAE684", "Quizartinib", "Sunitinib",
               "Navitoclax", "Ruxolitinib", "Fedratinib",
               "Afatinib", "Cetuximab")

# Build case study data
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

# Limit to ~12 drugs for readability
if (nrow(case_data) > 12) {
  # Prioritize: drugs with tissue-only recovery, then shared-only, then both
  case_data <- case_data %>%
    mutate(
      priority = case_when(
        dsen_tissue_any & !dsen_shared_any ~ 1,  # tissue only
        dsen_shared_any & dsen_tissue_any ~ 2,     # both
        dsen_shared_any & !dsen_tissue_any ~ 3,    # shared only
        en_any ~ 4,                                  # EN only
        TRUE ~ 5
      )
    ) %>%
    arrange(priority) %>%
    head(12) %>%
    select(-priority)
}

# Reshape for heatmap
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
    values = c("Recovered" = "#08519C", "Not found" = "gray90"),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = NULL,
    tag = "b"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9),
    plot.tag = element_text(face = "bold", size = 14)
  )

# ==============================================================================
# Combine panels
# ==============================================================================
fig <- panel_a + panel_b +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = NULL,
    caption = "Features selected in \u226550% of 200 bootstrap resamples"
  )

# ==============================================================================
# Save
# ==============================================================================
out_dir <- here("doc", "nature_methods", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
latex_dir <- here("doc", "nature_methods", "latex")

ggsave(file.path(out_dir, "fig4_target_recovery.pdf"), fig,
       width = 10, height = 6, device = "pdf")
cat("Saved:", file.path(out_dir, "fig4_target_recovery.pdf"), "\n")

# Copy to latex dir
file.copy(file.path(out_dir, "fig4_target_recovery.pdf"),
          file.path(latex_dir, "fig4_target_recovery.pdf"),
          overwrite = TRUE)
cat("Copied to:", file.path(latex_dir, "fig4_target_recovery.pdf"), "\n")

# Also save to results/figures
ggsave(file.path("results", "figures", "fig4_target_recovery_tissue.pdf"), fig,
       width = 10, height = 6, device = "pdf")
cat("Saved to results/figures/\n")

cat("\nDone.\n")
