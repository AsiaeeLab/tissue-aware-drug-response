#!/usr/bin/env Rscript
# =============================================================================
# Wrong-CV Target Recovery Analysis
# =============================================================================
#
# Computes target recovery under the WRONG (leaked) CV pipeline using existing
# bootstrap results from aggregate_cv_analysis.Rda. Compares with right-CV
# target recovery to quantify how leakage affects biomarker discovery.
#
# Output:
#   results/tables/wrong_cv_target_recovery.csv  (per-drug, per-threshold)
#   results/tables/wrong_vs_right_target_recovery_summary.csv  (aggregated)
#
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

source(here("methods", "00_paths.R"))

cat("=== Wrong-CV Target Recovery Analysis ===\n\n")

# =============================================================================
# 1. Build feature-to-gene mapping
# =============================================================================

cat("Loading feature names from CCLE data...\n")
load(file.path(paths$clean, "ccle.Rda"))
feature_names <- colnames(x)
cat("Total features:", length(feature_names), "\n")
rm(x, y, cls)

# Classify features by type and extract gene symbols
feature_info <- tibble(feature = feature_names) %>%
  mutate(
    type = case_when(
      grepl("_CN$", feature) ~ "CopyNumber",
      grepl("_mut\\.mis$", feature) ~ "Mutation_missense",
      grepl("_Expr$", feature) & !grepl("_\\d+_at_Expr$", feature) ~ "Expression",
      grepl("_mut\\.tr$", feature) ~ "Mutation_truncating",
      grepl("_mut\\.hs$", feature) ~ "Mutation_highsig",
      grepl("_mut\\.hsrec$", feature) ~ "Mutation_recurrent",
      grepl("_\\d+_at_Expr$", feature) ~ "Probe_expression",
      grepl("_Caution$", feature) ~ "QualityFlag",
      grepl("_ValidationFailed$|_ValidationUnavailable$", feature) ~ "QualityFlag",
      TRUE ~ "RPPA"
    ),
    gene = case_when(
      type %in% c("CopyNumber", "Expression",
                   "Mutation_missense", "Mutation_truncating",
                   "Mutation_highsig", "Mutation_recurrent") ~
        sub("_mut\\.hsrec$|_mut\\.hs$|_mut\\.tr$|_mut\\.mis$|_Expr$|_CN$", "", feature),
      TRUE ~ NA_character_
    )
  )

# Apply RPPA mapping
rppa_mapping <- read_csv(here("data", "rppa_gene_mapping.csv"), show_col_types = FALSE)
rppa_lookup <- setNames(rppa_mapping$hgnc_symbol, rppa_mapping$rppa_feature)

feature_info <- feature_info %>%
  mutate(
    gene = ifelse(type == "RPPA" & feature %in% names(rppa_lookup),
                  rppa_lookup[feature], gene)
  )

genes_in_data <- unique(feature_info$gene[!is.na(feature_info$gene)])
cat("Unique genes in feature set:", length(genes_in_data), "\n")

# =============================================================================
# 2. Load drug targets
# =============================================================================

cat("\nLoading drug targets...\n")
drug_targets <- read_csv(here("data", "drug_targets_curated.csv"), show_col_types = FALSE)
gene_mapping <- read_csv(here("data", "gene_name_mapping.csv"), show_col_types = FALSE)

standardize_gene <- function(gene, mapping_df) {
  match_idx <- match(toupper(gene), toupper(mapping_df$alias))
  if (!is.na(match_idx)) return(mapping_df$hgnc_symbol[match_idx])
  return(gene)
}

targets_known <- drug_targets %>%
  filter(!is.na(target_gene), source != "not_found") %>%
  mutate(
    target_gene = sapply(target_gene, standardize_gene, mapping_df = gene_mapping),
    target_in_features = target_gene %in% genes_in_data
  )

cat("Drugs with gene targets:", n_distinct(targets_known$drug_name), "\n")
cat("Drugs with targets in feature set:",
    n_distinct(targets_known$drug_name[targets_known$target_in_features]), "\n")

# =============================================================================
# 3. Load bootstrap results (wrong + right CV)
# =============================================================================

cat("\nLoading bootstrap results from aggregate_cv_analysis.Rda...\n")
load(here("results", "aggregate_cv_analysis.Rda"))
cat("Loaded results for", length(results_list), "drugs\n")

# Helper: map features to genes
map_features_to_genes <- function(features, feature_info_df) {
  feature_info_df %>%
    filter(feature %in% features, !is.na(gene)) %>%
    pull(gene) %>%
    unique()
}

# =============================================================================
# 4. Compute target recovery at multiple thresholds
# =============================================================================

THRESHOLDS <- c(0.005, 0.05, 0.10, 0.25, 0.50, 0.75)

cat("\nComputing target recovery across", length(THRESHOLDS), "thresholds...\n")

all_recovery <- list()

for (thresh in THRESHOLDS) {
  cat(sprintf("  Threshold: %.3f\n", thresh))

  drug_results <- map_dfr(results_list, function(r) {
    drug <- r$drug

    # Get targets for this drug
    drug_tgts <- targets_known %>%
      filter(drug_name == drug, target_in_features == TRUE) %>%
      pull(target_gene) %>%
      unique()

    if (length(drug_tgts) == 0) {
      return(tibble(
        drug = drug, threshold = thresh, n_targets = 0,
        wrong_any = NA, wrong_full = NA, wrong_n_found = NA_integer_,
        right_any = NA, right_full = NA, right_n_found = NA_integer_,
        wrong_n_selected = NA_integer_, right_n_selected = NA_integer_,
        wrong_targets_found = NA_character_, right_targets_found = NA_character_
      ))
    }

    # --- Wrong CV ---
    wrong_freq <- r$wrong_bs_cnt  # Already frequencies (0-1)
    wrong_selected <- names(wrong_freq)[wrong_freq >= thresh]
    wrong_genes <- map_features_to_genes(wrong_selected, feature_info)
    wrong_found <- intersect(drug_tgts, wrong_genes)

    # --- Right CV ---
    right_freq <- r$right_bs_cnt  # Already frequencies (0-1)
    right_selected <- names(right_freq)[right_freq >= thresh]
    right_genes <- map_features_to_genes(right_selected, feature_info)
    right_found <- intersect(drug_tgts, right_genes)

    tibble(
      drug = drug,
      threshold = thresh,
      n_targets = length(drug_tgts),
      wrong_any = length(wrong_found) > 0,
      wrong_full = length(wrong_found) == length(drug_tgts),
      wrong_n_found = length(wrong_found),
      right_any = length(right_found) > 0,
      right_full = length(right_found) == length(drug_tgts),
      right_n_found = length(right_found),
      wrong_n_selected = length(wrong_selected),
      right_n_selected = length(right_selected),
      wrong_targets_found = paste(wrong_found, collapse = ";"),
      right_targets_found = paste(right_found, collapse = ";")
    )
  })

  all_recovery[[as.character(thresh)]] <- drug_results
}

recovery_df <- bind_rows(all_recovery)

# Filter to drugs with targets in feature set
recovery_valid <- recovery_df %>%
  filter(n_targets > 0)

cat("\nDrugs with targets:", n_distinct(recovery_valid$drug), "\n")

# =============================================================================
# 5. Summary statistics
# =============================================================================

summary_df <- recovery_valid %>%
  group_by(threshold) %>%
  summarize(
    n_drugs = n_distinct(drug),
    wrong_any_pct = round(100 * mean(wrong_any, na.rm = TRUE), 1),
    wrong_any_n = sum(wrong_any, na.rm = TRUE),
    right_any_pct = round(100 * mean(right_any, na.rm = TRUE), 1),
    right_any_n = sum(right_any, na.rm = TRUE),
    wrong_full_pct = round(100 * mean(wrong_full, na.rm = TRUE), 1),
    right_full_pct = round(100 * mean(right_full, na.rm = TRUE), 1),
    mean_wrong_n_selected = round(mean(wrong_n_selected, na.rm = TRUE), 1),
    mean_right_n_selected = round(mean(right_n_selected, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("\n=== Target Recovery Summary (Any Target Found) ===\n\n")
cat(sprintf("%-10s  %6s  %12s  %12s  %12s  %12s\n",
            "Threshold", "Drugs", "Wrong %", "Wrong n",  "Right %", "Right n"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(summary_df))) {
  s <- summary_df[i, ]
  cat(sprintf("%-10.3f  %6d  %11.1f%%  %11d  %11.1f%%  %11d\n",
              s$threshold, s$n_drugs,
              s$wrong_any_pct, s$wrong_any_n,
              s$right_any_pct, s$right_any_n))
}

cat("\n=== Mean Features Selected ===\n\n")
cat(sprintf("%-10s  %15s  %15s\n", "Threshold", "Wrong (mean)", "Right (mean)"))
cat(strrep("-", 45), "\n")
for (i in seq_len(nrow(summary_df))) {
  s <- summary_df[i, ]
  cat(sprintf("%-10.3f  %15.1f  %15.1f\n",
              s$threshold, s$mean_wrong_n_selected, s$mean_right_n_selected))
}

# =============================================================================
# 6. Save results
# =============================================================================

table_dir <- file.path(paths$results, "tables")

write_csv(recovery_valid,
          file.path(table_dir, "wrong_cv_target_recovery.csv"))
cat("\nSaved:", file.path(table_dir, "wrong_cv_target_recovery.csv"), "\n")

write_csv(summary_df,
          file.path(table_dir, "wrong_vs_right_target_recovery_summary.csv"))
cat("Saved:", file.path(table_dir, "wrong_vs_right_target_recovery_summary.csv"), "\n")

cat("\nDone.\n")
