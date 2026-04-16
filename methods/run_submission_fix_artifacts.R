#!/usr/bin/env Rscript
# Purpose: Generate post-audit submission artifacts without model refits.
# Outputs:
#   - results/tables/threshold_sensitivity_jaccard.csv
#   - results/tables/threshold_sensitivity_target.csv
#   - results/tables/target_mapping_summary.csv

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(purrr)
})

source(here("methods", "00_paths.R"))

table_dir <- file.path(paths$results, "tables")
if (!dir.exists(table_dir)) {
  dir.create(table_dir, recursive = TRUE)
}

jaccard_thresholds <- c(0.50, 0.60, 0.70, 0.80, 0.90)
target_thresholds <- c(0.25, 0.50, 0.75)

jaccard_out <- file.path(table_dir, "threshold_sensitivity_jaccard.csv")
target_out <- file.path(table_dir, "threshold_sensitivity_target.csv")
mapping_out <- file.path(table_dir, "target_mapping_summary.csv")

get_freq <- function(bs_obj) {
  if (!is.null(bs_obj$freq)) {
    return(bs_obj$freq)
  }
  if (!is.null(bs_obj$cnt) && !is.null(bs_obj$nboot)) {
    return(bs_obj$cnt / bs_obj$nboot)
  }
  stop("Bootstrap object missing both freq and cnt/nboot.")
}

compute_jaccard <- function(feat_a, feat_b) {
  n_union <- length(union(feat_a, feat_b))
  if (n_union == 0) return(NA_real_)
  length(intersect(feat_a, feat_b)) / n_union
}

cat(strrep("=", 80), "\n")
cat("Generating submission fix artifacts\n")
cat(strrep("=", 80), "\n")

# ------------------------------------------------------------------------------
# Fix 3 (part 1): Jaccard threshold sensitivity
# ------------------------------------------------------------------------------
jaccard_source_candidates <- c(
  file.path(paths$scratch, "feature_stability_results.Rda"),
  file.path(paths$scratch, "feature_stability_checkpoint.Rda")
)
jaccard_source <- jaccard_source_candidates[file.exists(jaccard_source_candidates)][1]

if (is.na(jaccard_source) || !nzchar(jaccard_source)) {
  stop("No feature-stability bootstrap file found in scratch.")
}

load(jaccard_source)
if (!exists("all_results")) {
  stop("Expected object 'all_results' not found in: ", jaccard_source)
}

jaccard_inputs <- keep(all_results, function(x) {
  is.list(x) &&
    identical(x$status, "success") &&
    !is.null(x$bs_wrong) &&
    !is.null(x$bs_right)
})

if (length(jaccard_inputs) == 0) {
  stop("No successful bootstrap entries found for Jaccard sensitivity.")
}

jaccard_table <- map_dfr(jaccard_thresholds, function(thr) {
  per_drug <- map_dfr(jaccard_inputs, function(res) {
    freq_wrong <- get_freq(res$bs_wrong)
    freq_right <- get_freq(res$bs_right)
    feat_wrong <- names(freq_wrong)[freq_wrong >= thr]
    feat_right <- names(freq_right)[freq_right >= thr]

    tibble(
      jaccard = compute_jaccard(feat_wrong, feat_right),
      n_feat_wrong = length(feat_wrong),
      n_feat_right = length(feat_right)
    )
  })

  tibble(
    threshold = thr,
    mean_jaccard = mean(per_drug$jaccard, na.rm = TRUE),
    median_jaccard = median(per_drug$jaccard, na.rm = TRUE),
    pct_zero_overlap = mean(per_drug$jaccard == 0, na.rm = TRUE) * 100,
    mean_n_feat_wrong = mean(per_drug$n_feat_wrong, na.rm = TRUE),
    mean_n_feat_right = mean(per_drug$n_feat_right, na.rm = TRUE)
  )
})

write_csv(jaccard_table, jaccard_out)
cat("Wrote:", jaccard_out, "\n")

# ------------------------------------------------------------------------------
# Shared setup for target-recovery sensitivity + mapping summary
# ------------------------------------------------------------------------------
target_file <- here("data", "drug_targets_curated.csv")
mapping_file <- here("data", "gene_name_mapping.csv")
bootstrap_file <- here("results", "dsen_results", "dsen_bootstrap_results.Rda")

if (!file.exists(target_file)) {
  stop("Missing target curation file: ", target_file)
}
if (!file.exists(bootstrap_file)) {
  stop("Missing DSEN bootstrap file: ", bootstrap_file)
}

targets_raw <- read_csv(target_file, show_col_types = FALSE)
targets_known <- targets_raw %>%
  filter(!is.na(target_gene), source != "not_found")

if (file.exists(mapping_file)) {
  gene_mapping <- read_csv(mapping_file, show_col_types = FALSE)
} else {
  gene_mapping <- tibble(alias = character(), hgnc_symbol = character())
}

standardize_gene <- function(gene) {
  match_idx <- match(toupper(gene), toupper(gene_mapping$alias))
  if (!is.na(match_idx)) {
    return(gene_mapping$hgnc_symbol[match_idx])
  }
  gene
}

targets_known <- targets_known %>%
  mutate(target_gene = vapply(target_gene, standardize_gene, character(1)))

load(file.path(paths$clean, "ccle.Rda"))
feature_names <- colnames(x)
genes_in_data <- unique(sub(
  "_mut\\.hsrec$|_mut$|_Expr$|_CN$|_del$|_amp$",
  "",
  feature_names
))
rm(x, y, cls, feature_names)

targets_in_features <- targets_known %>%
  mutate(target_in_features = target_gene %in% genes_in_data)

# ------------------------------------------------------------------------------
# Fix 5: Target-mapping summary artifact
# ------------------------------------------------------------------------------
mapping_summary <- tibble(
  total_curated_pairs = nrow(targets_known),
  unique_drugs = n_distinct(targets_known$drug_name),
  unique_targets = n_distinct(targets_known$target_gene),
  mapped_pairs = sum(targets_in_features$target_in_features, na.rm = TRUE)
) %>%
  mutate(mapping_rate = mapped_pairs / total_curated_pairs)

write_csv(mapping_summary, mapping_out)
cat("Wrote:", mapping_out, "\n")

# ------------------------------------------------------------------------------
# Fix 3 (part 2): Target-recovery threshold sensitivity
# ------------------------------------------------------------------------------
load(bootstrap_file)
if (!exists("all_results")) {
  stop("Expected object 'all_results' not found in: ", bootstrap_file)
}

bootstrap_results <- all_results

compute_recovery_for_threshold <- function(thr) {
  per_drug <- map_dfr(bootstrap_results, function(res) {
    if (!is.list(res) || !identical(res$status, "success")) {
      return(tibble(
        n_targets = 0L,
        en_any = NA,
        en_full = NA,
        dsen_any = NA,
        dsen_full = NA
      ))
    }

    drug_targets <- targets_in_features %>%
      filter(drug_name == res$drug, target_in_features) %>%
      pull(target_gene) %>%
      unique()

    if (length(drug_targets) == 0) {
      return(tibble(
        n_targets = 0L,
        en_any = NA,
        en_full = NA,
        dsen_any = NA,
        dsen_full = NA
      ))
    }

    en_freq <- get_freq(res$en_bootstrap)
    dsen_freq <- get_freq(res$dsen_bootstrap)

    en_selected <- names(en_freq)[en_freq >= thr]
    dsen_selected <- names(dsen_freq)[dsen_freq >= thr]

    en_genes <- unique(sub(
      "_mut\\.hsrec$|_mut$|_Expr$|_CN$|_del$|_amp$",
      "",
      en_selected
    ))
    dsen_genes <- unique(sub(
      "_mut\\.hsrec$|_mut$|_Expr$|_CN$|_del$|_amp$",
      "",
      dsen_selected
    ))

    en_n_found <- length(intersect(drug_targets, en_genes))
    dsen_n_found <- length(intersect(drug_targets, dsen_genes))
    n_targets <- length(drug_targets)

    tibble(
      n_targets = n_targets,
      en_any = en_n_found > 0,
      en_full = en_n_found == n_targets,
      dsen_any = dsen_n_found > 0,
      dsen_full = dsen_n_found == n_targets
    )
  }) %>%
    filter(n_targets > 0)

  tibble(
    threshold = thr,
    en_any_recovery = mean(per_drug$en_any, na.rm = TRUE),
    en_full_recovery = mean(per_drug$en_full, na.rm = TRUE),
    dsen_any_recovery = mean(per_drug$dsen_any, na.rm = TRUE),
    dsen_full_recovery = mean(per_drug$dsen_full, na.rm = TRUE)
  )
}

target_table <- map_dfr(target_thresholds, compute_recovery_for_threshold)
write_csv(target_table, target_out)
cat("Wrote:", target_out, "\n")

cat("Done.\n")
