#!/usr/bin/env Rscript
#' =============================================================================
#' Process GDSC Drug Targets
#' =============================================================================
#'
#' Downloads and processes drug target annotations from GDSC to create a
#' curated drug targets file for feature selection validation.
#'
#' Input:  data/screened_compounds_gdsc.csv (downloaded from GDSC)
#' Output: data/drug_targets_curated.csv
#'
#' =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

source(here("methods", "00_paths.R"))

# Paths
data_dir <- here("data")
results_dir <- here("results")

gdsc_file <- file.path(data_dir, "screened_compounds_gdsc.csv")
output_file <- file.path(data_dir, "drug_targets_curated.csv")
drug_list_file <- file.path(results_dir, "drug_list.txt")

cat("=== Processing GDSC Drug Targets ===\n\n")

# Check if GDSC file exists
if (!file.exists(gdsc_file)) {
  stop("GDSC file not found. Run: curl -o data/screened_compounds_gdsc.csv ",
       "'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/screened_compounds_rel_8.5.csv'")
}

# Load GDSC drug annotations
gdsc <- read_csv(gdsc_file, show_col_types = FALSE)
cat("Loaded GDSC annotations:", nrow(gdsc), "compounds\n")

# Load our drug list
if (file.exists(drug_list_file)) {
  our_drugs <- read_lines(drug_list_file)
  cat("Loaded our drug list:", length(our_drugs), "drugs\n")
} else {
  # Generate from Sanger data
  load(file.path(paths$clean, "Sanger.Rdata"))
  our_drugs <- colnames(Sanger$dat)
  write_lines(our_drugs, drug_list_file)
  cat("Generated drug list from Sanger data:", length(our_drugs), "drugs\n")
}

# Create lookup for GDSC drugs (case-insensitive)
gdsc_lookup <- gdsc %>%
  mutate(
    drug_lower = tolower(DRUG_NAME),
    synonyms_lower = tolower(SYNONYMS)
  )

# Function to find a drug in GDSC by name or synonym
find_drug <- function(drug_name, gdsc_df) {
  name_lower <- tolower(drug_name)

  # Try exact match
  match <- gdsc_df %>% filter(drug_lower == name_lower)
  if (nrow(match) > 0) return(match[1, ])

  # Try synonym match
  match <- gdsc_df %>%
    filter(str_detect(synonyms_lower, fixed(name_lower, ignore_case = TRUE)))
  if (nrow(match) > 0) return(match[1, ])

  # Try partial match on name
  match <- gdsc_df %>%
    filter(str_detect(drug_lower, fixed(name_lower, ignore_case = TRUE)) |
           str_detect(name_lower, fixed(drug_lower, ignore_case = TRUE)))
  if (nrow(match) > 0) return(match[1, ])

  return(NULL)
}

# Gene symbol standardization
standardize_gene <- function(gene) {
  gene <- trimws(gene)
  gene <- toupper(gene)

  # Common mappings
  mappings <- c(
    "ABL" = "ABL1",
    "C-MET" = "MET",
    "HER2" = "ERBB2",
    "RAF" = "BRAF",  # Could also be RAF1
    "ERBB" = "ERBB2",
    "PDGFR" = "PDGFRA",  # or PDGFRB
    "VEGFR" = "KDR",     # VEGFR2 = KDR
    "FGFR" = "FGFR1",    # Could be 1-4
    "IGF1R" = "IGF1R",
    "AURORA" = "AURKA",  # or AURKB
    "MEK" = "MAP2K1",
    "ERK" = "MAPK1",
    "JAK" = "JAK1",      # Could be JAK1/2/3
    "SRC" = "SRC",
    "BCL2" = "BCL2",
    "MDM2" = "MDM2",
    "HDAC" = "HDAC1",    # Could be multiple
    "HSP90" = "HSP90AA1",
    "MTORC1" = "MTOR",
    "MTORC2" = "MTOR",
    "PI3K" = "PIK3CA",   # Could be multiple
    "AKT" = "AKT1"
  )

  if (gene %in% names(mappings)) {
    return(mappings[gene])
  }
  return(gene)
}

# Process each of our drugs
cat("\nProcessing drugs...\n")
results <- list()

for (drug in our_drugs) {
  match <- find_drug(drug, gdsc_lookup)

  if (is.null(match)) {
    # No match found
    results[[length(results) + 1]] <- tibble(
      drug_name = drug,
      target_gene = NA_character_,
      target_type = NA_character_,
      target_pathway = NA_character_,
      source = "not_found",
      gdsc_name = NA_character_
    )
  } else {
    # Parse targets (comma-separated)
    targets_raw <- match$TARGET
    pathway <- match$TARGET_PATHWAY
    gdsc_name <- match$DRUG_NAME

    if (is.na(targets_raw) || targets_raw == "") {
      results[[length(results) + 1]] <- tibble(
        drug_name = drug,
        target_gene = NA_character_,
        target_type = "unknown",
        target_pathway = pathway,
        source = "GDSC",
        gdsc_name = gdsc_name
      )
    } else {
      # Split targets
      targets <- str_split(targets_raw, "[,;]")[[1]]
      targets <- trimws(targets)
      targets <- targets[targets != ""]

      for (i in seq_along(targets)) {
        target_std <- standardize_gene(targets[i])
        target_type <- ifelse(i == 1, "primary", "secondary")

        results[[length(results) + 1]] <- tibble(
          drug_name = drug,
          target_gene = target_std,
          target_type = target_type,
          target_pathway = pathway,
          source = "GDSC",
          gdsc_name = gdsc_name
        )
      }
    }
  }
}

# Combine results
drug_targets <- bind_rows(results)

# Summary
cat("\n=== Summary ===\n")
cat("Total rows:", nrow(drug_targets), "\n")
cat("Unique drugs:", n_distinct(drug_targets$drug_name), "\n")
cat("Drugs with targets:", sum(!is.na(drug_targets$target_gene)), "\n")
cat("Drugs not found in GDSC:", sum(drug_targets$source == "not_found"), "\n")

drugs_with_targets <- drug_targets %>%
  filter(!is.na(target_gene)) %>%
  distinct(drug_name) %>%
  nrow()
cat("Unique drugs with at least one target:", drugs_with_targets, "\n")

# Show distribution
cat("\nTarget pathway distribution:\n")
drug_targets %>%
  filter(!is.na(target_pathway)) %>%
  count(target_pathway, sort = TRUE) %>%
  head(10) %>%
  print()

# Save
write_csv(drug_targets, output_file)
cat("\nSaved to:", output_file, "\n")

# Also save a simplified version with just primary targets
primary_targets <- drug_targets %>%
  filter(target_type == "primary" | is.na(target_type)) %>%
  select(drug_name, target_gene, target_pathway, source)

write_csv(primary_targets, file.path(data_dir, "drug_targets_primary.csv"))
cat("Primary targets saved to:", file.path(data_dir, "drug_targets_primary.csv"), "\n")
