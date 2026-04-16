#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(caret)
})

source(here("methods", "00_paths.R"))
source(here("methods", "01_data_utils.R"))
source(here("methods", "comparison_common.R"))

GLOBAL_SEED <- 42L
CV_NFOLDS <- 5L
MIN_SAMPLES_PER_TISSUE <- 25L

out_dir <- ensure_comparison_dir()
fold_file <- file.path(out_dir, "fold_assignments.rds")
index_file <- file.path(out_dir, "fold_index.csv")

message("Loading main data...")
data <- load_main_data("Sanger")
all_drugs <- get_drug_list(data)

fold_assignments <- list()
index_rows <- vector("list", length(all_drugs))

for (i in seq_along(all_drugs)) {
  drug <- all_drugs[i]
  message(sprintf("[%s] Drug %d/%d: %s", Sys.time(), i, length(all_drugs), drug))

  entry <- tryCatch({
    dd <- prepare_drug_data_filtered(data, drug, min_per_tissue = MIN_SAMPLES_PER_TISSUE)

    if (is.null(dd)) {
      list(
        fold_entry = list(
          drug_name = drug,
          status = "skipped",
          reason = "insufficient_valid_tissues"
        ),
        index = data.frame(
          drug = drug,
          n_samples = NA_integer_,
          n_tissues = NA_integer_,
          status = "skipped",
          stringsAsFactors = FALSE
        )
      )
    } else {
      set.seed(GLOBAL_SEED)
      folds <- caret::createFolds(dd$tissue, k = CV_NFOLDS, list = FALSE)

      list(
        fold_entry = list(
          drug_name = drug,
          status = "ok",
          cell_line_ids = rownames(dd$x),
          tissue = as.character(dd$tissue),
          foldid = as.integer(folds),
          n_samples = nrow(dd$x),
          n_tissues = length(unique(dd$tissue))
        ),
        index = data.frame(
          drug = drug,
          n_samples = nrow(dd$x),
          n_tissues = length(unique(dd$tissue)),
          status = "ok",
          stringsAsFactors = FALSE
        )
      )
    }
  }, error = function(e) {
    list(
      fold_entry = list(
        drug_name = drug,
        status = "error",
        reason = e$message
      ),
      index = data.frame(
        drug = drug,
        n_samples = NA_integer_,
        n_tissues = NA_integer_,
        status = "error",
        stringsAsFactors = FALSE
      )
    )
  })

  fold_assignments[[drug]] <- entry$fold_entry
  index_rows[[i]] <- entry$index
}

index_df <- do.call(rbind, index_rows)
saveRDS(fold_assignments, fold_file)
write.csv(index_df, index_file, row.names = FALSE)

ok_n <- sum(index_df$status == "ok")
skip_n <- sum(index_df$status == "skipped")
err_n <- sum(index_df$status == "error")

message("Fold assignment export complete")
message(sprintf("  OK: %d", ok_n))
message(sprintf("  Skipped: %d", skip_n))
message(sprintf("  Errors: %d", err_n))
message(sprintf("  Saved fold file: %s", fold_file))
message(sprintf("  Saved index file: %s", index_file))

append_progress(
  phase = "Phase 0.1 - Fold assignments",
  status = "completed",
  details = c(
    sprintf("Drugs total: %d", nrow(index_df)),
    sprintf("Drugs with valid folds: %d", ok_n),
    sprintf("Drugs skipped: %d", skip_n),
    sprintf("Drugs errored: %d", err_n)
  )
)
