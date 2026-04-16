#' =============================================================================
#' Drug Response Prediction - Data Utilities
#' =============================================================================
#'
#' This script provides functions for loading and preprocessing drug response
#' and genomic feature data.
#'
#' Usage:
#'   library(here)
#'   source(here("methods", "00_paths.R"))
#'   source(here("methods", "01_data_utils.R"))
#'
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Data Loading Functions
#' -----------------------------------------------------------------------------

#' Load the main drug response and predictor data
#'
#' @param data_source Character, either "Sanger" or "CCLE"
#' @return List containing response data (dat), tissue info, and predictors (x)
load_main_data <- function(data_source = "Sanger") {
  if (data_source == "Sanger") {
    # Load Sanger/GDSC data
    sanger_file <- file.path(paths$clean, "Sanger.Rdata")
    ccle_file <- file.path(paths$clean, "ccle.Rda")

    if (!file.exists(sanger_file)) {
      stop("Sanger.Rdata not found in: ", paths$clean)
    }
    if (!file.exists(ccle_file)) {
      stop("ccle.Rda not found in: ", paths$clean)
    }

    load(sanger_file)  # Loads 'Sanger' object
    load(ccle_file)    # Loads 'x', 'y', 'cls' objects

    message(sprintf("Loaded Sanger data: %d cell lines, %d drugs",
                    nrow(Sanger$dat), ncol(Sanger$dat)))
    message(sprintf("Loaded predictors: %d features", ncol(x)))

    return(list(
      response = Sanger,
      predictors = x,
      data_source = "Sanger"
    ))
  } else {
    stop("Data source not implemented: ", data_source)
  }
}

#' Get list of available drugs
#'
#' @param data List from load_main_data()
#' @return Character vector of drug names
get_drug_list <- function(data) {
  colnames(data$response$dat)
}

#' Get tissue types
#'
#' @param data List from load_main_data()
#' @return Character vector of unique tissue types
get_tissue_types <- function(data) {
  unique(data$response$tissue[!is.na(data$response$tissue)])
}

#' -----------------------------------------------------------------------------
#' Data Preparation Functions
#' -----------------------------------------------------------------------------

#' Prepare data for a single drug
#'
#' @param data List from load_main_data()
#' @param drug Character, drug name
#' @param min_samples Integer, minimum samples required
#' @return List with x (predictors), y (response), tissue (tissue labels)
prepare_drug_data <- function(data, drug, min_samples = 30) {
  # Get drug response
  y <- data$response$dat[, drug]
  names(y) <- rownames(data$response$dat)

  # Get tissue info
  tissue <- data$response$tissue
  names(tissue) <- rownames(data$response$dat)

  # Remove NAs
  valid_idx <- !is.na(y)
  y <- y[valid_idx]
  tissue <- tissue[valid_idx]

  # Match with predictors
  common_samples <- intersect(names(y), rownames(data$predictors))

  if (length(common_samples) < min_samples) {
    warning(sprintf("Drug %s has only %d samples (minimum: %d)",
                    drug, length(common_samples), min_samples))
    return(NULL)
  }

  # Subset and align
  y <- y[common_samples]
  tissue <- tissue[common_samples]
  x <- data$predictors[common_samples, ]

  # Order by sample name for reproducibility
  ord <- order(names(y))
  y <- y[ord]
  tissue <- tissue[ord]
  x <- x[ord, ]

  return(list(
    x = x,
    y = y,
    tissue = tissue,
    drug = drug,
    n_samples = length(y),
    n_features = ncol(x)
  ))
}

#' Prepare data with tissue filtering
#'
#' @param data List from load_main_data()
#' @param drug Character, drug name
#' @param min_per_tissue Integer, minimum samples per tissue
#' @return List with x, y, tissue (filtered to tissues with enough samples)
prepare_drug_data_filtered <- function(data, drug, min_per_tissue = 5) {
  # First get basic data
  drug_data <- prepare_drug_data(data, drug)

  if (is.null(drug_data)) return(NULL)

  # Count samples per tissue
  tissue_counts <- table(drug_data$tissue)
  valid_tissues <- names(tissue_counts)[tissue_counts >= min_per_tissue]

  if (length(valid_tissues) < 2) {
    warning(sprintf("Drug %s has fewer than 2 tissues with >= %d samples",
                    drug, min_per_tissue))
    return(NULL)
  }

  # Filter to valid tissues
  keep_idx <- drug_data$tissue %in% valid_tissues
  drug_data$x <- drug_data$x[keep_idx, ]
  drug_data$y <- drug_data$y[keep_idx]
  drug_data$tissue <- drug_data$tissue[keep_idx]
  drug_data$n_samples <- sum(keep_idx)
  drug_data$n_tissues <- length(valid_tissues)
  drug_data$tissues <- valid_tissues

  return(drug_data)
}

#' Prepare drug data with tissue filtering
#' (Note: This is different from the dsen_utils.R prepare_tissue_data function)
#'
#' @param data List from load_main_data()
#' @param drug Character, drug name
#' @param min_per_tissue Integer, minimum samples per tissue
#' @return List with x, y, tissue (filtered to tissues with enough samples)
prepare_drug_with_tissue <- function(data, drug, min_per_tissue = 15) {
  prepare_drug_data_filtered(data, drug, min_per_tissue)
}

#' -----------------------------------------------------------------------------
#' Feature Preprocessing Functions
#' -----------------------------------------------------------------------------

#' Convert matrix to numeric (handles character matrices)
#'
#' @param x Matrix
#' @return Numeric matrix
as_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  mode(x) <- "numeric"
  return(x)
}

#' Preprocess features: variance filter, correlation filter, scale
#'
#' @param x Feature matrix (samples x features)
#' @param y Response vector
#' @param var_thresh Variance threshold
#' @param cor_thresh Correlation threshold
#' @param scale_features Logical, whether to scale features
#' @return List with filtered x, selected feature names, and scaling parameters
preprocess_features <- function(x, y,
                                 var_thresh = VARIANCE_THRESHOLD,
                                 cor_thresh = CORRELATION_THRESHOLD,
                                 scale_features = TRUE) {

  x <- as_numeric_matrix(x)

  # Replace NAs with 0 (common for missing mutation data)
  x[is.na(x)] <- 0

  # Step 1: Variance filter
  col_vars <- apply(x, 2, var, na.rm = TRUE)
  high_var_idx <- which(col_vars > var_thresh)
  x <- x[, high_var_idx, drop = FALSE]

  # Step 2: Correlation filter
  cors <- apply(x, 2, function(col) cor(col, y, use = "pairwise.complete.obs"))
  high_cor_idx <- which(abs(cors) > cor_thresh)
  x <- x[, high_cor_idx, drop = FALSE]

  # Step 3: Remove duplicate columns
  x <- x[, !duplicated(t(x)), drop = FALSE]

  # Step 4: Scale (optional)
  scale_center <- NULL
  scale_scale <- NULL

  if (scale_features && ncol(x) > 0) {
    # Identify mutation columns (binary - don't scale these)
    is_mutation <- grepl("_mut(\\.|$)|^mut_", colnames(x), ignore.case = TRUE)

    if (any(!is_mutation)) {
      x_to_scale <- x[, !is_mutation, drop = FALSE]
      x_scaled <- scale(x_to_scale)
      scale_center <- attr(x_scaled, "scaled:center")
      scale_scale <- attr(x_scaled, "scaled:scale")
      x[, !is_mutation] <- x_scaled
    }
  }

  return(list(
    x = x,
    selected_features = colnames(x),
    n_features = ncol(x),
    scale_center = scale_center,
    scale_scale = scale_scale
  ))
}

#' Apply preprocessing from training data to test data
#'
#' @param x_test Test feature matrix
#' @param preproc List from preprocess_features() on training data
#' @return Preprocessed test matrix
apply_preprocessing <- function(x_test, preproc) {
  x_test <- as_numeric_matrix(x_test)
  x_test[is.na(x_test)] <- 0


  # Keep only features selected in training
  common_features <- intersect(colnames(x_test), preproc$selected_features)
  missing_features <- setdiff(preproc$selected_features, colnames(x_test))

  if (length(missing_features) > 0) {
    warning(sprintf("%d features from training not found in test data",
                    length(missing_features)))
  }

  x_test <- x_test[, common_features, drop = FALSE]

  # Apply scaling using training parameters
  if (!is.null(preproc$scale_center)) {
    is_mutation <- grepl("_mut(\\.|$)|^mut_", colnames(x_test), ignore.case = TRUE)

    if (any(!is_mutation)) {
      non_mut_features <- colnames(x_test)[!is_mutation]
      for (feat in non_mut_features) {
        if (feat %in% names(preproc$scale_center)) {
          x_test[, feat] <- (x_test[, feat] - preproc$scale_center[feat]) /
                             preproc$scale_scale[feat]
        }
      }
    }
  }

  return(x_test)
}

#' -----------------------------------------------------------------------------
#' Summary Functions
#' -----------------------------------------------------------------------------

#' Print summary of loaded data
#'
#' @param data List from load_main_data()
print_data_summary <- function(data) {
  cat("\n=== Data Summary ===\n")
  cat(sprintf("Data source: %s\n", data$data_source))
  cat(sprintf("Total cell lines: %d\n", nrow(data$response$dat)))
  cat(sprintf("Total drugs: %d\n", ncol(data$response$dat)))
  cat(sprintf("Total features: %d\n", ncol(data$predictors)))

  # Feature breakdown
  feature_types <- c(
    Expression = sum(grepl("_Expr$", colnames(data$predictors))),
    Mutation = sum(grepl("_mut(\\.|$)", colnames(data$predictors))),
    CopyNumber = sum(grepl("_CN$", colnames(data$predictors))),
    Protein = sum(grepl("_RPPA$|^[A-Z0-9]+_p[A-Z]", colnames(data$predictors)))
  )
  feature_types["Other"] <- ncol(data$predictors) - sum(feature_types)

  cat("\nFeature types:\n")
  for (type in names(feature_types)) {
    cat(sprintf("  %s: %d\n", type, feature_types[type]))
  }

  # Tissue distribution
  tissue_counts <- sort(table(data$response$tissue), decreasing = TRUE)
  cat(sprintf("\nTissue types: %d\n", length(tissue_counts)))
  cat("Top 5 tissues:\n")
  for (i in 1:min(5, length(tissue_counts))) {
    cat(sprintf("  %s: %d\n", names(tissue_counts)[i], tissue_counts[i]))
  }

  invisible(data)
}

message("Data utilities loaded.")
