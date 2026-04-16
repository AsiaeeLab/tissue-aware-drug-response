#' =============================================================================
#' Drug Response Prediction - Parallel Processing Utilities
#' =============================================================================
#'
#' This script provides functions for running analyses in parallel.
#'
#' Usage:
#'   library(here)
#'   source(here("methods", "05_parallel_utils.R"))
#'
#' =============================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(foreach)
  library(doParallel)
})

#' -----------------------------------------------------------------------------
#' Parallel Setup
#' -----------------------------------------------------------------------------

#' Initialize parallel backend
#'
#' @param ncores Number of cores to use
#' @return Cluster object
setup_parallel <- function(ncores = NCORES) {
  # Check available cores
  max_cores <- detectCores()
  ncores <- min(ncores, max_cores - 1)  # Leave one core free

  message(sprintf("Setting up parallel backend with %d cores (max available: %d)",
                  ncores, max_cores))

  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  return(cl)
}

#' Stop parallel backend
#'
#' @param cl Cluster object
stop_parallel <- function(cl) {
  stopCluster(cl)
  message("Parallel backend stopped.")
}

#' -----------------------------------------------------------------------------
#' Parallel CV Functions
#' -----------------------------------------------------------------------------

#' Run CV analysis for multiple drugs in parallel
#'
#' @param data Data list from load_main_data()
#' @param drugs Vector of drug names
#' @param output_dir Directory to save results
#' @param ncores Number of cores
#' @param overwrite Whether to overwrite existing results
#' @return Data frame with summary for all drugs
run_cv_parallel <- function(data, drugs, output_dir = results_dirs$cv_results,
                             ncores = NCORES, overwrite = FALSE) {

  # Filter drugs that need processing
  if (!overwrite) {
    existing <- sapply(drugs, function(d) {
      file.exists(file.path(output_dir, paste0("cv_", gsub("/", "_", d), ".Rda")))
    })
    drugs_to_run <- drugs[!existing]
    message(sprintf("%d drugs already processed, %d to run",
                    sum(existing), length(drugs_to_run)))
  } else {
    drugs_to_run <- drugs
  }

  if (length(drugs_to_run) == 0) {
    message("All drugs already processed.")
    return(NULL)
  }

  # Set up parallel backend
  cl <- setup_parallel(ncores)

  # Export required objects and functions to workers
  clusterExport(cl, c("data", "output_dir", "paths",
                      "VARIANCE_THRESHOLD", "CORRELATION_THRESHOLD",
                      "CV_NFOLDS", "BOOTSTRAP_NREPS", "ALPHA_GRID", "GLOBAL_SEED"),
                envir = environment())

  # Export functions
  clusterEvalQ(cl, {
    source("R/00_config.R")
    source("R/01_data_utils.R")
    source("R/02_cv_utils.R")
  })

  # Process drugs in parallel
  results <- foreach(drug = drugs_to_run,
                      .combine = rbind,
                      .errorhandling = "pass",
                      .packages = c("glmnet", "caret")) %dopar% {

    tryCatch({
      # Prepare data for this drug
      drug_data <- prepare_drug_data(data, drug)

      if (is.null(drug_data)) {
        return(data.frame(drug = drug, status = "skipped_insufficient_data",
                          mse_wrong = NA, mse_right = NA, stringsAsFactors = FALSE))
      }

      # Run wrong CV
      wrong_result <- wrong_cv_tune_alpha(drug_data$x, drug_data$y, seed = GLOBAL_SEED)

      # Run right CV
      right_result <- cv_tune_alpha(drug_data$x, drug_data$y, seed = GLOBAL_SEED)

      # Bootstrap for wrong CV
      bs_wrong <- bootstrap_features(drug_data$x, drug_data$y,
                                      alpha = wrong_result$alpha.min,
                                      lambda = wrong_result$lambda.min,
                                      seed = GLOBAL_SEED)

      # Bootstrap for right CV
      bs_right <- bootstrap_features(drug_data$x, drug_data$y,
                                      alpha = right_result$alpha.min,
                                      lambda = right_result$lambda.min,
                                      seed = GLOBAL_SEED)

      # Compare results
      comparison <- compare_cv_results(wrong_result, right_result,
                                        bs_wrong, bs_right)

      # Save results
      filename <- paste0("cv_", gsub("/", "_", drug), ".Rda")
      save(drug_data, wrong_result, right_result, bs_wrong, bs_right, comparison,
           file = file.path(output_dir, filename))

      # Return summary
      data.frame(
        drug = drug,
        status = "completed",
        n_samples = drug_data$n_samples,
        mse_wrong = comparison$mse_wrong,
        mse_right = comparison$mse_right,
        mse_pct_increase = comparison$mse_pct_increase,
        n_features_wrong = comparison$n_features_wrong,
        n_features_right = comparison$n_features_right,
        jaccard = comparison$jaccard,
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      data.frame(drug = drug, status = paste("error:", e$message),
                 mse_wrong = NA, mse_right = NA, stringsAsFactors = FALSE)
    })
  }

  # Stop parallel backend
  stop_parallel(cl)

  return(results)
}

#' Run DSEN analysis for multiple drugs in parallel
#'
#' @param data Data list from load_main_data()
#' @param drugs Vector of drug names
#' @param output_dir Directory to save results
#' @param ncores Number of cores
#' @param overwrite Whether to overwrite existing results
#' @return Data frame with summary
run_dsen_parallel <- function(data, drugs, output_dir = results_dirs$dsen_results,
                               ncores = NCORES, overwrite = FALSE) {

  # Filter drugs that need processing
  if (!overwrite) {
    existing <- sapply(drugs, function(d) {
      file.exists(file.path(output_dir, paste0("dsen_", gsub("/", "_", d), ".Rda")))
    })
    drugs_to_run <- drugs[!existing]
    message(sprintf("%d drugs already processed, %d to run",
                    sum(existing), length(drugs_to_run)))
  } else {
    drugs_to_run <- drugs
  }

  if (length(drugs_to_run) == 0) {
    message("All drugs already processed.")
    return(NULL)
  }

  # Set up parallel backend
  cl <- setup_parallel(ncores)

  clusterExport(cl, c("data", "output_dir", "paths",
                      "VARIANCE_THRESHOLD", "CORRELATION_THRESHOLD",
                      "CV_NFOLDS", "BOOTSTRAP_NREPS", "ALPHA_GRID", "GLOBAL_SEED"),
                envir = environment())

  clusterEvalQ(cl, {
    source("R/00_config.R")
    source("R/01_data_utils.R")
    source("R/02_cv_utils.R")
    source("R/03_dsen_utils.R")
  })

  results <- foreach(drug = drugs_to_run,
                      .combine = rbind,
                      .errorhandling = "pass",
                      .packages = c("glmnet", "caret", "Matrix")) %dopar% {

    tryCatch({
      # Prepare data with tissue filtering
      drug_data <- prepare_drug_data_filtered(data, drug, min_per_tissue = 5)

      if (is.null(drug_data)) {
        return(data.frame(drug = drug, status = "skipped_insufficient_tissues",
                          mse_dsen = NA, stringsAsFactors = FALSE))
      }

      # Run DSEN CV
      dsen_result <- dsen_tune_alpha(drug_data$x, drug_data$y, drug_data$tissue,
                                      seed = GLOBAL_SEED)

      # Run DSEN bootstrap
      dsen_bs <- dsen_bootstrap(drug_data$x, drug_data$y, drug_data$tissue,
                                 alpha = dsen_result$alpha.min,
                                 lambda = dsen_result$lambda.min,
                                 seed = GLOBAL_SEED)

      # Save results
      filename <- paste0("dsen_", gsub("/", "_", drug), ".Rda")
      save(drug_data, dsen_result, dsen_bs,
           file = file.path(output_dir, filename))

      # Return summary
      data.frame(
        drug = drug,
        status = "completed",
        n_samples = drug_data$n_samples,
        n_tissues = drug_data$n_tissues,
        mse_dsen = dsen_result$best_cvm,
        alpha = dsen_result$alpha.min,
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      data.frame(drug = drug, status = paste("error:", e$message),
                 mse_dsen = NA, stringsAsFactors = FALSE)
    })
  }

  stop_parallel(cl)

  return(results)
}

#' -----------------------------------------------------------------------------
#' Progress Monitoring
#' -----------------------------------------------------------------------------

#' Check progress of parallel run
#'
#' @param output_dir Directory with results
#' @param pattern File pattern to match
#' @return Summary of completed files
check_progress <- function(output_dir, pattern = "cv_*.Rda") {
  files <- list.files(output_dir, pattern = glob2rx(pattern))

  message(sprintf("Completed: %d files in %s", length(files), output_dir))

  if (length(files) > 0) {
    # Get file modification times
    file_paths <- file.path(output_dir, files)
    mod_times <- file.info(file_paths)$mtime

    message(sprintf("Most recent: %s at %s",
                    files[which.max(mod_times)],
                    format(max(mod_times), "%Y-%m-%d %H:%M:%S")))
  }

  return(length(files))
}

#' Combine results from parallel runs
#'
#' @param output_dir Directory with results
#' @param pattern File pattern
#' @return Combined data frame
combine_parallel_results <- function(output_dir, pattern = "cv_*.Rda") {
  files <- list.files(output_dir, pattern = glob2rx(pattern), full.names = TRUE)

  if (length(files) == 0) {
    message("No result files found.")
    return(NULL)
  }

  message(sprintf("Combining %d result files...", length(files)))

  all_results <- list()

  for (f in files) {
    tryCatch({
      load(f)  # Loads comparison object

      drug_name <- gsub("^cv_|^dsen_|\\.Rda$", "", basename(f))

      if (exists("comparison")) {
        all_results[[drug_name]] <- as.data.frame(comparison)
        all_results[[drug_name]]$drug <- drug_name
      }
    }, error = function(e) {
      warning(sprintf("Failed to load %s: %s", f, e$message))
    })
  }

  if (length(all_results) == 0) {
    return(NULL)
  }

  combined <- do.call(rbind, all_results)
  rownames(combined) <- NULL

  return(combined)
}

message("Parallel processing utilities loaded.")
