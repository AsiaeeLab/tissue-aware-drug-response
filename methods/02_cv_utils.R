#' =============================================================================
#' Drug Response Prediction - Cross-Validation Utilities
#' =============================================================================
#'
#' This script provides cross-validation functions with CORRECT standard error
#' calculation. This fixes the SE vs SD bug in the original implementation.
#'
#' Key fix: cvsd = sd(cv_errors) / sqrt(K)  [SE, not SD]
#'
#' Usage:
#'   library(here)
#'   source(here("methods", "00_paths.R"))
#'   source(here("methods", "01_data_utils.R"))
#'   source(here("methods", "02_cv_utils.R"))
#'
#' =============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(caret)
})

#' -----------------------------------------------------------------------------
#' Custom Cross-Validation with Correct SE
#' -----------------------------------------------------------------------------

#' Cross-validation for elastic net with proper per-fold preprocessing
#'
#' This function performs CV with feature selection INSIDE each fold,
#' avoiding data leakage. It returns standard ERROR (not SD).
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param alpha Elastic net mixing parameter
#' @param nfolds Number of CV folds
#' @param lambda Optional lambda sequence
#' @param var_thresh Variance threshold for feature filtering
#' @param cor_thresh Correlation threshold for feature filtering
#' @param seed Random seed for fold creation
#' @return List with cvm (mean CV error), cvse (standard error), lambda, etc.
my_cv_glmnet <- function(x, y, alpha = 0.5, nfolds = CV_NFOLDS,
                          lambda = NULL, var_thresh = VARIANCE_THRESHOLD,
                          cor_thresh = CORRELATION_THRESHOLD,
                          seed = NULL,
                          foldid = NULL) {

  n <- nrow(x)

  # Create folds (allow externally provided fold IDs)
  if (!is.null(foldid)) {
    folds <- foldid
  } else {
    if (!is.null(seed)) set.seed(seed)
    folds <- createFolds(y, k = nfolds, list = FALSE)
  }

  # Storage for results
  cv_errors <- NULL
  internal_lambda <- lambda

  for (k in 1:nfolds) {
    # Split data
    train_idx <- folds != k
    test_idx <- folds == k

    x_train <- x[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    x_test <- x[test_idx, , drop = FALSE]
    y_test <- y[test_idx]

    # Preprocess TRAINING data only (feature selection + scaling)
    preproc <- preprocess_features(x_train, y_train,
                                    var_thresh = var_thresh,
                                    cor_thresh = cor_thresh,
                                    scale_features = TRUE)

    x_train_proc <- preproc$x

    # Skip if too few features
    if (ncol(x_train_proc) < 2) {
      warning(sprintf("Fold %d: Only %d features after filtering, skipping",
                      k, ncol(x_train_proc)))
      next
    }

    # Apply same preprocessing to TEST data
    x_test_proc <- apply_preprocessing(x_test, preproc)

    # Ensure same features in train and test
    common_features <- intersect(colnames(x_train_proc), colnames(x_test_proc))
    x_train_proc <- x_train_proc[, common_features, drop = FALSE]
    x_test_proc <- x_test_proc[, common_features, drop = FALSE]

    # Fit model on training data
    if (k == 1 && is.null(lambda)) {
      # First fold: let glmnet determine lambda sequence
      fit <- glmnet(x_train_proc, y_train, alpha = alpha, nlambda = 100)
      internal_lambda <- fit$lambda
      cv_errors <- matrix(NA, nrow = nfolds, ncol = length(internal_lambda))
    } else {
      fit <- glmnet(x_train_proc, y_train, alpha = alpha,
                    lambda = internal_lambda)
    }

    # Predict on test data
    y_pred <- predict(fit, newx = x_test_proc, s = internal_lambda)

    # Compute MSE for each lambda
    fold_mse <- colMeans((y_test - y_pred)^2)

    # Store results (handle potential lambda mismatch)
    if (length(fold_mse) == ncol(cv_errors)) {
      cv_errors[k, ] <- fold_mse
    } else {
      # Lambda sequences might differ slightly - interpolate if needed
      cv_errors[k, 1:min(length(fold_mse), ncol(cv_errors))] <-
        fold_mse[1:min(length(fold_mse), ncol(cv_errors))]
    }
  }

  # Remove any all-NA rows (failed folds)
  valid_folds <- rowSums(!is.na(cv_errors)) > 0
  cv_errors <- cv_errors[valid_folds, , drop = FALSE]
  n_valid_folds <- sum(valid_folds)

  if (n_valid_folds < 2) {
    stop("Fewer than 2 valid folds completed")
  }

  # Compute mean and STANDARD ERROR (not SD!)
  cvm <- colMeans(cv_errors, na.rm = TRUE)
  cvsd <- apply(cv_errors, 2, sd, na.rm = TRUE)
  cvse <- cvsd / sqrt(n_valid_folds)  # THIS IS THE KEY FIX

  # Find optimal lambda
  min_idx <- which.min(cvm)
  lambda_min <- internal_lambda[min_idx]

  # Find lambda.1se (largest lambda within 1 SE of minimum)
  threshold <- cvm[min_idx] + cvse[min_idx]
  lambda_1se_idx <- max(which(cvm <= threshold))
  lambda_1se <- internal_lambda[lambda_1se_idx]

  return(list(
    lambda = internal_lambda,
    cvm = cvm,
    cvse = cvse,              # Standard ERROR
    cvsd = cvsd,              # Also keep SD for reference
    nfolds = n_valid_folds,
    lambda.min = lambda_min,
    lambda.1se = lambda_1se,
    alpha = alpha,
    cv_errors = cv_errors     # Raw fold errors for diagnostics
  ))
}

#' Grid search over alpha values
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param alphas Vector of alpha values to try
#' @param nfolds Number of CV folds
#' @param seed Random seed
#' @return List with results for each alpha and optimal parameters
cv_tune_alpha <- function(x, y, alphas = ALPHA_GRID, nfolds = CV_NFOLDS,
                           seed = GLOBAL_SEED, foldid = NULL) {

  results <- list()
  min_errors <- numeric(length(alphas))

  for (i in seq_along(alphas)) {
    alpha <- alphas[i]
    message(sprintf("  Testing alpha = %.1f", alpha))

    cv_result <- my_cv_glmnet(x, y, alpha = alpha, nfolds = nfolds,
                              seed = seed, foldid = foldid)
    results[[i]] <- cv_result
    min_errors[i] <- min(cv_result$cvm, na.rm = TRUE)
  }

  # Find best alpha
  best_idx <- which.min(min_errors)
  best_alpha <- alphas[best_idx]
  best_result <- results[[best_idx]]

  return(list(
    alphas = alphas,
    all_results = results,
    min_errors = min_errors,
    alpha.min = best_alpha,
    lambda.min = best_result$lambda.min,
    lambda.1se = best_result$lambda.1se,
    best_cvm = min(best_result$cvm),
    best_cvse = best_result$cvse[which.min(best_result$cvm)]
  ))
}

#' -----------------------------------------------------------------------------
#' Bootstrap Feature Importance
#' -----------------------------------------------------------------------------

#' Bootstrap resampling for feature selection stability
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param alpha Optimal alpha from CV (default 0.5)
#' @param lambda Optimal lambda from CV (default NULL, will use CV)
#' @param nboot Number of bootstrap replicates
#' @param method "right" (per-bootstrap preprocessing) or "wrong" (global preprocessing)
#' @param var_thresh Variance threshold
#' @param cor_thresh Correlation threshold
#' @param seed Random seed
#' @return List with feature counts (freq) and selection frequencies
bootstrap_features <- function(x, y, alpha = 0.5, lambda = NULL,
                                nboot = BOOTSTRAP_NREPS,
                                method = "right",
                                var_thresh = VARIANCE_THRESHOLD,
                                cor_thresh = CORRELATION_THRESHOLD,
                                seed = GLOBAL_SEED) {

  set.seed(seed)
  n <- nrow(x)
  all_features <- colnames(x)

  # Storage
  feature_counts <- rep(0, length(all_features))
  names(feature_counts) <- all_features
  coef_sums <- rep(0, length(all_features))
  names(coef_sums) <- all_features

  # For "wrong" method, preprocess globally ONCE
  if (method == "wrong") {
    global_preproc <- preprocess_features(x, y,
                                           var_thresh = var_thresh,
                                           cor_thresh = cor_thresh,
                                           scale_features = TRUE)
    x_global <- global_preproc$x
  }

  for (b in 1:nboot) {
    # Bootstrap sample
    boot_idx <- sample(n, n, replace = TRUE)

    if (method == "wrong") {
      # Wrong: use globally preprocessed data
      x_boot <- x_global[boot_idx, , drop = FALSE]
      y_boot <- y[boot_idx]
      x_proc <- x_boot
    } else {
      # Right: preprocess within each bootstrap
      x_boot <- x[boot_idx, , drop = FALSE]
      y_boot <- y[boot_idx]

      preproc <- preprocess_features(x_boot, y_boot,
                                      var_thresh = var_thresh,
                                      cor_thresh = cor_thresh,
                                      scale_features = TRUE)

      if (ncol(preproc$x) < 2) next
      x_proc <- preproc$x
    }

    if (ncol(x_proc) < 2) next

    # Fit model with CV if lambda not provided
    tryCatch({
      if (is.null(lambda)) {
        cv_fit <- cv.glmnet(x_proc, y_boot, alpha = alpha, nfolds = 5)
        fit <- glmnet(x_proc, y_boot, alpha = alpha)
        use_lambda <- cv_fit$lambda.min
      } else {
        fit <- glmnet(x_proc, y_boot, alpha = alpha, lambda = lambda)
        use_lambda <- lambda
      }

      coefs <- as.vector(coef(fit, s = use_lambda))[-1]  # Remove intercept
      names(coefs) <- colnames(x_proc)

      # Update counts for features that passed preprocessing AND have non-zero coefs
      selected <- names(coefs)[coefs != 0]
      feature_counts[selected] <- feature_counts[selected] + 1
      coef_sums[names(coefs)] <- coef_sums[names(coefs)] + coefs

    }, error = function(e) {
      # Silent failure for bootstrap iteration
    })

    if (b %% 50 == 0) {
      message(sprintf("  Bootstrap %d/%d complete", b, nboot))
    }
  }

  # Compute selection frequency
  selection_freq <- feature_counts / nboot

  # Sort by selection frequency
  sorted_idx <- order(selection_freq, decreasing = TRUE)

  return(list(
    freq = selection_freq,
    cnt = feature_counts,
    mean_coef = coef_sums / nboot,
    sorted_features = names(selection_freq)[sorted_idx],
    nboot = nboot,
    method = method,
    top_features = names(selection_freq)[sorted_idx][1:min(20, length(sorted_idx))]
  ))
}

#' -----------------------------------------------------------------------------
#' Wrong CV (for comparison - replicates published methods)
#' -----------------------------------------------------------------------------

#' "Wrong" CV: Feature selection BEFORE splitting (data leakage)
#'
#' This function replicates the flawed methodology used in published papers.
#' Feature selection is done on the FULL dataset before CV.
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param alpha Elastic net mixing parameter
#' @param nfolds Number of CV folds
#' @param var_thresh Variance threshold
#' @param cor_thresh Correlation threshold
#' @param seed Random seed
#' @return List with CV results (biased due to leakage)
wrong_cv_glmnet <- function(x, y, alpha = 0.5, nfolds = CV_NFOLDS,
                             var_thresh = VARIANCE_THRESHOLD,
                             cor_thresh = CORRELATION_THRESHOLD,
                             seed = NULL) {

  # WRONG: Preprocess on FULL data before CV
  preproc <- preprocess_features(x, y,
                                  var_thresh = var_thresh,
                                  cor_thresh = cor_thresh,
                                  scale_features = TRUE)

  x_proc <- preproc$x

  if (ncol(x_proc) < 2) {
    stop("Too few features after filtering")
  }

  # Now run standard CV (but features already selected on full data = LEAKAGE)
  if (!is.null(seed)) set.seed(seed)

  cv_fit <- cv.glmnet(x_proc, y, alpha = alpha, nfolds = nfolds)

  # Note: cv.glmnet returns cvsd which is actually SE (divides by sqrt(nfolds))
  return(list(
    lambda = cv_fit$lambda,
    cvm = cv_fit$cvm,
    cvse = cv_fit$cvsd,  # glmnet's cvsd is actually SE
    nfolds = nfolds,
    lambda.min = cv_fit$lambda.min,
    lambda.1se = cv_fit$lambda.1se,
    alpha = alpha,
    n_features = ncol(x_proc),
    selected_features = colnames(x_proc)
  ))
}

#' Grid search with wrong CV (for comparison)
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param alphas Vector of alpha values
#' @param nfolds Number of CV folds
#' @param seed Random seed
#' @return List with results
wrong_cv_tune_alpha <- function(x, y, alphas = ALPHA_GRID, nfolds = CV_NFOLDS,
                                 seed = GLOBAL_SEED) {

  # WRONG: Preprocess on FULL data first
  preproc <- preprocess_features(x, y, scale_features = TRUE)
  x_proc <- preproc$x

  results <- list()
  min_errors <- numeric(length(alphas))

  for (i in seq_along(alphas)) {
    alpha <- alphas[i]
    set.seed(seed)

    cv_fit <- cv.glmnet(x_proc, y, alpha = alpha, nfolds = nfolds)

    results[[i]] <- list(
      lambda = cv_fit$lambda,
      cvm = cv_fit$cvm,
      cvse = cv_fit$cvsd,
      lambda.min = cv_fit$lambda.min,
      alpha = alpha
    )
    min_errors[i] <- min(cv_fit$cvm)
  }

  best_idx <- which.min(min_errors)

  return(list(
    alphas = alphas,
    all_results = results,
    min_errors = min_errors,
    alpha.min = alphas[best_idx],
    lambda.min = results[[best_idx]]$lambda.min,
    best_cvm = min_errors[best_idx],
    n_features = ncol(x_proc),
    selected_features = colnames(x_proc)
  ))
}

#' -----------------------------------------------------------------------------
#' Comparison Metrics
#' -----------------------------------------------------------------------------

#' Compute Jaccard similarity between two feature sets
#'
#' @param features1 Character vector of feature names
#' @param features2 Character vector of feature names
#' @return Numeric Jaccard index (0 to 1)
jaccard_similarity <- function(features1, features2) {
  intersection <- length(intersect(features1, features2))
  union <- length(union(features1, features2))

  if (union == 0) return(NA)
  return(intersection / union)
}

#' Compare wrong vs right CV results
#'
#' @param wrong_result Result from wrong_cv_tune_alpha()
#' @param right_result Result from cv_tune_alpha()
#' @param bs_wrong Bootstrap result from wrong CV
#' @param bs_right Bootstrap result from right CV
#' @param threshold Selection frequency threshold for Jaccard
#' @return List with comparison metrics
compare_cv_results <- function(wrong_result, right_result,
                                bs_wrong = NULL, bs_right = NULL,
                                threshold = 0.8) {

  # MSE comparison
  mse_wrong <- wrong_result$best_cvm
  mse_right <- right_result$best_cvm
  mse_diff <- mse_wrong - mse_right
  mse_pct_increase <- (mse_right - mse_wrong) / mse_wrong * 100

  # Feature count comparison
  n_feat_wrong <- if (!is.null(bs_wrong))
    sum(bs_wrong$freq >= threshold) else if (!is.null(wrong_result$n_features))
      wrong_result$n_features else NA
  n_feat_right <- if (!is.null(bs_right))
    sum(bs_right$freq >= threshold) else NA

  # Jaccard similarity (if bootstrap results available)
  jaccard <- NA
  if (!is.null(bs_wrong) && !is.null(bs_right)) {
    feat_wrong <- names(bs_wrong$freq)[bs_wrong$freq >= threshold]
    feat_right <- names(bs_right$freq)[bs_right$freq >= threshold]
    jaccard <- jaccard_similarity(feat_wrong, feat_right)
  }

  return(list(
    mse_wrong = mse_wrong,
    mse_right = mse_right,
    mse_diff = mse_diff,
    mse_pct_increase = mse_pct_increase,
    n_features_wrong = n_feat_wrong,
    n_features_right = n_feat_right,
    jaccard = jaccard,
    alpha_wrong = wrong_result$alpha.min,
    alpha_right = right_result$alpha.min,
    lambda_wrong = wrong_result$lambda.min,
    lambda_right = right_result$lambda.min
  ))
}

message("Cross-validation utilities loaded (with correct SE calculation).")
