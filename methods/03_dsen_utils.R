#' =============================================================================
#' Drug Response Prediction - Data Shared Elastic Net (DSEN) Utilities
#' =============================================================================
#'
#' This script provides functions for Data Shared Elastic Net, which models
#' drug response with tissue-specific and shared coefficients.
#'
#' Based on: Gross & Tibshirani (2016) "Data Shared Lasso"
#'
#' Key fixes in this implementation:
#' 1. Per-fold scaling (no global scaling before CV)
#' 2. Bootstrap within tissue type
#' 3. Correct standard error calculation
#'
#' Usage:
#'   library(here)
#'   source(here("methods", "00_paths.R"))
#'   source(here("methods", "01_data_utils.R"))
#'   source(here("methods", "02_cv_utils.R"))
#'   source(here("methods", "03_dsen_utils.R"))
#'
#' =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(glmnet)
  library(caret)
})

#' -----------------------------------------------------------------------------
#' DSEN Matrix Construction
#' -----------------------------------------------------------------------------

#' Build block-diagonal design matrix for DSEN
#'
#' The DSEN formulation creates a matrix Z where:
#' - First p columns: shared coefficients (same across tissues)
#' - Next p*T columns: tissue-specific deviations (block diagonal)
#'
#' @param data_list List of lists, each with $x (sparse matrix) and $y (vector)
#' @return Sparse matrix Z
build_dsen_matrix <- function(data_list) {
  n_tasks <- length(data_list)
  if (n_tasks == 0) stop("data_list is empty")

  p <- ncol(data_list[[1]]$x)

  # Get sample sizes per task
  ns <- sapply(data_list, function(item) nrow(item$x))
  n_total <- sum(ns)

  # Cumulative sample counts for indexing (1-based)
  start_rows <- c(1, cumsum(ns)[-n_tasks] + 1)

  # Build sparse matrix indices
  is <- integer(0)
  js <- integer(0)
  xs <- numeric(0)

  for (t in seq_len(n_tasks)) {
    x_t <- data_list[[t]]$x

    if (!inherits(x_t, "sparseMatrix")) {
      x_t <- as.matrix(x_t)
    }

    nz_idx <- which(x_t != 0, arr.ind = TRUE)
    if (nrow(nz_idx) == 0) next

    row_offset <- start_rows[t] - 1
    i_vals <- row_offset + nz_idx[, 1]
    j_vals <- nz_idx[, 2]
    vals <- as.numeric(x_t[nz_idx])

    # Shared coefficients (first p columns)
    is <- c(is, i_vals, i_vals)
    js <- c(js, j_vals, p * t + j_vals)
    xs <- c(xs, vals, vals)
  }

  n_cols <- p * (n_tasks + 1)
  if (length(is) == 0) {
    return(sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                        dims = c(n_total, n_cols)))
  }

  # Guard against out-of-bounds indices
  valid <- is >= 1 & is <= n_total & js >= 1 & js <= n_cols
  if (!all(valid)) {
    is <- is[valid]
    js <- js[valid]
    xs <- xs[valid]
  }

  Z <- sparseMatrix(i = is, j = js, x = xs, dims = c(n_total, n_cols))
  return(Z)
}

#' Get column index from sparse matrix
#' (Helper for efficient sparse matrix manipulation)
getj <- function(mat) {
  dp <- diff(mat@p)
  rep(seq_along(dp), dp) - 1
}

#' Compute penalty factors for DSEN
#'
#' Creates penalty factor vector for glmnet where:
#' - Shared coefficients get penalty = lamrat (controls shared vs tissue-specific tradeoff)
#' - Tissue-specific coefficients get penalty based on tissue sample size
#'   (smaller tissues get higher penalty to prevent overfitting)
#'
#' @param data_list List of tissue data (each with $x matrix)
#' @param lamrat Lambda ratio for shared coefficients
#' @param tissue_penalty_mode How to compute tissue penalties:
#'   - "uniform": All tissues get penalty = 1 (original behavior)
#'   - "sample_size": Penalty = sqrt(max_n / n_tissue) - smaller tissues penalized more
#'   - "sample_size_mild": Penalty = (max_n / n_tissue)^0.25 - milder adjustment
#' @return Numeric vector of penalty factors
compute_dsen_penalties <- function(data_list, lamrat,
                                    tissue_penalty_mode = c("sample_size", "uniform", "sample_size_mild")) {
  tissue_penalty_mode <- match.arg(tissue_penalty_mode)

  n_tasks <- length(data_list)
  p <- ncol(data_list[[1]]$x)

 # Get sample sizes per tissue
  tissue_sizes <- sapply(data_list, function(item) nrow(item$x))
  max_n <- max(tissue_sizes)

  # Compute tissue-specific penalties
  if (tissue_penalty_mode == "uniform") {
    # Original behavior: all tissues get penalty = 1
    tissue_penalties <- rep(1, n_tasks)
  } else if (tissue_penalty_mode == "sample_size") {
    # Penalty inversely proportional to sqrt(sample size)
    # Larger tissues can estimate more reliably -> lower penalty
    tissue_penalties <- sqrt(max_n / tissue_sizes)
  } else if (tissue_penalty_mode == "sample_size_mild") {
    # Milder adjustment: fourth root instead of square root
    tissue_penalties <- (max_n / tissue_sizes)^0.25
  }

  # Build penalty factor vector:
  # [lamrat repeated p times for shared, then tissue penalties repeated p times each]
  pf <- c(
    rep(lamrat, p),                           # Shared block
    rep(tissue_penalties, each = p)           # Tissue-specific blocks
  )

  return(pf)
}

#' -----------------------------------------------------------------------------
#' DSEN Cross-Validation (with per-fold scaling)
#' -----------------------------------------------------------------------------

#' Prepare data list by tissue with per-fold preprocessing
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param tissue Tissue labels
#' @param train_idx Training sample indices
#' @param var_thresh Variance threshold
#' @param cor_thresh Correlation threshold
#' @return List of tissue-specific data (sparse matrices)
prepare_tissue_data <- function(x, y, tissue, train_idx,
                                 var_thresh = VARIANCE_THRESHOLD,
                                 cor_thresh = CORRELATION_THRESHOLD,
                                 screening = c("pooled", "union_threshold", "union_topk"),
                                 topk = 200L) {

  screening <- match.arg(screening)

  # Subset to training data
  x_train <- x[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  tissue_train <- tissue[train_idx]

  # Preprocess on TRAINING data only
  x_train <- as_numeric_matrix(x_train)
  x_train[is.na(x_train)] <- 0

  # ---------------------------------------------------------------------------
  # Feature screening (DSEN requires same features across tissues/tasks).
  # - pooled: pooled variance + pooled correlation thresholding (current default)
  # - union_threshold: union over tissues of (var>thresh & |cor|>thresh)
  # - union_topk: union of top-k |cor| features per tissue (after pooled variance)
  # ---------------------------------------------------------------------------

  # Pooled variance filter (always apply first to reduce dimensionality)
  pooled_vars <- apply(x_train, 2, var, na.rm = TRUE)
  keep_var_pooled <- pooled_vars > var_thresh
  x_train <- x_train[, keep_var_pooled, drop = FALSE]

  tissues <- unique(tissue_train[!is.na(tissue_train)])

  if (screening == "pooled") {
    pooled_cors <- suppressWarnings(cor(x_train, y_train, use = "pairwise.complete.obs"))
    pooled_cors[is.na(pooled_cors)] <- 0
    keep_cor <- abs(pooled_cors) > cor_thresh
    x_train <- x_train[, keep_cor, drop = FALSE]
  } else if (screening == "union_threshold") {
    keep_any <- rep(FALSE, ncol(x_train))
    for (tis in tissues) {
      tis_idx <- which(tissue_train == tis)
      if (length(tis_idx) < 3) next

      x_tis <- x_train[tis_idx, , drop = FALSE]
      y_tis <- y_train[tis_idx]

      tis_vars <- apply(x_tis, 2, var, na.rm = TRUE)
      keep_var <- tis_vars > var_thresh

      tis_cors <- suppressWarnings(cor(x_tis, y_tis, use = "pairwise.complete.obs"))
      tis_cors[is.na(tis_cors)] <- 0
      keep_cor <- abs(tis_cors) > cor_thresh

      keep_any <- keep_any | (keep_var & keep_cor)
    }

    if (any(keep_any)) {
      x_train <- x_train[, keep_any, drop = FALSE]
    } else {
      pooled_cors <- suppressWarnings(cor(x_train, y_train, use = "pairwise.complete.obs"))
      pooled_cors[is.na(pooled_cors)] <- 0
      keep_cor <- abs(pooled_cors) > cor_thresh
      x_train <- x_train[, keep_cor, drop = FALSE]
    }
  } else if (screening == "union_topk") {
    topk <- as.integer(topk)
    if (is.na(topk) || topk < 1) topk <- 200L

    keep_any <- rep(FALSE, ncol(x_train))
    for (tis in tissues) {
      tis_idx <- which(tissue_train == tis)
      if (length(tis_idx) < 3) next

      x_tis <- x_train[tis_idx, , drop = FALSE]
      y_tis <- y_train[tis_idx]

      tis_cors <- suppressWarnings(cor(x_tis, y_tis, use = "pairwise.complete.obs"))
      tis_abs <- abs(tis_cors)
      tis_abs[is.na(tis_abs)] <- 0

      if (all(tis_abs == 0)) next

      ord <- order(tis_abs, decreasing = TRUE)
      take <- ord[seq_len(min(topk, length(ord)))]
      keep_any[take] <- TRUE
    }

    if (any(keep_any)) {
      x_train <- x_train[, keep_any, drop = FALSE]
    }
  }

  # Remove duplicates
  x_train <- x_train[, !duplicated(t(x_train)), drop = FALSE]

  # Scale (excluding mutations)
  is_mutation <- grepl("_mut(\\.|$)|^mut_", colnames(x_train), ignore.case = TRUE)
  scale_center <- NULL
  scale_scale <- NULL

  if (any(!is_mutation)) {
    x_to_scale <- x_train[, !is_mutation, drop = FALSE]
    scale_center <- colMeans(x_to_scale, na.rm = TRUE)
    scale_scale <- apply(x_to_scale, 2, sd, na.rm = TRUE)
    scale_scale[is.na(scale_scale) | scale_scale == 0] <- 1

    x_train[, !is_mutation] <- sweep(
      sweep(x_to_scale, 2, scale_center, "-"),
      2, scale_scale, "/"
    )
  }

  selected_features <- colnames(x_train)

  # Build tissue-specific data list
  data_list <- list()
  tissues <- unique(tissue_train[!is.na(tissue_train)])

  for (tis in tissues) {
    tis_idx <- which(tissue_train == tis)

    if (length(tis_idx) < 2) next

    y_tis <- y_train[tis_idx]
    x_tis <- x_train[tis_idx, , drop = FALSE]

    # Convert to sparse matrix
    nz_idx <- which(x_tis != 0, arr.ind = TRUE)
    if (nrow(nz_idx) == 0) {
      x_sparse <- sparseMatrix(i = 1, j = 1, x = 0,
                               dims = c(nrow(x_tis), ncol(x_tis)))
    } else {
      x_sparse <- sparseMatrix(i = nz_idx[, 1], j = nz_idx[, 2],
                               x = x_tis[nz_idx],
                               dims = c(nrow(x_tis), ncol(x_tis)))
    }

    data_list[[tis]] <- list(y = y_tis, x = x_sparse)
  }

  return(list(
    data_list = data_list,
    selected_features = selected_features,
    scale_center = scale_center,
    scale_scale = scale_scale
  ))
}

#' DSEN Cross-Validation with proper per-fold preprocessing
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param tissue Tissue labels
#' @param alpha Elastic net mixing parameter
#' @param nfolds Number of CV folds
#' @param lamrat Lambda ratio (default: sqrt(n_tissues))
#' @param tissue_penalty_mode How to penalize tissue-specific coefficients:
#'   "sample_size" (default) - smaller tissues get higher penalty
#'   "uniform" - all tissues get same penalty (original behavior)
#'   "sample_size_mild" - milder sample-size adjustment
#' @param seed Random seed
#' @return List with CV results
dsen_cv <- function(x, y, tissue, alpha = 0.5, nfolds = CV_NFOLDS,
                     lamrat = NULL, seed = GLOBAL_SEED,
                     foldid = NULL,
                     tissue_penalty_mode = c("sample_size", "uniform", "sample_size_mild"),
                     screening = c("pooled", "union_threshold", "union_topk"),
                     topk = 200L,
                     nlambda = 100L) {

  screening <- match.arg(screening)
  tissue_penalty_mode <- match.arg(tissue_penalty_mode)
  nlambda <- as.integer(nlambda)
  if (is.na(nlambda) || nlambda < 10) nlambda <- 100L

  set.seed(seed)

  # Create stratified folds by tissue (or use provided fold assignment)
  if (is.null(foldid)) {
    folds <- createFolds(tissue, k = nfolds, list = FALSE)
  } else {
    folds <- foldid
  }

  cv_errors <- NULL
  internal_lambda <- NULL

  for (k in 1:nfolds) {
    train_idx <- which(folds != k)
    test_idx <- which(folds == k)

    # Prepare training data with per-fold preprocessing
    train_data <- prepare_tissue_data(
      x, y, tissue, train_idx,
      screening = screening,
      topk = topk
    )
    data_list <- train_data$data_list

    if (length(data_list) < 2) {
      warning(sprintf("Fold %d: fewer than 2 tissues, skipping", k))
      next
    }

    n_tasks <- length(data_list)
    p <- ncol(data_list[[1]]$x)

    # Set lambda ratio if not provided
    if (is.null(lamrat)) {
      lamrat_use <- sqrt(n_tasks)
    } else {
      lamrat_use <- lamrat
    }

    # Build DSEN matrix
    Z <- build_dsen_matrix(data_list)
    yz <- unlist(lapply(data_list, function(item) item$y))

    # Set penalty factors (tissue-specific penalties based on sample size)
    pf <- compute_dsen_penalties(data_list, lamrat_use, tissue_penalty_mode)

    # Fit DSEN
    if (k == 1) {
      fit <- glmnet(Z, yz, alpha = alpha, nlambda = nlambda,
                    penalty.factor = pf, intercept = TRUE,
                    standardize = FALSE)
      internal_lambda <- fit$lambda
      cv_errors <- matrix(NA, nrow = nfolds, ncol = length(internal_lambda))
    } else {
      fit <- glmnet(Z, yz, alpha = alpha, lambda = internal_lambda,
                    penalty.factor = pf, intercept = TRUE,
                    standardize = FALSE)
    }

    # Prepare test data
    x_test <- as_numeric_matrix(x[test_idx, train_data$selected_features, drop = FALSE])
    x_test[is.na(x_test)] <- 0

    # Apply scaling from training
    is_mutation <- grepl("_mut(\\.|$)|^mut_", colnames(x_test), ignore.case = TRUE)
    if (!is.null(train_data$scale_center)) {
      for (feat in names(train_data$scale_center)) {
        if (feat %in% colnames(x_test) && !is_mutation[colnames(x_test) == feat]) {
          x_test[, feat] <- (x_test[, feat] - train_data$scale_center[feat]) /
                             train_data$scale_scale[feat]
        }
      }
    }

    y_test <- y[test_idx]
    tissue_test <- tissue[test_idx]

    # Build test data list
    test_data_list <- list()
    for (tis in names(data_list)) {
      tis_idx <- which(tissue_test == tis)
      if (length(tis_idx) == 0) next

      y_tis <- y_test[tis_idx]
      x_tis <- x_test[tis_idx, , drop = FALSE]

      nz_idx <- which(x_tis != 0, arr.ind = TRUE)
      if (nrow(nz_idx) == 0) {
        x_sparse <- sparseMatrix(i = 1, j = 1, x = 0,
                                 dims = c(nrow(x_tis), ncol(x_tis)))
      } else {
        x_sparse <- sparseMatrix(i = nz_idx[, 1], j = nz_idx[, 2],
                                 x = x_tis[nz_idx],
                                 dims = c(nrow(x_tis), ncol(x_tis)))
      }

      test_data_list[[tis]] <- list(y = y_tis, x = x_sparse)
    }

    # Predict and compute error
    for (l_idx in seq_along(internal_lambda)) {
      coefs <- as.vector(coef(fit, s = internal_lambda[l_idx]))
      intercept <- coefs[1]
      beta <- coefs[-1]

      # Predict using DSEN structure
      y_pred <- dsen_predict(test_data_list, beta, p) + intercept
      y_true <- unlist(lapply(test_data_list, function(item) item$y))

      cv_errors[k, l_idx] <- mean((y_true - y_pred)^2)
    }
  }

  # Remove failed folds
  valid_folds <- rowSums(!is.na(cv_errors)) > 0
  cv_errors <- cv_errors[valid_folds, , drop = FALSE]
  n_valid_folds <- sum(valid_folds)

  # Compute mean and SE
  cvm <- colMeans(cv_errors, na.rm = TRUE)
  cvsd <- apply(cv_errors, 2, sd, na.rm = TRUE)
  cvse <- cvsd / sqrt(n_valid_folds)

  # Find optimal lambda
  min_idx <- which.min(cvm)

  return(list(
    lambda = internal_lambda,
    cvm = cvm,
    cvse = cvse,
    cvsd = cvsd,
    nfolds = n_valid_folds,
    lambda.min = internal_lambda[min_idx],
    alpha = alpha,
    cv_errors = cv_errors
  ))
}

#' DSEN prediction from coefficients
#'
#' @param data_list List of tissue data
#' @param beta Coefficient vector (length p * (n_tasks + 1))
#' @param p Number of features
#' @return Predicted values
dsen_predict <- function(data_list, beta, p) {
  n_tasks <- length(data_list)
  beta_shared <- beta[1:p]

  predictions <- c()

  for (t in seq_along(data_list)) {
    x_t <- data_list[[t]]$x
    beta_t <- beta[(p * t + 1):(p * (t + 1))]

    # Total coefficient = shared + tissue-specific
    pred_t <- as.numeric(x_t %*% beta_shared + x_t %*% beta_t)
    predictions <- c(predictions, pred_t)
  }

  return(predictions)
}

#' Grid search over alpha for DSEN
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param tissue Tissue labels
#' @param alphas Vector of alpha values
#' @param nfolds Number of folds
#' @param tissue_penalty_mode How to penalize tissue-specific coefficients
#' @param seed Random seed
#' @return List with optimal parameters
dsen_tune_alpha <- function(x, y, tissue, alphas = ALPHA_GRID,
                             nfolds = CV_NFOLDS,
                             tissue_penalty_mode = c("sample_size", "uniform", "sample_size_mild"),
                             seed = GLOBAL_SEED, foldid = NULL) {

  tissue_penalty_mode <- match.arg(tissue_penalty_mode)
  results <- list()
  min_errors <- numeric(length(alphas))

  for (i in seq_along(alphas)) {
    message(sprintf("  DSEN alpha = %.1f", alphas[i]))
    cv_result <- dsen_cv(x, y, tissue, alpha = alphas[i],
                          nfolds = nfolds, seed = seed, foldid = foldid,
                          tissue_penalty_mode = tissue_penalty_mode)
    results[[i]] <- cv_result
    min_errors[i] <- min(cv_result$cvm, na.rm = TRUE)
  }

  best_idx <- which.min(min_errors)

  return(list(
    alphas = alphas,
    all_results = results,
    min_errors = min_errors,
    alpha.min = alphas[best_idx],
    lambda.min = results[[best_idx]]$lambda.min,
    best_cvm = min_errors[best_idx],
    best_cvse = results[[best_idx]]$cvse[which.min(results[[best_idx]]$cvm)]
  ))
}

#' Grid search over alpha and lamrat for DSEN
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param tissue Tissue labels
#' @param alphas Vector of alpha values
#' @param lamrats Vector of shared-penalty multipliers
#' @param nfolds Number of folds
#' @param tissue_penalty_mode How to penalize tissue-specific coefficients
#' @param screening Feature screening mode (see prepare_tissue_data)
#' @param topk If screening == "union_topk", number of features per tissue
#' @param nlambda Number of lambdas in glmnet path
#' @param seed Random seed
#' @param foldid Optional fold assignments to reuse across models
#' @return List with optimal parameters
dsen_tune_alpha_lamrat <- function(x, y, tissue,
                                  alphas = ALPHA_GRID,
                                  lamrats = c(0.5, 1, 2, 4, 8, 16),
                                  nfolds = CV_NFOLDS,
                                  tissue_penalty_mode = c("sample_size", "uniform", "sample_size_mild"),
                                  screening = c("pooled", "union_threshold", "union_topk"),
                                  topk = 200L,
                                  nlambda = 100L,
                                  seed = GLOBAL_SEED,
                                  foldid = NULL) {

  screening <- match.arg(screening)
  tissue_penalty_mode <- match.arg(tissue_penalty_mode)
  lamrats <- as.numeric(lamrats)
  lamrats <- lamrats[is.finite(lamrats) & lamrats > 0]
  if (length(lamrats) == 0) stop("lamrats must contain positive finite values")

  grid <- expand.grid(
    alpha = as.numeric(alphas),
    lamrat = lamrats,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$best_cvm <- NA_real_
  grid$lambda_min <- NA_real_
  grid$best_cvse <- NA_real_

  results <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    a <- grid$alpha[i]
    r <- grid$lamrat[i]
    message(sprintf("  DSEN alpha = %.2f, lamrat = %.2f (%s)", a, r, screening))

    cv_result <- dsen_cv(
      x, y, tissue,
      alpha = a,
      nfolds = nfolds,
      lamrat = r,
      seed = seed,
      foldid = foldid,
      tissue_penalty_mode = tissue_penalty_mode,
      screening = screening,
      topk = topk,
      nlambda = nlambda
    )

    results[[i]] <- cv_result
    grid$best_cvm[i] <- min(cv_result$cvm, na.rm = TRUE)
    min_idx <- which.min(cv_result$cvm)
    grid$lambda_min[i] <- cv_result$lambda[min_idx]
    grid$best_cvse[i] <- cv_result$cvse[min_idx]
  }

  best_idx <- which.min(grid$best_cvm)
  best_cv <- results[[best_idx]]
  best_min_idx <- which.min(best_cv$cvm)

  list(
    grid = grid,
    all_results = results,
    alpha.min = grid$alpha[best_idx],
    lamrat.min = grid$lamrat[best_idx],
    lambda.min = best_cv$lambda[best_min_idx],
    best_cvm = grid$best_cvm[best_idx],
    best_cvse = best_cv$cvse[best_min_idx],
    screening = screening,
    topk = topk,
    nlambda = nlambda
  )
}

#' -----------------------------------------------------------------------------
#' DSEN Bootstrap (within-tissue resampling)
#' -----------------------------------------------------------------------------

#' Bootstrap for DSEN with within-tissue resampling
#'
#' @param x Feature matrix
#' @param y Response vector
#' @param tissue Tissue labels
#' @param alpha Optimal alpha
#' @param lambda Optimal lambda
#' @param lamrat Lambda ratio for shared coefficients (default: sqrt(n_tissues))
#' @param nboot Number of bootstrap replicates
#' @param tissue_penalty_mode How to penalize tissue-specific coefficients
#' @param seed Random seed
#' @return List with bootstrap results
dsen_bootstrap <- function(x, y, tissue, alpha, lambda,
                            lamrat = NULL,
                            nboot = BOOTSTRAP_NREPS,
                            tissue_penalty_mode = "sample_size",
                            seed = GLOBAL_SEED) {

  set.seed(seed)

  tissues <- unique(tissue[!is.na(tissue)])
  all_features <- colnames(x)

  # Storage for coefficient estimates
  # Columns: shared + one per tissue
  coef_matrix <- matrix(0, nrow = length(all_features),
                        ncol = length(tissues) + 1)
  rownames(coef_matrix) <- all_features
  colnames(coef_matrix) <- c("Shared", tissues)

  feature_counts <- rep(0, length(all_features))
  names(feature_counts) <- all_features

  # Tissue-specific selection counts (features x tissues)
  tissue_counts <- matrix(0, nrow = length(all_features), ncol = length(tissues))
  rownames(tissue_counts) <- all_features
  colnames(tissue_counts) <- tissues

  for (b in 1:nboot) {
    # Bootstrap WITHIN each tissue
    boot_idx <- c()
    for (tis in tissues) {
      tis_idx <- which(tissue == tis)
      if (length(tis_idx) < 2) next

      # Resample within this tissue
      boot_tis_idx <- sample(tis_idx, length(tis_idx), replace = TRUE)
      boot_idx <- c(boot_idx, boot_tis_idx)
    }

    # Prepare data
    train_data <- prepare_tissue_data(x, y, tissue, boot_idx)
    data_list <- train_data$data_list

    if (length(data_list) < 2) next

    n_tasks <- length(data_list)
    p <- ncol(data_list[[1]]$x)

    # Build DSEN matrix
    Z <- build_dsen_matrix(data_list)
    yz <- unlist(lapply(data_list, function(item) item$y))

    # Use provided lamrat or default to sqrt(n_tasks)
    lamrat_use <- if (is.null(lamrat)) sqrt(n_tasks) else lamrat
    pf <- compute_dsen_penalties(data_list, lamrat_use, tissue_penalty_mode)

    # Fit model
    tryCatch({
      fit <- glmnet(Z, yz, alpha = alpha, lambda = lambda,
                    penalty.factor = pf, intercept = TRUE,
                    standardize = FALSE)

      coefs <- as.vector(coef(fit, s = lambda))[-1]

      # Extract shared and tissue-specific coefficients
      beta_shared <- coefs[1:p]
      names(beta_shared) <- train_data$selected_features

      # Update shared coefficients
      coef_matrix[names(beta_shared), "Shared"] <-
        coef_matrix[names(beta_shared), "Shared"] + beta_shared

      # Update tissue-specific coefficients
      task_names <- names(data_list)
      for (t in seq_along(task_names)) {
        beta_t <- coefs[(p * t + 1):(p * (t + 1))]
        names(beta_t) <- train_data$selected_features
        tis_name <- task_names[t]

        if (tis_name %in% colnames(coef_matrix)) {
          coef_matrix[names(beta_t), tis_name] <-
            coef_matrix[names(beta_t), tis_name] + beta_t
        }
      }

      # Count selected features (shared)
      selected <- names(beta_shared)[beta_shared != 0]
      feature_counts[selected] <- feature_counts[selected] + 1

      # Count selected features (tissue-specific)
      for (t in seq_along(task_names)) {
        beta_t <- coefs[(p * t + 1):(p * (t + 1))]
        tis_name <- task_names[t]
        selected_t <- train_data$selected_features[beta_t != 0]
        if (tis_name %in% colnames(tissue_counts)) {
          tissue_counts[selected_t, tis_name] <-
            tissue_counts[selected_t, tis_name] + 1
        }
      }

    }, error = function(e) {
      warning(sprintf("Bootstrap %d failed: %s", b, e$message))
    })

    if (b %% 50 == 0) {
      message(sprintf("  DSEN Bootstrap %d/%d complete", b, nboot))
    }
  }

  # Average coefficients
  coef_matrix <- coef_matrix / nboot

  return(list(
    mean_coefs = coef_matrix,
    cnt = feature_counts,  # Raw counts (0-nboot) for shared block
    tissue_cnt = tissue_counts,  # Raw counts (0-nboot) per tissue, features x tissues
    nboot = nboot,
    tissues = tissues
  ))
}

message("DSEN utilities loaded (with per-fold scaling and within-tissue bootstrap).")
