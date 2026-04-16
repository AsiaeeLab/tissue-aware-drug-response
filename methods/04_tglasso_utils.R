# =============================================================================
# TG-LASSO Utilities
# =============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(caret)
})

# Ensure test matrix has all training features in the same order.
align_to_features <- function(x_mat, feature_names) {
  x_mat <- as_numeric_matrix(x_mat)
  if (length(feature_names) == 0) {
    return(matrix(0, nrow = nrow(x_mat), ncol = 0))
  }

  out <- matrix(0, nrow = nrow(x_mat), ncol = length(feature_names))
  colnames(out) <- feature_names
  rownames(out) <- rownames(x_mat)

  common <- intersect(colnames(x_mat), feature_names)
  if (length(common) > 0) {
    out[, common] <- x_mat[, common, drop = FALSE]
  }
  out
}

fit_tglasso_tissue_model <- function(x_train, y_train, tissue_train, target_tissue,
                                     var_thresh = VARIANCE_THRESHOLD,
                                     cor_thresh = CORRELATION_THRESHOLD,
                                     nlambda = 100L) {
  val_idx <- which(tissue_train == target_tissue)
  fit_idx <- which(tissue_train != target_tissue)

  if (length(val_idx) < 2 || length(fit_idx) < 10) {
    return(NULL)
  }

  x_fit <- x_train[fit_idx, , drop = FALSE]
  y_fit <- y_train[fit_idx]
  x_val <- x_train[val_idx, , drop = FALSE]
  y_val <- y_train[val_idx]

  pre <- preprocess_features(
    x_fit, y_fit,
    var_thresh = var_thresh,
    cor_thresh = cor_thresh,
    scale_features = TRUE
  )

  if (is.null(pre) || ncol(pre$x) < 2) return(NULL)

  x_fit_proc <- pre$x
  x_val_proc <- apply_preprocessing(x_val, pre)
  x_val_proc <- align_to_features(x_val_proc, colnames(x_fit_proc))

  fit <- glmnet(
    x_fit_proc, y_fit,
    alpha = 1,
    nlambda = as.integer(nlambda)
  )

  pred_val <- predict(fit, newx = x_val_proc, s = fit$lambda)
  val_mse <- colMeans((y_val - pred_val)^2)
  best_idx <- which.min(val_mse)
  best_lambda <- fit$lambda[best_idx]

  coef_vec <- as.vector(coef(fit, s = best_lambda))[-1]
  names(coef_vec) <- colnames(x_fit_proc)
  selected_features <- names(coef_vec)[coef_vec != 0]

  list(
    target_tissue = target_tissue,
    fit = fit,
    preproc = pre,
    feature_names = colnames(x_fit_proc),
    lambda = best_lambda,
    val_mse = val_mse[best_idx],
    selected_features = selected_features
  )
}

fit_tglasso_models_for_training <- function(x_train, y_train, tissue_train,
                                            var_thresh = VARIANCE_THRESHOLD,
                                            cor_thresh = CORRELATION_THRESHOLD,
                                            nlambda = 100L) {
  tissues <- sort(unique(tissue_train[!is.na(tissue_train)]))
  models <- list()

  for (tis in tissues) {
    model <- tryCatch(
      fit_tglasso_tissue_model(
        x_train, y_train, tissue_train,
        target_tissue = tis,
        var_thresh = var_thresh,
        cor_thresh = cor_thresh,
        nlambda = nlambda
      ),
      error = function(e) NULL
    )

    if (!is.null(model)) {
      models[[tis]] <- model
    }
  }

  # Fallback global model for tissues absent in training splits.
  fallback <- tryCatch({
    pre <- preprocess_features(
      x_train, y_train,
      var_thresh = var_thresh,
      cor_thresh = cor_thresh,
      scale_features = TRUE
    )
    if (is.null(pre) || ncol(pre$x) < 2) {
      NULL
    } else {
      fit <- glmnet(pre$x, y_train, alpha = 1, nlambda = as.integer(nlambda))
      pred_train <- predict(fit, newx = pre$x, s = fit$lambda)
      mse_train <- colMeans((y_train - pred_train)^2)
      lambda <- fit$lambda[which.min(mse_train)]

      coef_vec <- as.vector(coef(fit, s = lambda))[-1]
      names(coef_vec) <- colnames(pre$x)

      list(
        fit = fit,
        preproc = pre,
        feature_names = colnames(pre$x),
        lambda = lambda,
        selected_features = names(coef_vec)[coef_vec != 0]
      )
    }
  }, error = function(e) NULL)

  list(models = models, fallback = fallback)
}

predict_tglasso_model <- function(model_obj, x_new) {
  x_proc <- apply_preprocessing(x_new, model_obj$preproc)
  x_aligned <- align_to_features(x_proc, model_obj$feature_names)
  as.numeric(predict(model_obj$fit, newx = x_aligned, s = model_obj$lambda))
}

#' TG-LASSO cross-validation with optional externally supplied outer folds.
#'
#' @return list(cv_errors, cvm, predictions, foldid, lambda_by_fold_tissue,
#'              selected_features)
tglasso_cv <- function(x, y, tissue, nfolds = 5, foldid = NULL, seed = 42,
                       cor_threshold = 0.1, var_threshold = 0.01,
                       nlambda = 100L) {
  x <- as_numeric_matrix(x)
  x[is.na(x)] <- 0

  if (is.null(foldid)) {
    set.seed(seed)
    folds <- caret::createFolds(tissue, k = nfolds, list = FALSE)
  } else {
    folds <- as.integer(foldid)
    nfolds <- length(unique(folds))
  }

  predictions <- rep(NA_real_, length(y))
  cv_errors <- rep(NA_real_, nfolds)
  lambda_by_fold_tissue <- list()
  selected_features <- list()

  for (k in seq_len(nfolds)) {
    train_idx <- which(folds != k)
    test_idx <- which(folds == k)

    x_train <- x[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    tissue_train <- tissue[train_idx]

    x_test <- x[test_idx, , drop = FALSE]
    y_test <- y[test_idx]
    tissue_test <- tissue[test_idx]

    model_pack <- fit_tglasso_models_for_training(
      x_train, y_train, tissue_train,
      var_thresh = var_threshold,
      cor_thresh = cor_threshold,
      nlambda = nlambda
    )

    models <- model_pack$models
    fallback <- model_pack$fallback

    lambda_info <- lapply(models, function(m) m$lambda)
    lambda_by_fold_tissue[[as.character(k)]] <- lambda_info

    feat_info <- lapply(models, function(m) m$selected_features)
    selected_features[[as.character(k)]] <- feat_info

    y_pred_fold <- rep(NA_real_, length(test_idx))
    test_tissues <- unique(tissue_test)

    for (tis in test_tissues) {
      tis_idx <- which(tissue_test == tis)
      x_slice <- x_test[tis_idx, , drop = FALSE]

      if (tis %in% names(models)) {
        y_pred_fold[tis_idx] <- predict_tglasso_model(models[[tis]], x_slice)
      } else if (!is.null(fallback)) {
        y_pred_fold[tis_idx] <- predict_tglasso_model(fallback, x_slice)
      } else {
        y_pred_fold[tis_idx] <- mean(y_train, na.rm = TRUE)
      }
    }

    predictions[test_idx] <- y_pred_fold
    cv_errors[k] <- mean((y_test - y_pred_fold)^2, na.rm = TRUE)
  }

  list(
    cv_errors = cv_errors,
    cvm = mean(cv_errors, na.rm = TRUE),
    predictions = predictions,
    foldid = folds,
    lambda_by_fold_tissue = lambda_by_fold_tissue,
    selected_features = selected_features
  )
}

# Bootstrap TG-LASSO feature stability by within-tissue resampling.
tglasso_bootstrap_features <- function(x, y, tissue, nboot = 200,
                                       seed = 42,
                                       var_threshold = 0.01,
                                       cor_threshold = 0.1,
                                       nlambda = 100L) {
  x <- as_numeric_matrix(x)
  x[is.na(x)] <- 0
  features <- colnames(x)

  counts <- setNames(rep(0, length(features)), features)
  tissues <- sort(unique(tissue[!is.na(tissue)]))

  set.seed(seed)

  for (b in seq_len(nboot)) {
    boot_idx <- unlist(lapply(tissues, function(tis) {
      idx <- which(tissue == tis)
      sample(idx, length(idx), replace = TRUE)
    }), use.names = FALSE)

    xb <- x[boot_idx, , drop = FALSE]
    yb <- y[boot_idx]
    tb <- tissue[boot_idx]

    model_pack <- tryCatch(
      fit_tglasso_models_for_training(
        xb, yb, tb,
        var_thresh = var_threshold,
        cor_thresh = cor_threshold,
        nlambda = nlambda
      ),
      error = function(e) NULL
    )

    if (is.null(model_pack)) next

    selected <- unique(unlist(lapply(model_pack$models, function(m) m$selected_features),
                              use.names = FALSE))
    if (!is.null(model_pack$fallback)) {
      selected <- unique(c(selected, model_pack$fallback$selected_features))
    }

    if (length(selected) > 0) {
      selected <- intersect(selected, names(counts))
      counts[selected] <- counts[selected] + 1
    }

    if (b %% 50 == 0) {
      message(sprintf("  TG-LASSO bootstrap %d/%d", b, nboot))
    }
  }

  list(
    cnt = counts,
    freq = counts / nboot,
    nboot = nboot,
    sorted_features = names(sort(counts, decreasing = TRUE))
  )
}

message("TG-LASSO utilities loaded.")
