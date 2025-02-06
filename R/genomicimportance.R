# Genomic Importance Functions

###########################
### Functions
###########################
#
# FitGIModelsCNtoS
# Jitter
# FitGIModelsStoCN
# MakeGIVectors
#

#' Fit a model to each bin in a CN-dataset predicting signature exposure
#'
#' @description
#' Fit a model for each signature in Y given each bin of copy-numbers in CN dataset X. \cr
#' The two model options are a linear model and a beta model.
#'
#' @param X *S4* A QDNAseq object that contains \code{copynumber}, and/or \code{segmented} assay slots. \cr
#' OR alternatively a matrix of bins (x) by samples (y) may be provided. The values being the copy-numbers.
#' @param Y *matrix* A matrix of samples (x) by signatures (y). Values are signature exposures. \cr
#' `dim(Y)[1]` must equal `dim(X)[2]`
#' @param model_option (optional) *character.* Which model to use. Options are \code{lm} or \code{beta}. \cr
#' For more modelling details see the `betareg` or `stats` R packages.
#' @param bin_indices (optional) *numeric/integer vector.* A single integer or a vector of integers. \cr
#' Subsets the copy-number object to only run on the bins in these indices.
#' @param data_slot (optional) *character.* If passing a QDNAseq object, the CN data slot to use. \cr
#' Common options are "copynumber" (corrected or uncorrected read counts), "segmented", or "calls".
#' @param run_parallel (optional) *logical.* If TRUE run job in parallel. Will automatically use 1 less than total cores available.
#' @param cores (optional) *integer.* Alternatively/in-addition, the user may choose the number of cores to use. \cr
#' Must be <= 1 less than the number available.
#' @param epsilon (optional) *numeric.* Applies to the beta-modelling option only. \cr
#' For response values too near-to 0 or 1, and or subtract this amount +/- a small jitter.
#'
#' @return A list of lists. One entry for each bin. Each bin then has a list of: \cr
#' `model_predvals` \cr
#' `predictor_index`
#'
#' @export
FitGIModelsCNtoS <- function (X, Y,
                              model_option = 'beta',
                              bin_indices = NULL,
                              data_slot = NULL,
                              run_parallel = FALSE,
                              cores = 1,
                              epsilon = 5e-4) {

  # Validate inputs
  stopifnot(inherits(Y, 'matrix'))
  stopifnot(dim(Y)[1] == dim(X)[2])
  if (inherits(X, c("QDNAseqCopyNumbers", "QDNAseqReadCounts"))) {
    if (!is.null(data_slot)) {
      stopifnot(data_slot %in% Biobase::assayDataElementNames(X))
    }
    X <- Biobase::assayDataElement(X, data_slot)
  } else {
    stopifnot(inherits(X, 'matrix'))
  }
  if (cores != 1) {
    require(foreach)
    available_cores <- parallel::detectCores() - 1
    stopifnot(available_cores >= cores && cores > 0)
    doMC::registerDoMC(cores)
  } else if (run_parallel) {
    available_cores <- parallel::detectCores() - 1
    doMC::registerDoMC(available_cores)
  }
  if (!is.null(bin_indices)) {
    stopifnot(all(bin_indices <= dim(X)[1]))
  } else {
    bin_indices <- 1:nrow(X)
  }

  # Run modelling
  start.time <- Sys.time()
  if (model_option == 'lm') {
    message("Begin per-bin LM runs...")
    results <- foreach(i = bin_indices, .packages = "stats") %dopar% {
      predictor <- X[i, ]
      na_mask <- !is.na(predictor)
      predictor <- predictor[na_mask]
      if (length(predictor) < dim(Y)[2]) {
        return(list(model_predvals = list(), predictor_index = i))
      } else {
        Y <- Y[na_mask,]
      }
      models <- list()
      for (j in seq_len(ncol(Y))) {
        response <- Y[, j]
        model <- lm(response ~ predictor)
        model_sum <- summary(model)
        models[[colnames(Y)[j]]] <- c(model_sum$coefficients[,1][2],
                                      model_sum$coefficients[,4][2])
      }
      return(list(model_predvals = models,
                  predictor_index = i))
    }
  } else if (model_option == 'beta') {
    message("Begin per-bin beta-model runs...")
    results <- foreach(i = bin_indices, .export = "Jitter", .packages = "betareg") %dopar% {
      predictor <- X[i, ]
      na_mask <- !is.na(predictor)
      predictor <- predictor[na_mask]
      if (length(predictor) < dim(Y)[2]) {
        return(list(model_predvals = list(), predictor_index = i))
      } else {
        Y <- Y[na_mask,]
      }
      models <- list()
      for (j in seq_len(ncol(Y))) {
        response <- Y[, j]
        response <- pmax(
          pmin(response, 1 - epsilon + Jitter(response)),
          epsilon - Jitter(response)
        )
        model <- betareg::betareg(response ~ predictor)
        model_sum <- summary(model)
        models[[colnames(Y)[j]]] <- c(model_sum$coefficients$mean[,1][2],
                                      model_sum$coefficients$mean[,4][2])
      }
      return(list(model_predvals = models,
                  predictor_index = i))
    }
  } else {
    stop("Invalid value for model_option parameter.")
  }

  time.taken <- Sys.time() - start.time
  message("Done")
  message("Time elapsed: ", round(time.taken, digits = 2),
          " ", attr(time.taken, "units"))

  return(results)
}


# Define jitter as a deterministic small perturbation based on value
Jitter <- function(x, factor = 5e-5) {
  runif(length(x), min = -factor, max = factor)
}


#' Fit a model using CN-signature exposures to predict each bin in a CN-dataset
#'
#' @description
#' Fit a model for each copy-number bin in Y given the signature exposures matrix X. \cr
#' The rows of X and columns of Y should match. \cr
#'
#' Because signature exposures for a sample sum to 1, fitting a model as described encounters the issue of multi-collinearity.
#' FitGIModelsStoCN overcomes that issue by first performing lasso regression and then fitting a linear model (LM).
#' Lasso regression makes use of the \pkg{glmnet} package while the \pkg{stats} is called to fit each LM.
#'
#' This function is a complement to \link{FitGIModelsCNtoS}.
#'
#' @param X *matrix.* A matrix of samples (x) by signatures (y). Values are signature exposures. \cr
#' `dim(X)[1]` must equal `dim(Y)[2]`
#' @param Y *S4/matrix.* A QDNAseq object that contains \code{copynumber}, and/or \code{segmented} assay slots. \cr
#' OR alternatively a matrix of bins (x) by samples (y) may be provided. The values being the copy-numbers.
#' @param bin_indices (optional) *numeric/integer vector.* A single integer or a vector of integers. \cr
#' Subsets the copy-number object to only run on the bins in these indices.
#' @param data_slot (optional) *character.* If passing a QDNAseq object, the CN data slot to use. \cr
#' Common options are "copynumber" (corrected or uncorrected read counts), "segmented", or "calls".
#' @param run_parallel (optional) *logical.* If TRUE run job in parallel. Will automatically use 1 less than total cores available.
#' @param cores (optional) *integer.* Alternatively/in-addition, the user may choose the number of cores to use. \cr
#' Must be <= 1 less than the number available.
#'
#' @return A list of lists. One entry for each bin. Each bin then has a list of: \cr
#' `model_predvals` \cr
#' `predictor_index`
#'
#' @export
FitGIModelsStoCN <- function (X, Y,
                              bin_indices = NULL,
                              data_slot = NULL,
                              run_parallel = FALSE,
                              cores = 1) {

  # Validate inputs
  stopifnot(inherits(X, 'matrix'))
  stopifnot(dim(X)[1] == dim(Y)[2])
  if (inherits(Y, c("QDNAseqCopyNumbers", "QDNAseqReadCounts"))) {
    if (!is.null(data_slot)) {
      stopifnot(data_slot %in% Biobase::assayDataElementNames(Y))
    }
    Y <- Biobase::assayDataElement(Y, data_slot)
  } else {
    stopifnot(inherits(Y, 'matrix'))
  }
  if (cores != 1) {
    require(foreach)
    available_cores <- parallel::detectCores() - 1
    stopifnot(available_cores >= cores && cores > 0)
    doMC::registerDoMC(cores)
  } else if (run_parallel) {
    available_cores <- parallel::detectCores() - 1
    doMC::registerDoMC(available_cores)
  }
  if (!is.null(bin_indices)) {
    stopifnot(all(bin_indices <= dim(Y)[1]))
  } else {
    bin_indices <- 1:nrow(Y)
  }

  # Run things in parallel if possible
  n_sigs <- ncol(X)
  message("Begin per-bin Lasso + LM runs...")
  start.time <- Sys.time()
  results <- foreach(i = bin_indices,
                     .packages = c("stats", "glmnet")) %dopar% {

                       # Extract bin-level values and prep for modelling
                       na_mask <- !is.na(Y[i, ])
                       if (sum(na_mask) < (n_sigs - 1)) {
                         return(list(model_predvals = list(), bin_index = i))
                       } else {
                         cn_outcome <- Y[i, na_mask]
                         X <- X[na_mask,]
                       }
                       coefs <- rep(NA, n_sigs)
                       pvalues <- rep(NA, n_sigs)

                       # Run lasso + lm
                       cv_fit <- glmnet::cv.glmnet(X, cn_outcome, alpha = 1, standardize = TRUE)
                       best_lambda <- cv_fit$lambda.min
                       lasso_model <- glmnet::glmnet(X, cn_outcome, alpha = 1, lambda = best_lambda)
                       lasso_coefs <- as.vector(coef(lasso_model))[-1]
                       predictor_to_drop <- which(abs(lasso_coefs) == min(abs(lasso_coefs)))
                       if (length(predictor_to_drop) == n_sigs) {
                         return(list(model_predvals = list(), bin_index = i))
                       } else {
                         x <- X[, -predictor_to_drop]
                       }
                       final_model <- lm(cn_outcome ~ x)

                       # Extract coefs and pvalues
                       good_pred_idx <- setdiff(1:dim(X)[2], predictor_to_drop)
                       good_pred_coef <- coef(final_model)[2:(1 + length(good_pred_idx))]
                       coefs[good_pred_idx] <- good_pred_coef
                       good_pred_pvals <- summary(final_model)$coefficients[2:(1 + length(good_pred_idx)), 4]
                       pvalues[good_pred_idx] <- good_pred_pvals

                       results <- list()
                       for (j in 1:n_sigs) {
                         results[[colnames(X)[j]]] <- c(coefs[j], pvalues[j])
                       }

                       return(list(model_predvals = results,
                                   bin_index = i))
                     }

  time.taken <- Sys.time() - start.time
  message("Done")
  message("Time elapsed: ", round(time.taken, digits = 2),
          " ", attr(time.taken, "units"))

  return(results)
}


#' Reformat genomic importance modelling output
#'
#' @description
#' Reformat genomic importance modelling output into a list of vectors.
#'
#' @param results *list of lists.* The output from FitGIModelsCNtoS.
#'
#' @return A list of vectors. Two vectors (pvals and coefficients) for each signature. \cr
#' All vectors have the same length as the provided number of bins to FitGIModelsCNtoS.
#'
#' @export
MakeGIVectors <- function (results) {

  sigs <- unique(unlist(lapply(results, function(x) names(x$model_predvals))))

  # Initialize empty list to store vectors for each matrix
  matrix_list <- list()
  for (sig in sigs) {
    matrix_list[[paste0(sig, "_bincoef")]] <- numeric(length(results))
    matrix_list[[paste0(sig, "_binpvals")]] <- numeric(length(results))
  }

  # Loop through 'results' to populate the vectors
  for (i in seq_along(results)) {
    model_predvals <- results[[i]]$model_predvals

    # Fill with 0 if model_predvals is empty
    if (length(model_predvals) == 0) {
      for (sig in sigs) {
        matrix_list[[paste0(sig, "_bincoef")]][i] <- NA
        matrix_list[[paste0(sig, "_binpvals")]][i] <- NA
      }

      # Otherwise extract values from list
    } else {
      for (sig in sigs) {
        if (!is.null(model_predvals[[sig]])) {
          matrix_list[[paste0(sig, "_bincoef")]][i] <- model_predvals[[sig]][1]
          matrix_list[[paste0(sig, "_binpvals")]][i] <- model_predvals[[sig]][2]
        } else {
          # Fill missing VS with 0
          matrix_list[[paste0(sig, "_bincoef")]][i] <- NA
          matrix_list[[paste0(sig, "_binpvals")]][i] <- NA
        }
      }
    }
  }
  return(matrix_list)
}


