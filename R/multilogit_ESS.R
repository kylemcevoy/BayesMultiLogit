#' Multinomial Logistic Regression using Data Augmentation (Elliptical Slice Sampler)
#' 
#' @description This function is an R wrapper for the function \code{multilogit_C_ESS}. The wrapper
#' contains error-checking that the C++ function lacks. \code{multilogit_C_ESS} implements a
#' data augmentation method for multinomial logistic regression using MCMC.
#' The data augmentation method is outlined in a paper by Jared D. Fisher and Kyle R. McEvoy titled
#' "Bayesian Multinomial Logistic Regression for Numerous Categories", currently a work in progress.
#' 
#' The sampler uses the elliptical slice sampling algorithm
#' from the Murray, Adams, Mackay (2010) paper Elliptical Slice Sampling in the Journal of
#' Machine Learning Research. This sampler requires a multivariate normal prior on the betas.
#'
#' @param Y An N by C numeric matrix where the ith row is a set of
#' indicators for observation i of N total observations giving which
#' of the C categories the observation is classified into.
#' @param X An N by P numeric matrix where N is the total number of
#' observations and P is the total number of predictor variables the ith row
#' gives the values of the predictor variables for the ith outcome observation. 
#' The first column should be an intercept column of 1s. Non-intercept X columns
#' should be centered around 0 and scaled by their standard deviations for best results. 
#' @param n_sample positive integer giving the number of samples to draw as
#'  output after burn-in.
#' @param n_burn non-negative integer giving the number of samples of burn-in
#'  before the chain output is saved.
#' @param prior_mean a vector of length equal to P the number of predictors.
#' Giving the mean for the normal prior on the coefficient vector.
#' @param prior_var a P by P matrix giving the covariance matrix of the prior
#' coefficients. Note that covariance structure is assumed to be constant across
#' categories. Must be a positive-definite matrix.
#' @param reference_cat Either NULL or an integer between 1 and C where C is
#' the total number of categories. If left NULL, coefficients for all
#' categories will be generated, leading to a loss of identifiability in the
#' coefficients. If an integer, all of the coefficients for that category will
#' be held constant at 0.
#' @param probs logical If TRUE, for each non-burn chain sample of coefficients
#' the categorical probabilities of each observation will be calculated and returned.
#' @param progress logical If TRUE, the function will report progress at each thousandth
#' iteration.
#' @return List object containing posterior_coef, the chain of coefficient
#' values as an P by C by n_sample array. If probs are TRUE, it will also contain
#' posterior_prob a N by C by n_sample array containing the calculated probabilities of the observations
#' being classified into each of the C categories.
#' @examples 
#' Y <- matrix(0, nrow = 150, ncol = 3)
#' Y <- sapply(c(1,2,3), function(x) Y[, x] <- as.numeric((as.numeric(iris$Species) == x) )) 
#' X <- scale(iris[ , 1:4])
#' X <- cbind(1, X)
#' out <- multilogit_ESS(Y, X, n_sample = 3000, n_burn = 1500, reference_cat = 1)
#'  
#' out_2 <- multilogit_ESS(Y, X, n_sample = 2000, n_burn = 500,
#'                          prior_var = diag(x = 1, nrow = 5))
#'                          
#' out_3 <- multilogit_ESS(Y, X, n_sample = 5000, n_burn = 2000,
#'                        prior_var = diag(x = 1, nrow = 5),
#'                        reference_cat = 1)
#' 
multilogit_ESS <- function(Y, X, n_sample = 1000, n_burn = 200, prior_mean = NULL,
                           prior_var = NULL, reference_cat = NULL, probs = FALSE, progress = TRUE){
  data <- .validate_data(Y, X)
  Y <- data$Y
  X <- data$X
  n_pred <- ncol(X)

  n_sample <- .validate_scalar_integer(n_sample, "n_sample", 1L)
  n_burn <- .validate_scalar_integer(n_burn, "n_burn", 0L)
  .validate_iteration_total(n_sample, n_burn)
  probs <- .validate_flag(probs, "probs")
  progress <- .validate_flag(progress, "progress")
  .warn_missing_intercept(X)

  if (!is.null(prior_mean)) {
    if (!is.numeric(prior_mean) || length(prior_mean) != n_pred ||
        anyNA(prior_mean) || any(!is.finite(prior_mean))) {
      stop("prior_mean must be a finite numeric vector with one value per predictor.", call. = FALSE)
    }
    prior_mean <- as.numeric(prior_mean)
  }
  if (!is.null(prior_var)) {
    prior_var <- .validate_covariance(prior_var, n_pred, "prior_var")
  }
  reference_cat <- .validate_reference(reference_cat, ncol(Y))
  
  
  out <- multilogit_C_ESS(Y, X, n_sample = n_sample, n_burn = n_burn, 
                      prior_mean = prior_mean,
                      prior_var = prior_var, reference_cat = reference_cat, probs = probs, progress = progress)
  
  out
}
