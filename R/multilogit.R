#' Multinomial Logistic Regression using Data Augmentation
#' 
#' @description This function implements a data augmentation method for
#' multinomial logistic regression using MCMC. The data augmentation method is
#' outlined in a paper by Jared D. Fisher and Kyle R. McEvoy titled
#' "Bayesian Multinomial Logistic Regression for Numerous Categories",
#' currently a work in progress. The sampler uses a Random-Walk Metropolis
#' algorithm for choosing proposals and acceptances.
#' 
#'
#' @param Y An N by C numeric matrix where the ith row is a set of
#' indicators for observation i of N total observations giving which
#' of the C categories the observation is classified into.
#' @param X An N by P numeric matrix where N is the total number of
#' observations and P is the total number of predictor variables the ith row
#' gives the values of the predictor variables for the ith outcome observation. 
#' @param n_sample positive integer giving the number of samples to draw as
#'  output after burn-in.
#' @param n_burn non-negative integer giving the number of samples of burn-in
#'  before the chain output is saved.
#' @param n_sigma_check positive integer giving the period for the number of
#' samples on which to automatically tune the step size of the random walk.
#' @param step_size finite positive standard deviation for the random-walk
#' proposal distribution.
#' @param prior character with either values "flat" or "normal"
#' @param prior_mean a vector of length equal to P the number of predictors.
#' Giving the mean for the normal prior on the coefficient vector. Only use if
#' prior = "normal".
#' @param prior_var a P by P matrix giving the covariance matrix of the prior
#' coefficients. Note that covariance structure is assumed to be constant across
#' categories. Must be a positive-definite matrix. Only use with prior =
#' "normal".
#' @param reference_cat Either NULL or an integer between 1 and C where C is
#' the total number of categories. If left NULL, coefficients for all
#' categories will be generated, leading to a loss of identifiability in the
#' coefficients. If an integer, all of the coefficients for that category will
#' be held constant at 0.
#' @param probs logical; if TRUE, return fitted category probabilities for each
#' retained draw.
#' @param progress logical; if TRUE, print periodic progress updates.
#' @return List object containing posterior_coef, the chain of coefficient
#' values as an P by C by n_sample array. And posterior_prob a N by C by
#' n_sample array containing the calculated probabilities of the observations
#' being classified into each of the C categories.
#' @examples
#' Y <- matrix(0, nrow = 150, ncol = 3)
#' Y[1:50, 1] <- 1
#' Y[51:100, 2] <- 1
#' Y[101:150, 3] <- 1
#' X <- cbind(1, iris[ , 1:4])
#' X <- as.matrix(X)
#' out <- multilogit(Y, X, n_sample = 3000, n_burn = 1500, n_sigma_check = 20,
#'  step_size = 0.1, prior = "normal", reference_cat = 1)
#'  
#' out_2 <- multilogit(Y, X, n_sample = 2000, n_burn = 500, prior = "flat")
#' 

multilogit <- function(Y, X, n_sample = 1000, n_burn = 200, n_sigma_check = 20, step_size = 0.1, prior = "flat",
                      prior_mean = NULL, prior_var = NULL, reference_cat = NULL, probs = FALSE, progress = TRUE){
  data <- .validate_data(Y, X)
  Y <- data$Y
  X <- data$X
  n_cat <- ncol(Y)
  n_pred <- ncol(X)

  n_sample <- .validate_scalar_integer(n_sample, "n_sample", 1L)
  n_burn <- .validate_scalar_integer(n_burn, "n_burn", 0L)
  .validate_iteration_total(n_sample, n_burn)
  n_sigma_check <- .validate_scalar_integer(n_sigma_check, "n_sigma_check", 1L)
  probs <- .validate_flag(probs, "probs")
  progress <- .validate_flag(progress, "progress")

  if (length(prior) != 1L || !is.character(prior) || is.na(prior) ||
      !prior %in% c("flat", "normal")) {
    stop("prior must be either 'flat' or 'normal'.", call. = FALSE)
  }
  if (length(step_size) != 1L || !is.numeric(step_size) || is.na(step_size) ||
      !is.finite(step_size) || step_size <= 0) {
    stop("step_size must be a finite positive number.", call. = FALSE)
  }
  if (prior == "flat" && (!is.null(prior_mean) || !is.null(prior_var))) {
    stop("prior_mean and prior_var may only be supplied when prior = 'normal'.", call. = FALSE)
  }
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
  reference_cat <- .validate_reference(reference_cat, n_cat)
  
  
  out <- multilogit_C(Y, X, n_sample = n_sample, n_burn = n_burn, n_sigma_check = n_sigma_check,
                  prior = prior, step_size = step_size, prior_mean = prior_mean,
                  prior_var = prior_var, reference_cat = reference_cat, probs = probs, progress = progress)
  
  out
}
