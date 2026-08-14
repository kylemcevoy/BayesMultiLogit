#' Multinomial Logistic Regression using the Polya-Gamma Method
#' 
#' 
#' @description This is an R wrapper to the function \code{multilogit_PG_C}.
#' This function implements the Polya-Gamma method for
#' multinomial logistic regression. Rewritten for C++ by Jared Fisher
#' and Kyle McEvoy, but code originally written by Jesse Windle, James Scott and Nick Polson.
#' 
#' Copyright 2013 Nick Polson, James Scott, and Jesse Windle.
#' This file is part of BayesLogit, distributed under the GNU General Public
#' License version 3 or later and without ANY warranty, implied or otherwise.
#' @param Y An N by C numeric matrix where the ith row is a set of
#' indicators for observation i of N giving which of the C categories the
#' observation is classified into.
#' @param X An N by P numeric matrix where the ith row gives
#' the values of the predictor variables for the ith outcome observation. 
#' The first column of X should be an intercept column of 1's.
#' Non-intercept X columns should be centered and scaled by their standard deviations for best results.  
#' @param n_sample positive integer giving the number of samples to draw as
#'  output after burn-in.
#' @param n_burn non-negative integer giving the number of samples of burn-in
#'  before the chain output is saved.
#' @param probs logical If TRUE probabilities are calculated and returned.
#' @param progress logical If TRUE, the function reports progress every thousandth iteration.
#' @return List object containing posterior_coef, the chain of coefficient
#' values as an P by C by n_sample array. If probs are TRUE, posterior_prob a N by C by
#' n_sample array containing the calculated probabilities of the observations
#' being classified into each of the C categories is also returned in the list.
#' 
#' @examples Y <- matrix(0, nrow = 150, ncol = 3)
#' Y <- sapply(c(1,2,3), function(x) Y[, x] <- as.numeric((as.numeric(iris$Species) == x) )) 
#' X <- scale(iris[ , 1:4])
#' X <- cbind(1, X)
#' out2 <- multilogit_PG(Y, X, n_sample = 2000, n_burn = 1000, probs = TRUE, progress = TRUE)
multilogit_PG <- function(Y, X, n_sample = 1000L, n_burn = 200L, probs = FALSE, progress = TRUE){
  data <- .validate_data(Y, X)
  Y <- data$Y
  X <- data$X
  n_sample <- .validate_scalar_integer(n_sample, "n_sample", 1L)
  n_burn <- .validate_scalar_integer(n_burn, "n_burn", 0L)
  .validate_iteration_total(n_sample, n_burn)
  probs <- .validate_flag(probs, "probs")
  progress <- .validate_flag(progress, "progress")
  .warn_missing_intercept(X)
  
  out <- multilogit_PG_C(Y = Y, X = X, n_sample = n_sample, n_burn = n_burn, probs = probs, progress = progress)
  
  out
}
