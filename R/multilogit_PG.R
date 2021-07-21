#' Multinomial Logistic Regression using the Polya-Gamma Method
#' 
#' 
#' @description WARNING: This function can result in R freezing in a non-interruptable state.
#' 
#' This is an R wrapper to the function \code{multilogit_PG_C}.
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
multilogit_PG <- function(Y, X, n_sample = 1000L, n_burn = 200L, probs = TRUE, progress = TRUE){
  # Here, we consider Y to be a matrix of counts with number of rows equal to the number of subjects, and the number of columns equal to the number of categories. 
  # In other words, Y's dimensions are n_sub x n_cat
  Y = as.matrix(Y)
  X = as.matrix(X)
  
  n_sub = nrow(Y)
  n_cat = ncol(Y)
  # X is a n_sub x n_pred matrix of covariates, aka the design matrix. 
  n_pred = ncol(X)
  
  if(nrow(X) != n_sub){
    stop("unequal number of rows in X and Y")
  }
  
  
  if(!identical(X[,1], rep(1, nrow(X)))){
    warning("Function expects first column of the design matrix to be an intercept column.")
  }
  
  out <- multilogit_PG_C(Y = Y, X = X, n_sample = n_sample, n_burn = n_burn, probs = probs, progress = progress)
  
  return(out)
}