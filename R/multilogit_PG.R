#' Multinomial Logistic Regression using the Polya-Gamma Method
#' 
#' @description This is an R wrapper to the function multilogit_PG_C.
#' This function implements the Polya-Gamma method for
#' multinomial logistic regression. Rewritten for C++ by Jared Fisher
#' and Kyle McEvoy, but code originally written by Jesse Windle.
#' 
#' Copyright 2013 Nick Polson, James Scott, and Jesse Windle.
#' This file is part of BayesLogit, distributed under the GNU General Public
#' License version 3 or later and without ANY warranty, implied or otherwise.
#' @param Y An N by C numeric matrix where the ith row is a set of
#' indicators for observation i of N giving which of the C categories the
#' observation is classified into.
#' @param X An N by P numeric matrix where the ith row gives
#' the values of the predictor variables for the ith outcome observation. 
#' @param n_sample positive integer giving the number of samples to draw as
#'  output after burn-in.
#' @param n_burn non-negative integer giving the number of samples of burn-in
#'  before the chain output is saved.
#' @return List object containing posterior_coef, the chain of coefficient
#' values as an P by C by n_sample array. And posterior_prob a N by C by
#' n_sample array containing the calculated probabilities of the observations
#' being classified into each of the C categories.
#' 
#' @examples Y <- matrix(0, nrow = 150, ncol = 3)
#' Y[1:50, 1] <- 1
#' Y[51:100, 2] <- 1
#' Y[101:150, 3] <- 1
#' X <- cbind(1, iris[ , 1:4])
#' X <- as.matrix(X)
#' out2 <- multilogit_PG(Y, X, n_sample = 4000, n_burn = 2000)
multilogit_PG <- function(Y, X, n_sample = 1000L, n_burn = 200L){
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
  
  out <- multilogit_PG_C(Y = Y, X = X, n_sample = n_sample, n_burn = n_burn)
  
  return(out)
}