#' Multinomial Logistic Regression using the Holmes-Held Method
#' 
#' @description This is an R wrapper to the C++ functions multilogit_holmesheld_C and multilogit_hh_inv_C.
#' This function implements the Holmes and Held method for
#' multinomial logistic regression described in their 2006 paper
#' 'Bayesian auxiliary variable models for binary and multinomial regression'
#' in Bayesian Analysis. The C++ code was written by Kyle R. McEvoy using the pseudo-code in
#' this paper as a template.
#' 
#'
#' @param Y An N by C numeric matrix where the ith row is a set of
#' indicators for observation i giving which
#' of the C categories the observation is classified into.
#' @param X An N by P numeric matrix where the ith row gives
#' the values of the predictor variables for the ith observation. The first column of X
#' should be an intercept column of 1s. Non-intercept X columns
#' should be centered and scaled by their standard deviations for best results. 
#' @param v a P by P numeric matrix giving the covariance matrix of coefficients.
#' The function only accepts one matrix for all categories.
#' @param n_sample positive integer giving the number of samples to draw as
#'  output after burn-in.
#' @param n_burn non-negative integer giving the number of samples of burn-in
#'  before the chain output is saved.
#' @param probs logical If TRUE categorical probabilities will be calculated and returned for
#' non-burn samples.
#' @param progress logical If TRUE the function will report its progress each thousandth iteration.
#' @param fast logical If TRUE the wrapper function will call \code{multilogit_hh_inv_C} as opposed to 
#' \code{mulitlogit_holmesheld_C}. \code{multilogit_hh_inv_C} uses faster matrix inverse functions so is faster, 
#' but is more prone to errors due to numerical tolerance issues when inverting the matrices.
#' @family Holmes-Held methods.
#' @seealso \code{multilogit_holmesheld_C} and \code{multilogit_hh_inv_C} for the C++ functions without
#' error checking.
#' @return List object containing posterior_coef, the chain of coefficient 
#' values as an P by C by n_sample array. If probs is set to TRUE, the 
#' posterior_prob a N by C by n_sample array containing the calculated 
#' probabilities of the N observations being classified into each of the C categories 
#' will also be returned.
#' @examples 
#' Y <- matrix(0, nrow = 150, ncol = 3)
#' Y <- sapply(c(1,2,3), function(x) Y[, x] <- as.numeric((as.numeric(iris$Species) == x) )) 
#' X <- cbind(1, iris[ , 1:4])
#' X <- as.matrix(X)
#' v <- diag(1, ncol(X))
#' out <- multilogit_holmesheld(Y, X, v, n_sample = 2000, n_burn = 1000)
#' 
#' 
multilogit_holmesheld <- function(Y,
                                  X,
                                  v = diag(1, nrow = ncol(X)),
                                  n_sample = 1000L,
                                  n_burn = 200L,
                                  probs = TRUE,
                                  progress = TRUE,
                                  fast = FALSE) {
  # Here, we consider Y to be a matrix of indicators with number of rows equal to the number of subjects, and the number of columns equal to the number of categories. 
  # In other words, Y's dimensions are n_sub x n_cat
  Y = as.matrix(Y)
  X = as.matrix(X)
  v = as.matrix(v)
  
  n_sub = nrow(Y)
  n_cat = ncol(Y)
  # X is a n_sub x n_pred matrix of covariates, aka the design matrix. 
  n_pred = ncol(X)
  
  if(nrow(X) != n_sub){
    stop("unequal number of rows in X and Y")
  }
  
  if(!isSymmetric(v) | any(eigen(v)$values < 0)   ){
    stop("Variance matrix v is not a symmetric positive semi-definite matrix.")
  }
  
  if(!identical(X[,1], rep(1, nrow(X)))){
    warning("Function expects first column of the design matrix to be an intercept column.")
  }
  
  if(!all.equal(dim(v), c(n_pred, n_pred))){
    stop("v should be a square matrix with each dimension equal to the number of predictors (including intercept).")
  }
  
  
  if (fast == TRUE){
    out <- multilogit_hh_inv_C(Y = Y, X = X, v = v, n_sample = n_sample, n_burn = n_burn,
                               probs = probs, progress = progress)
    
  }
  
  else {
    out <- multilogit_holmesheld_C(Y = Y, X = X, v = v, n_sample = n_sample, n_burn = n_burn,
                                   probs = probs, progress = progress)
  }
 
  
  return(out)
}