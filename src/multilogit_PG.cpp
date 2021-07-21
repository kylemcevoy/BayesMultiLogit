// This file was written by Kyle R. McEvoy and Jared Fisher adapting the file MultLogitPG.R which was
// part of the R package BayesLogit as found at https://github.com/jwindle/BayesLogit-Thesis,
// distributed under the GNU General Public License version 3 or later and without ANY warranty,
// implied or otherwise.

// Original code copyright 2013 Nick Polson, James Scott, and Jesse Windle. 
// Rewritten for C++ by Kyle McEvoy and Jared Fisher 2021.

//

#include "RcppArmadillo.h"
#include "helloPG.h"


using namespace Rcpp;
using namespace arma;


//' Multinomial Logistic Regression using the Polya-Gamma Method
//' 
//' @description WARNING: This function can result in R freezing in a non-interruptable state.
//' 
//' This function implements the Polya-Gamma method for
//' multinomial logistic regression. Rewritten for C++ by Jared Fisher
//' and Kyle McEvoy, but code originally written by Jesse Windle.
//' 
//' Copyright 2013 Nick Polson, James Scott, and Jesse Windle.
//' This file is part of BayesLogit, distributed under the GNU General Public
//' License version 3 or later and without ANY warranty, implied or otherwise.
//' 
//' @param Y numeric Matrix An N by C numeric matrix where the ith row is a set of
//' indicators for observation i of N giving which of the C categories the
//' observation is classified into.
//' @param X numeric Matrix An N by P numeric matrix where the ith row gives
//' the values of the predictor variables for the ith outcome observation. The first column of 
//' X should be an intercept column of 1's.
//' Non-intercept X columns should be centered and scaled by their standard deviations for best results. 
//' @param n_sample positive integer giving the number of samples to draw as
//'  output after burn-in.
//' @param n_burn non-negative integer giving the number of samples of burn-in
//'  before the chain output is saved.
//' @param progress logical If TRUE, the function will report its progress at every
//' thousandth iteration.
//' @param probs logical If TRUE, the function will calculate and return the probabilities
//' for each chain sample of each observation falling into the categories
//' based on the fitted coefficients.
//' @return List object containing posterior_coef, the chain of coefficient
//' values as an P by C by n_sample array. And, if probs is set to TRUE, 
//' posterior_prob a N by C by n_sample array containing the calculated 
//' probabilities of the observations being classified into each of the C categories.
//' 
//' @examples Y <- matrix(0, nrow = 150, ncol = 3)
//' Y <- sapply(c(1,2,3), function(x) Y[, x] <- as.numeric((as.numeric(iris$Species) == x) )) 
//' X <- scale(iris[ , 1:4])
//' X <- cbind(1, X)
//' out2 <- multilogit_PG_C(Y, X, n_sample = 2000, n_burn = 500, probs = TRUE, progress = TRUE)
//' 
// [[Rcpp::export]]
List multilogit_PG_C(arma::mat const &Y,
                     arma::mat const &X, 
                     size_t n_sample = 1000,
                     size_t n_burn = 200,
                     bool probs = true,
                     bool progress = true) {
  

  size_t N = Y.n_rows;
  size_t Q = Y.n_cols;
  size_t P = X.n_cols;
  
  arma::mat m_0(P, Q-1, fill::zeros);
  arma::cube P_0(P, P, Q - 1, fill::zeros);
  
  arma::mat w(N, Q, fill::zeros);
  arma::mat beta(P, Q, fill::zeros);
  
  arma::cube beta_out(P, Q, n_sample, fill::zeros);
  arma::cube prob_out(N, Q, n_sample, fill::zeros);
  
  arma::vec n(N, fill::ones);
  
  arma::mat y_sub = Y.submat(0, 0, N-1, Q-2);
  
  arma::mat kappa(N, Q - 1, fill::zeros); 
   
  for (size_t i = 0; i < Q - 1; i++)
    {
      kappa.col(i) = (y_sub.col(i) - 0.5) % n;
    
    }
  
  arma::mat b_0(P, Q - 1, fill::zeros);
  
  for (size_t j = 0; j < (Q - 1); j++)
  {
    
    b_0.col(j) = P_0.slice(j) * m_0.col(j);
    
  }
  
  for (size_t i = 0; i < (n_burn + n_sample); i++)
  {
    if(progress == true && ((i%1000) == 0 || (i + 1) == n_burn + n_sample)){
      Rcout << "iteration " << i + 1 << " of " << n_sample + n_burn << "\n";
    }
    
    for (size_t j = 0; j < Q - 1; j++)
    {
      
      arma::mat beta_woj(P, Q - 1, fill::zeros);
      
      // There's gotta be a better way to do this, but .shed_col() was throwing an error.
      
      arma::vec zeros(P, fill::zeros);
      
      beta_woj = beta;
      
      beta_woj.col(j) = zeros;
      
      arma::mat exp_probs = exp(X * beta_woj);  
      
      arma::mat A = sum(exp_probs, 1);
      
      arma::vec c_j = log(A);
      
      arma::mat eta_j = (X * beta.col(j)) - c_j;
    
      w.col(j) = rpg(n, eta_j);
      
      // Need element-wise product of X and w. C++ doesn't vectorize product
      
      arma::mat X_omega(N, P, fill::zeros);
      
      for (size_t p = 0; p < P; p++)
      {
        X_omega.col(p) = X.col(p) % w.col(j);
        
      }
      
      arma::mat PL_j = X.t() * (X_omega);
      arma::mat bl_j = X.t() * (kappa.col(j) + (c_j % w.col(j)));
      
      arma::mat P1_j = PL_j + P_0.slice(j);
      
      arma::mat v1_j = inv(P1_j);
      
      arma::mat m1_j = v1_j * (bl_j + b_0.col(j));
      
      arma::mat sqrtv1_j = chol(v1_j, "lower");
      
      arma::vec Tp = Rcpp::as<arma::vec>(Rcpp::rnorm(P));
      
      beta.col(j) = m1_j + sqrtv1_j * Tp;
      
      if (i >= n_burn)
      {
        beta_out.subcube(0, j, (i - n_burn), (P - 1), j, (i - n_burn)) = beta.col(j);
        
        
        
        
        
      }
      
      
      
    }

    if (probs == true && (i >= n_burn))
    {
      for(size_t row = 0; row < N ; row++)
      {
        prob_out.subcube(row, 0, i - n_burn, row, Q - 1, i - n_burn) =
          exp(X.row(row) * beta_out.slice(i - n_burn)) / sum(  exp(X.row(row) * beta_out.slice(i - n_burn))  );
      }
    }

  }
  


if (probs == true) {
  
  return(List::create(
      _["posterior_prob"] = prob_out,
      _["posterior_coef"] = beta_out
  ));
}
else {
  
  return(List::create(
    _["posterior_coef"] = beta_out
  ));
}
  
}

