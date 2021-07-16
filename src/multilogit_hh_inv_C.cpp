// This code was written by Kyle McEvoy working off the pseudocode included in a 2006 paper
// by Chris C. Holmes and Leonhard Held and with corrections suggested by
// by Ralf van der Lans as found in the 2011 response:

// 1. Holmes,  C.  C.  and  Held,  L.  (2006). Bayesian  auxiliary  variable  models
// for  binary  and  multinomialregression. Bayesian Analysis, 1(1):145–168.

// 2. Holmes, C. C. and Held, L. (2011). Response to van der Lans.
// Bayesian Analysis, 6(2):357-358.

#include "RcppArmadillo.h"
#include "lambda_sampler.h"
#include "trunc_logis.h"
using namespace Rcpp;
using namespace arma;


//' Multinomial Logistic Regression using the Holmes-Held Method
//' 
//' @description This function implements the Holmes-Held method for
//' multinomial logistic regression described in their 2006 paper
//' 'Bayesian auxiliary variable models for binary and multinomial regression'
//' in Bayesian Analysis. The C++ code was written using the pseudo-code in
//' this paper as a template.
//' 
//'
//' @param Y An N by C numeric matrix where the ith row is a set of
//' indicators for observation i of N total observations giving which
//' of the C categories the observation is classified into.
//' @param X An N by P numeric matrix where the ith row gives
//' the values of the predictor variables for the ith outcome observation. 
//' @param v a P by P numeric matrix giving the covariance matrix of coefficients.
//' The function only accepts one matrix for all categories.
//' @param n_sample positive integer giving the number of samples to draw as
//'  output after burn-in.
//' @param n_burn non-negative integer giving the number of samples of burn-in
//'  before the chain output is saved.
//' @return List object containing posterior_coef, the chain of coefficient
//' values as an P by C by n_sample array. And posterior_prob a N by C by
//' n_sample array containing the calculated probabilities of the observations
//' being classified into each of the C categories.
//' @examples 
//' Y <- matrix(0, nrow = 150, ncol = 3)
//' Y[1:50, 1] <- 1
//' Y[51:100, 2] <- 1
//' Y[101:150, 3] <- 1
//' X <- cbind(1, iris[ , 1:4])
//' X <- as.matrix(X)
//' v <- diag(10, ncol(X))
//' out <- multilogit_holmesheld_C(Y, X, v, n_sample = 4000, n_burn = 2000)
//' 
// [[Rcpp::export]]
List multilogit_hh_inv_C(
    NumericMatrix Y_, 
    NumericMatrix X_, 
    NumericMatrix v_, 
    size_t n_sample = 1000,
    size_t n_burn = 200,
    bool probs = true,
    bool progress = true
){
  
  /**
   * Sizes
   */
  size_t N = Y_.nrow();
  size_t Q = Y_.ncol();
  size_t P = X_.ncol();
  
  /**
   * Move rcpp types to armadillo for linear algebra. ???Reconsider this for big data with pointers, pass through, etc.
   */
  //updated code to use the Rcpp as function which converts from R to C++ data types.
  
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat v = Rcpp::as<arma::mat>(v_);
  
  
  /**
   * Create storage objects
   */
  arma::mat Z(N, Q - 1, fill::zeros);
  arma::cube lambda(N,N,Q-1,fill::zeros);
  
  //output objects
  
  arma::cube beta(P, Q, n_sample + n_burn, fill::zeros);
  arma::cube prob(N, Q, n_sample + n_burn, fill::zeros); 
  
  arma::cube beta_out(P, Q, n_sample, fill::zeros);
  arma::cube prob_out(N, Q, n_sample, fill::zeros);
  
  
  // R code asks for this vectorized identity matrix approach:
  // Lambda <- array(data = as.vector(diag(N)), dim = c(N, N, (Q - 1)))
  // will just 0's work? or will the following work? 
  arma::mat identity(N,N,fill::eye);
  for(size_t i=0;i<Q-1;i++){
    lambda.slice(i) = identity; 
  }
  
  //creates a NxQ NumericMatrix filled with 0s. This will be interpreted as logical in the right argument of trunc_logis.
  
  NumericMatrix right_switch(N,Q);
  
  for(size_t i=0;i<N;i++){
    for(size_t j=0;j<(Q-1);j++){
      if(Y(i,j)!=0 ){
        right_switch(i,j) = 1; 
      }
    }
  }
  
  
  arma::mat Y_sub = Y.submat(0,0, N-1, Q-2);
  double y_sum = arma::accu(Y_sub);
  
  // arma uses unsigned int uvecs generated by arma::find as element pointers for objects.
  
  arma::uvec ids = arma::find(Y_sub == 1);
  arma::uvec not_ids = arma::find(Y_sub == 0);
  
  //trunc_logis takes arguments (location, scale, right) when right = TRUE truncation: (0, inf), False: (-inf, 0)
  
  for(size_t i = 0; i < y_sum; i++)
  {
    Z(ids(i)) = trunc_logis(0, 1, TRUE);
    
  }
  
  for(size_t i = 0; i < (N * (Q - 1)) - y_sum; i++)
  {
    
    Z(not_ids(i)) = trunc_logis(0, 1, FALSE);
    
  }
  
  
  /**
   * Cache frequently used values
   */
  arma::mat v_inv = inv(v);
  
  
  /**
   * Temporary objects needed inside MCMC
   */
  arma::mat V(P,P, fill::zeros);
  arma::mat L(P,P, fill::zeros);
  arma::vec B(P, fill::zeros);
  arma::mat lambda_inv(N,N, fill::zeros);
  arma::vec Tp(P, fill::zeros);
  double m{ 0 };
  double R{ 0 };
  double C{ 0 };
  arma::mat beta_slice(P, Q, fill::zeros);
  
  /** 
   * MCMC
   */
  for(size_t i=0;i<n_sample + n_burn;i++){
    
    // status report
    if( progress == true && ((i%1000) == 0 || (i + 1) == n_burn + n_sample)){
      Rcout << "iteration " << i + 1 << " of " << n_sample + n_burn << "\n";
    }
    
    
    for(size_t q=0;q<(Q-1);q++){
      
      lambda_inv = inv(diagmat(lambda.slice(q)));
      
      V = inv_sympd( X.t() * lambda_inv * X + v_inv); 
      L = chol(V,"lower"); // we want lower triangular right? 
      
      B = V * X.t() * lambda_inv * Z.col(q);
      
      Tp = Rcpp::rnorm(P, 0, 1);
      
      beta.subcube(0, q, i, (P - 1), q, i) = B + (L * Tp);
      
      arma::mat beta_slice = beta.slice(i);
      
      
      for (size_t j = 0; j < N-1;j++){
        
        m = as_scalar(X.row(j) * beta_slice.col(q));
        
        C = sum(arma::exp(X.row(j) * beta_slice)) - arma::as_scalar(exp(X.row(j) * beta_slice.col(q)));
        
        double loc = m - std::log(C);
        
        
        if (right_switch(j,q))
        {
          
          Z(j,q) = trunc_logis(loc, 1, TRUE) + log(C);
          
        }
        
        else
        {
          Z(j,q) = trunc_logis(loc, 1, FALSE) + log(C);
        }
        
        R = Z(j,q) - m;
        
        // #abs(R) is necessary as the parameter r must be a positive real number.
        lambda(j,j,q) = lambda_sampler(R);
        
      } // end j loop through observations
      
    } // end q loop through categories
    
    beta_slice = beta.slice(i);
    
    if (probs == true){
      
      prob.slice(i) = exp(X * beta_slice);
      
      arma::vec prob_slice_sums = sum(prob.slice(i), 1);
      
      for (size_t j = 0; j <= N - 1; j++)
      {
        for (size_t q = 0; q <= Q - 1; q++)
        {
          prob(j, q, i) = prob(j, q, i) / as_scalar(prob_slice_sums(j)); 
          
        }
        
      }
      
    }
    //}
    
    
  } //end i loop over mcmc iterations 
  
  beta_out = beta.subcube(0, 0, n_burn, P - 1, Q - 1, n_burn + n_sample - 1);
  
  if (probs == true){
    
    prob_out = prob.subcube(0, 0, n_burn, N - 1, Q - 1, n_burn + n_sample - 1);
    
  }
  
  
  if (probs == true){
    
    return(List::create(
        _["posterior_prob"] = prob_out,
        _["posterior_coef"] = beta_out
    ));
    
  }
  
  else{
    return(List::create(
        _["posterior_coef"] = beta_out
    ));
    
  }
  
  
  
}


