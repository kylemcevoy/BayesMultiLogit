#include "dmvnrm_arma.h"
#include "RcppArmadillo.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//' Multinomial Logistic Regression using Data Augmentation (Metropolis-Hastings)
//' 
//' This function implements a data augmentation method for
//' multinomial logisitic regression using MCMC. The data augmentation method is
//' outlined in a paper by Jared D. Fisher and Kyle R. McEvoy titled
//' "Bayesian Multinomial Logistic Regression for Numerous Categories",
//' currently a work in progress. The sampler uses a Random-Walk Metropolis
//' algorithm for choosing proposals and acceptances.
//' 
//'
//' @param Y An N by C numeric matrix where the ith row is a set of
//' indicators for observation i of N total observations giving which
//' of the C categories the observation is classified into.
//' @param X An N by P numeric matrix where N is the total number of
//' observations and P is the total number of predictor variables the ith row
//' gives the values of the predictor variables for the ith outcome observation.
//' The first column of X should be an intercept column of 1s.
//' Non-intercept X columns should be centered and scaled by their standard deviations for best results. 
//' @param n_sample positive integer giving the number of samples to draw as
//'  output after burn-in.
//' @param n_burn non-negative integer giving the number of samples of burn-in
//'  before the chain output is saved.
//' @param n_sigma_check non-negative integer gives the period for the number of
//' samples on which to automatically tune the step size of the random walk.
//' @param prior character with either values "flat" or "normal"
//' @param prior_mean a vector of length equal to P the number of predictors.
//' Giving the mean for the normal prior on the coefficient vector. Only use if
//' prior = "normal".
//' @param prior_var a P by P matrix giving the covariance matrix of the prior
//' coefficients. Note that covariance structure is assumed to be constant across
//' categories. Must be a positive semi-definite matrix. Only use with prior = 
//' "normal".
//' @param reference_cat Either NULL or an integer between 1 and C where C is
//' the total number of categories. If left NULL, coefficients for all
//' categories will be generated, leading to a loss of identifiability in the
//' coefficients. If an integer, all of the coefficents for that category will
//' be held constant at 0.
//' @param probs logical If TRUE, for each non-burn chain sample of coefficients
//' the categorical probabilities of each observation will be calculated and returned.
//' @param progress logical If TRUE, the function will report progress at each thousandth
//' iteration.
//' @return List object containing posterior_coef, the chain of coefficient
//' values as an P by C by n_sample array. If probs are TRUE, it will also contain
//' posterior_prob a N by C by n_sample array containing the calculated probabilities of the observations
//' being classified into each of the C categories.
//' @examples 
//' Y <- matrix(0, nrow = 150, ncol = 3)
//' Y <- sapply(c(1,2,3), function(x) Y[, x] <- as.numeric((as.numeric(iris$Species) == x) )) 
//' X <- scale(iris[ , 1:4])
//' X <- cbind(1, X)
//' out <- multilogit_C(Y, X, n_sample = 3000, n_burn = 1500, n_sigma_check = 20,
//'  step_size = 0.1, prior = "normal", reference_cat = 1)
//'  
//' out_2 <- multilogit_C(Y, X, n_sample = 2000, n_burn = 500, prior = "flat")
//' 
// [[Rcpp::export]]
List multilogit_C(
  arma::mat const &Y,
  arma::mat const &X, 
  size_t n_sample = 1000,
  size_t n_burn = 200,
  size_t n_sigma_check = 20,
  String prior = "flat",
  double step_size = 0.1,
  Nullable<NumericVector> prior_mean = R_NilValue,
  Nullable<NumericMatrix> prior_var = R_NilValue,
  Nullable<IntegerVector> reference_cat = R_NilValue,
  bool probs = true,
  bool progress = true
){
   
  size_t nCat = Y.n_cols;
  size_t nSub = Y.n_rows;
  
  size_t nPred = X.n_cols;
  
  arma::vec nObs(nSub);
  
  for(size_t i = 0; i < nSub; i++){
    
    nObs(i) = sum(Y.row(i));
  }
  
  /**
   * 
   * Priors
   */
  
  arma::rowvec beta_mean;
  arma::mat beta_var;
  size_t reference_beta = NA_INTEGER;
  
  if(prior_mean.isNotNull())
  {
    NumericVector mean2(prior_mean);
    
    beta_mean = Rcpp::as<arma::rowvec>(mean2);
    
  }
  
  else
  {
    
    
    beta_mean = rowvec(nPred, fill::zeros);
    
    
  }
  
  if(prior_var.isNotNull())
  {
    NumericMatrix v2(prior_var);
    
    beta_var = Rcpp::as<arma::mat>(v2);
    
    
  }
  
  
  else
  {
    
    beta_var = mat(nPred, nPred, fill::eye);
    
    if(prior == "normal")
      {
      Rcout << "Default prior variance is I_P\n";
      }
    
  }
  
  
  if(reference_cat.isNotNull())
  {
   
   IntegerVector temp_ref(reference_cat);
    
    reference_beta = temp_ref(0) - 1;
    
  }
  
  else
  {
    
    Rcout << "No reference category selected. Beta values for all categories generated.\n";
    
  }
  
  /**
   * Data structures
   */
  arma::cube beta_out(nPred, nCat, n_sample);
  arma::mat beta(nPred, nCat); beta.zeros();
  arma::mat candidate_sigmas(nPred, nCat); candidate_sigmas.fill(step_size);
  arma::mat acceptance(nPred, nCat); acceptance.zeros();
  arma::cube prob_out(nSub, nCat, n_sample);
  arma::mat prob(nSub, nCat);
  arma::vec phi(nSub, fill::ones);
  arma::mat phi_out(nSub, n_sample);
  double rate;
  size_t counter=0;
  
  // needed for random walk Metropolis:
  double proposal; 
  // double numerator;
  // double denominator;
  double lnumerator;
  double ldenominator;
  arma::vec betaWITHproposal(nPred);
  
  
  /**
   * MCMC
   */
  for(size_t iter=0; iter < (n_burn + n_sample); iter++) {
    
    if( progress == true && ((iter%1000) == 0 || (iter + 1) == n_burn + n_sample)) {
      
      Rcout << "iteration " << iter + 1 << " of " << n_sample + n_burn << "\n";
      
    }
    
    
    /**
     * Data augmentation
     */
    for(size_t i = 0; i < nSub; i++) { //can parallelize this in the future???
      
      rate = sum(exp(X.row(i) * beta)); //colptr for future speed up???
      
      phi(i) = (R::rgamma(nObs[i],1)) / rate;
      
    }
    
    /** 
     * Sample betas - right now it's a simple univariate random walk MH sampler. Lots of room to improve. 
     */
    
    if(prior == "flat") {
    
 
    for(size_t j = 0; j < nCat; j++) { //can parallelize this in the future
      
      if(reference_beta == j) {
        
        continue;
        
      }
      
      for(size_t k = 0; k < nPred; k++) {
        
        // generate proposal
        proposal = beta(k,j) + R::rnorm(0, candidate_sigmas(k,j));
        
        // proposal = beta(k,j) + R::rnorm(0, step_size);
        betaWITHproposal = beta.col(j); // can instead be cached and updated elementwise
        
        betaWITHproposal(k) = proposal;
        
        // calculate proposal's posterior
        // numerator = exp(dot(Y.col(j), X * betaWITHproposal) - dot(phi, exp(X * betaWITHproposal)));
        lnumerator = (dot(Y.col(j), X * betaWITHproposal) - dot(phi, exp(X * betaWITHproposal)));
        // calculate current value's posterior 
        //denominator = exp(dot(Y.col(j), X * beta.col(j)) - dot(phi, exp(X * beta.col(j))));
        ldenominator = (dot(Y.col(j), X * beta.col(j)) - dot(phi, exp(X * beta.col(j))));
        
        // accept/reject
        // if((numerator / denominator) > R::runif(0,1)) {
        if(( lnumerator - ldenominator) > log(R::runif(0,1))) {
          
          beta(k,j) = proposal;
          
          acceptance(k,j) += 1;
        }
      }
      
      
    }
    }
    
    if(prior == "normal") {
      
      for(size_t j = 0; j < nCat; j++) { //can parallelize this in the future
        
        if(reference_beta == j) {
          
          continue;
          
        }
        
        
        
        for(size_t k = 0;k < nPred; k++) {
          // generate proposal
          proposal = beta(k, j) + R::rnorm(0, candidate_sigmas(k, j));
          // proposal = beta(k,j) + R::rnorm(0, step_size);
          betaWITHproposal = beta.col(j); // can instead be cached and updated elementwise
          
          betaWITHproposal(k) = proposal;
          
          arma::rowvec row_beta_wp = trans(betaWITHproposal);
            
          arma::rowvec row_beta = trans(beta.col(j));
          
          
          
          // calculate proposal's posterior
//          numerator = exp(dot(Y.col(j), X * betaWITHproposal) - dot(phi, exp(X * betaWITHproposal))) * 
  //          dmvnrm_arma(row_beta_wp, beta_mean, beta_var);
          lnumerator = (dot(Y.col(j), X * betaWITHproposal) - dot(phi, exp(X * betaWITHproposal))) +
            dmvnrm_arma(row_beta_wp, beta_mean, beta_var,  true);
          
          // calculate current value's posterior 
          //denominator = exp(dot(Y.col(j), X * beta.col(j)) - dot(phi, exp(X * beta.col(j)))) *
            //dmvnrm_arma(row_beta, beta_mean, beta_var);
          ldenominator = (dot(Y.col(j), X * beta.col(j)) - dot(phi, exp(X * beta.col(j)))) +
            dmvnrm_arma(row_beta, beta_mean, beta_var, true);

          
          // accept/reject
          //if((numerator / denominator) > R::runif(0,1)) {
          if((lnumerator - ldenominator) > log(R::runif(0,1))) {
            
            beta(k,j) = proposal;
            
            acceptance(k,j) += 1;
            
          }
          
        }
        
      }
      
    }
    
    // Calculate probabilities 
    
    if (probs == true) {
      
      for(size_t i = 0; i < nSub; i++) {
        
        double temp_rate = sum(exp(X.row(i) * beta));
        
        prob.row(i) = exp(X.row(i) * beta) / temp_rate;
        
      }
      
    }
    
    /**
     * Adjust candidate sigmas
     **/
    
    if(iter < n_burn) {
      
      if((iter % n_sigma_check) == 0) {  
        
        for(size_t j = 0; j < nCat; j++) { //can parallelize this in the future???
          
          if(reference_beta == j) {
            
            continue;
          }
          for(size_t k = 0; k < nPred; k++) {
            // if accepting too many, then increase candidate sigma a lot
            if(acceptance(k, j) > (n_sigma_check * .40)) {
              
              candidate_sigmas(k, j) = candidate_sigmas(k, j) * 2;
            }
            // if accepting too few, then increase candidate sigma a little
            if(acceptance(k, j) < (n_sigma_check * .20)) {
              
              candidate_sigmas(k, j) = candidate_sigmas(k, j) * .9;
            }
            
            acceptance(k, j) = 0;
          }
        }
      }
    }
    
    
    /**
     * Write out 
     */
    
  if(probs == true) {    
    if(iter >= n_burn) {
      beta_out.slice(counter) = beta;
      prob_out.slice(counter) = prob;
      phi_out.col(counter) = phi;
      counter += 1;
    }
    
    
  } 
  
  else {
  if(iter >= n_burn) {
    beta_out.slice(counter) = beta;
    phi_out.col(counter) = phi;
    counter += 1;
  }
  
  }
  } // end iter loop, i.e. end MCMC
  
  if(probs == true) {
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
