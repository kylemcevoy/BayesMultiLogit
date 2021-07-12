// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "dmvnrm_arma.h"
#include "RcppArmadillo.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//' Multinomial Logistic Regression using Data Augmentation
//' 
//' @description This function implements a data augmentation method for
//' multinomial logisitic regression using MCMC. The data augmentation method is
//' outlined in a paper by Jared D. Fisher and Kyle R. McEvoy titled
//' "Bayesian Multinomial Logistic Regression for Numerous Categories",
//' currently a work in progress. The sampler uses a Random-Walk Metropolis
//' algorithm for choosing proposals and acceptances.
//' 
//'
//' @param Y_ An N by C numeric matrix where the ith row is a set of
//' indicators for observation i of N total observations giving which
//' of the C categories the observation is classified into.
//' @param X_ An N by P numeric matrix where N is the total number of
//' observations and P is the total number of predictor variables the ith row
//' gives the values of the predictor variables for the ith outcome observation. 
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
//' out <- multilogit_C(Y, X, n_sample = 3000, n_burn = 1500, n_sigma_check = 20,
//'  step_size = 0.1, prior = "normal", reference_cat = 1)
//'  
//' out_2 <- multilogit_C(Y, X, n_sample = 2000, n_burn = 500, prior = "flat")
// [[Rcpp::export]]
List multilogit_C(
  NumericMatrix Y_,
  NumericMatrix X_, 
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
   
   
  /**
   * Y 
   */
  
  size_t nCat = Y_.ncol();
  size_t nSub = Y_.nrow();
  
  
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  
  arma::vec nObs(nSub);
  // nObs = cumsum(Y,1); doesn't work???
  for(size_t i=0;i<nSub;i++){
    nObs(i) = sum(Y.row(i));
  }
  
  
  /**
   * X
   */
  size_t nPred = X_.ncol();
  
  arma::mat X = Rcpp::as<arma::mat>(X_);
  
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
  arma::cube beta_out(nPred,nCat,n_sample);
  arma::mat beta(nPred,nCat); beta.zeros();
  arma::mat candidate_sigmas(nPred,nCat); candidate_sigmas.fill(step_size);
  arma::mat acceptance(nPred,nCat); acceptance.zeros();
  arma::cube prob_out(nSub,nCat,n_sample);
  arma::mat prob(nSub,nCat);
  arma::vec phi(nSub, fill::ones);
  arma::mat phi_out(nSub,n_sample);
  double rate;
  size_t counter=0;
  
  // needed for random walk Metropolis:
  double proposal; 
  double numerator;
  double denominator;
  arma::vec betaWITHproposal(nPred);
  
  
  /**
   * MCMC
   */
  for(size_t iter=0;iter<(n_burn+n_sample);iter++){
    
    if( progress == true && ((iter%1000) == 0 || (iter + 1) == n_burn + n_sample)){
      Rcout << "iteration " << iter + 1 << " of " << n_sample + n_burn << "\n";
    }
    
    // Rcout << "iteration " << iter << std::endl;
    
    /**
     * Data augmentation
     */
    for(size_t i=0;i<nSub;i++){ //can parallelize this in the future???
      rate = sum(exp(X.row(i)*beta)); //colptr for future speed up???
      // Rcout << "Xbeta"<< i<< " = "<< X.row(i)*beta << std::endl;
      // Rcout << "rate=" << rate << std::endl;
      
      phi(i) = (R::rgamma(nObs[i],1))/rate;
      
      // while denominator of probabilities is calculated, I figured we'd calculate the probs because it's not done elsewhere right now
      if(probs == true){
        prob.row(i) = exp(X.row(i)*beta)/rate;
        }
      
    }
    // Rcout << "last rate=" << rate << std::endl;
    
    /** 
     * Sample betas - right now it's a too-simple univariate random walk MH sampler. Lots of room to improve. 
     */
    
    if(prior == "flat")
    {
    
 
    for(size_t j=0; j<nCat; j++){ //can parallelize this in the future???
      // can add option here if a reference category is desired s.t. betas are identified
      // beta_sampler(Y.col(j),X,beta.col(j),phi); //add prior info in the future??? currently just a flat prior
      if(reference_beta == j){
        
        continue;
        
      }
      
      for(size_t k=0;k<nPred;k++){
        // generate proposal
        proposal = beta(k,j) + R::rnorm(0, candidate_sigmas(k,j));
        // proposal = beta(k,j) + R::rnorm(0, step_size);
        betaWITHproposal = beta.col(j); // can instead be cached and updated elementwise
        betaWITHproposal(k) = proposal;
        
        // calculate proposal's posterior
        numerator = exp( dot(Y.col(j),X*betaWITHproposal) - dot(phi,exp(X*betaWITHproposal)));
        
        // calculate current value's posterior 
        denominator = exp( dot(Y.col(j),X*beta.col(j)) - dot(phi,exp(X*beta.col(j))));
        
        // accept/reject
        if((numerator/denominator) > R::runif(0,1)){
          beta(k,j) = proposal;
          acceptance(k,j) += 1;
        }
      }
      
      
    }
    }
    
    if(prior == "normal")
    {
      
      for(size_t j=0; j<nCat; j++){ //can parallelize this in the future???
        
        if(reference_beta == j){
          
          continue;
          
        }
        
        
        
        for(size_t k=0;k<nPred;k++){
          // generate proposal
          proposal = beta(k,j) + R::rnorm(0, candidate_sigmas(k,j));
          // proposal = beta(k,j) + R::rnorm(0, step_size);
          betaWITHproposal = beta.col(j); // can instead be cached and updated elementwise
          betaWITHproposal(k) = proposal;
          
          arma::rowvec row_beta_wp = trans(betaWITHproposal);
            
          arma::rowvec row_beta = trans(beta.col(j));
          
          
          
          // calculate proposal's posterior
          numerator = exp( dot(Y.col(j),X*betaWITHproposal) - dot(phi,exp(X*betaWITHproposal))) * dmvnrm_arma(row_beta_wp, beta_mean, beta_var)  ;
          
          // calculate current value's posterior 
          denominator = exp( dot(Y.col(j),X*beta.col(j)) - dot(phi,exp(X*beta.col(j)))) * dmvnrm_arma(row_beta, beta_mean, beta_var);
          
          // accept/reject
          if((numerator/denominator) > R::runif(0,1)){
            beta(k,j) = proposal;
            acceptance(k,j) += 1;
          }
        }
        
        
      }
    }
    
    /**
     * Adjust candidate sigmas
     **/
    
    if(iter<n_burn){
      if((iter%n_sigma_check)==0){  
        for(size_t j=0; j<nCat; j++){ //can parallelize this in the future???
          if(reference_beta == j){
            continue;
          }
          for(size_t k=0;k<nPred;k++){
            // if accepting too many, then increase candidate sigma a lot
            if(acceptance(k,j) > (n_sigma_check *.40)){
              candidate_sigmas(k,j) = candidate_sigmas(k,j)*2;
            }
            // if accepting too few, then increase candidate sigma a little
            if(acceptance(k,j) < (n_sigma_check *.20)){
              candidate_sigmas(k,j) = candidate_sigmas(k,j)*.9;
            }
            acceptance(k,j) = 0;
          }
        }
      }
    }
    
    
    /**
     * Write out 
     */
    
  if(probs == true){    
    if(iter>=n_burn){
      beta_out.slice(counter) = beta;
      prob_out.slice(counter) = prob;
      phi_out.col(counter) = phi;
      counter += 1;
    }
    
    
  } 
  
  else{
  if(iter>=n_burn){
    beta_out.slice(counter) = beta;
    phi_out.col(counter) = phi;
    counter += 1;
  }
  
  }
  } // end iter loop, i.e. end MCMC
  
  if(probs == true){
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



/* just plugging straight into main function.
arma::vec beta_sampler(arma::vec& yj, arma::mat& X, arma::vec& betaj, arma::vec& phi){
  // This is a simple random walk metropolis sampler, only as a placeholder for now. 
  double proposal; 
  double numerator;
  double denominator;
  arma::vec betaWITHproposal(betaj.n_elem);
  
  for(size_t k=0;k<betaj.n_elem;k++){
    // generate proposal
    proposal = betaj(k) + R::rnorm(0,0.1);
    betaWITHproposal = betaj; // can instead be cached and updated elementwise
    betaWITHproposal(k) = proposal;
    
    // calculate proposal's posterior
    numerator = exp( dot(yj,X*betaWITHproposal) - dot(phi,exp(X*betaWITHproposal)));
    
    // calculate current value's posterior 
    denominator = exp( dot(yj,X*betaj) - dot(phi,exp(X*betaj)));
    
    // accept/reject
    if((numerator/denominator) > R::runif(0,1)){
      betaj(k) = proposal;
    }
  }
  
  return(betaj);
} 
*/