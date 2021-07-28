#include "dmvnrm_arma.h"
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List metropolis_sampler(
    arma::mat const &Y,
    arma::mat const &X, 
    size_t n_sample = 1000,
    size_t n_burn = 200,
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
  arma::mat acceptance(nPred, nCat); acceptance.zeros();
  arma::cube prob_out(nSub, nCat, n_sample);
  arma::mat prob(nSub, nCat);
  arma::vec phi(nSub, fill::ones);
  arma::mat phi_out(nSub, n_sample);
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
  for(size_t iter=0; iter < (n_burn + n_sample); iter++) {
    
    if( progress == true && ((iter%1000) == 0 || (iter + 1) == n_burn + n_sample)) {
      
      Rcout << "iteration " << iter + 1 << " of " << n_sample + n_burn << "\n";
      
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
          proposal = beta(k,j) + R::rnorm(0, step_size);
          
          // proposal = beta(k,j) + R::rnorm(0, step_size);
          arma::mat beta_with_proposal = beta; // can instead be cached and updated elementwise
          
          beta_with_proposal(k,j) = proposal;
          
          // calculate proposal's posterior: prior is flat so we can ignore it.
          
          // likelihood numerator:
          
          arma::mat exp_X_beta = arma::exp(X * beta);
          
          arma::mat exp_X_beta_prop = arma::exp(X * beta_with_proposal);
          
          arma::vec likelihood_denom = arma::sum(exp_X_beta, 1);
          
          arma::vec likelihood_denom_prop = arma::sum(exp_X_beta_prop, 1);
          
          arma::vec likelihood(nSub, fill::zeros);
          
          arma::vec likelihood_prop(nSub, fill::zeros);
          
          
          for (size_t m = 0; m < nSub; m++) {
            
            likelihood_prop(m) = arma::dot(Y.row(m), exp_X_beta_prop.row(m)) / likelihood_denom_prop(m);
            
            likelihood(m) = arma::dot(Y.row(m), exp_X_beta.row(m)) / likelihood_denom(m);
            
          }
          
          double proposal_posterior = arma::sum(arma::log(likelihood_prop));
          
          
          double current_posterior = arma::sum(arma::log(likelihood));
          
          
          // accept/reject
          if((proposal_posterior - current_posterior) > log(R::runif(0,1))) {
            
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
          proposal = beta(k,j) + R::rnorm(0, step_size);
          
          // proposal = beta(k,j) + R::rnorm(0, step_size);
          arma::mat beta_with_proposal = beta; // can instead be cached and updated elementwise
          
          beta_with_proposal(k,j) = proposal;
          
          arma::rowvec beta_j = trans(beta.col(j));
          
          arma::rowvec beta_j_prop = trans(beta_with_proposal.col(j));
          
          // calculate proposal's posterior: prior is flat so we can ignore it.
          
          // likelihood numerator:
          
          arma::mat exp_X_beta = arma::exp(X * beta);
          
          arma::mat exp_X_beta_prop = arma::exp(X * beta_with_proposal);
          
          arma::vec likelihood_denom = arma::sum(exp_X_beta, 1);
          
          arma::vec likelihood_denom_prop = arma::sum(exp_X_beta_prop, 1);
          
          arma::vec likelihood(nSub, fill::zeros);
          
          arma::vec likelihood_prop(nSub, fill::zeros);
          
          
          for (size_t m = 0; m < nSub; m++) {
            
            likelihood_prop(m) = arma::dot(Y.row(m), exp_X_beta_prop.row(m)) / likelihood_denom_prop(m);
            
            likelihood(m) = arma::dot(Y.row(m), exp_X_beta.row(m)) / likelihood_denom(m);
            
          }
          
          double proposal_posterior = arma::prod(likelihood_prop) * dmvnrm_arma(beta_j_prop, beta_mean, beta_var);
          
          double current_posterior = arma::prod(likelihood) * dmvnrm_arma(beta_j, beta_mean, beta_var);
          
          // accept/reject
          if((proposal_posterior / current_posterior) > R::runif(0,1)) {
            
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
     * Write out 
     */
    
    if(probs == true) {    
      if(iter >= n_burn) {
        beta_out.slice(counter) = beta;
        prob_out.slice(counter) = prob;
        counter += 1;
      }
      
      
    } 
    
    else {
      if(iter >= n_burn) {
        beta_out.slice(counter) = beta;
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