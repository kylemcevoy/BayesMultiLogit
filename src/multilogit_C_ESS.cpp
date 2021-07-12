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


// [[Rcpp::export]]
List multilogit_C_ESS(
    NumericMatrix Y_,
    NumericMatrix X_, 
    size_t n_sample = 1000,
    size_t n_burn = 200,
    String prior = "normal",
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
    
    // beta_var = beta_var ;
    // Rcout << "Default prior variance is 10 * I_P\n";
    
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
  
  // For ESS
  // for random normal generation
  arma::mat chol_var = arma::chol(beta_var, "lower");
  arma::vec nu(nPred);
  // double u;
  double theta; // issue, so making this not "double"
  double theta_min;
  double theta_max;
  double proposal_likelihood;
  double loglike_threshold;
  // arma::vec denominator(nCat);
  arma::vec beta_proposal(nPred);
  // arma::vec u(1); //??? or how to scalar?
  
  /**
   * Data structures
   */
  arma::cube beta_out(nPred,nCat,n_sample);
  arma::mat beta(nPred,nCat); beta.zeros();
  // arma::mat candidate_sigmas(nPred,nCat); candidate_sigmas.fill(step_size);
  arma::mat acceptance(nPred,nCat); acceptance.zeros();
  arma::cube prob_out(nSub,nCat,n_sample);
  arma::mat prob(nSub,nCat);
  arma::vec phi(nSub, fill::ones);
  arma::mat phi_out(nSub,n_sample);
  double rate;
  size_t counter=0;
  
  // needed for ESS 

  // For debugging
  int while_iters; 
  // Setting initial values of beta to the linear probability model least squares...
  if(reference_beta != 0) beta.col(0) = arma::inv(arma::trans(X)*X) * arma::trans(X) * Y.col(0);
  if(reference_beta != 1) beta.col(1) = arma::inv(arma::trans(X)*X) * arma::trans(X) * Y.col(1);
  if(reference_beta != 2) beta.col(2) = arma::inv(arma::trans(X)*X) * arma::trans(X) * Y.col(2);
  
  
  /**
   * MCMC
   */
  for(size_t iter=0;iter<(n_burn+n_sample);iter++){
    
    if( progress == true && ((iter%1000) == 0 || (iter + 1) == n_burn + n_sample) ){
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
     * Sample betas - Elliptical Slice sampling. Assumes MVN prior on beta's.  
     */
  
    
    if(prior == "normal")
    {
      
      for(size_t j=0; j<nCat; j++){ //can parallelize this in the future???
        // Rcout << "j=" << j << std::endl;
        
        if(reference_beta == j){
          
          continue;
          
        }
        
        /**
         * ESS from Murray, Adams, MacKay (2010)
         */
        
        // I added this line because I'm worried that arma::randn() is not randomizing properly.

        arma::vec var_weight = Rcpp::rnorm(nPred);
        
        // 1 // Choose ellipse
        nu = chol_var * var_weight;
        // Rcout << "nu=" << nu << std::endl;
        
        // 2 // Log-likelihood threshold
        // u = arma::randu(1);
        loglike_threshold =  dot(Y.col(j),X*beta.col(j)) - dot(phi,exp(X*beta.col(j))) + log(R::runif(0,1));
        // Rcout << "denominator=" << denominator << std::endl;
        
        // 3 // Draw an initial proposal, also defining a bracket
        // u = arma::randu(1);
        theta = R::runif(0,2*datum::pi); 
        theta_min = theta - 2*datum::pi;
        theta_max = theta;
        // Rcout << "theta=" << theta << std::endl;
        
        // "slicing"
        proposal_likelihood = loglike_threshold +1;
        while_iters = 0;
        do{
          // Rcout << "while_iters=" << while_iters << std::endl;
          // Rcout << "iter,j = " << iter << " " << j << std::endl;
          
          // 4 // f' <- f*cos(theta) + nu*sin(theta)
          beta_proposal = beta.col(j)*cos(theta) + nu*sin(theta);
          // Rcout << "beta_proposal=" << beta_proposal << std::endl;
          
          // 5 //
          proposal_likelihood =  dot(Y.col(j),X*beta_proposal) - dot(phi,exp(X*beta_proposal));
          
          if (proposal_likelihood > loglike_threshold){
            beta.col(j) = beta_proposal;
            
            break;
            }
          // Rcout << "numerator=" << numerator << std::endl;
          
          
          // 7 // Else shrink the bracket and try again
          // 8 // if theta< 0, then theta_min = theta, else theta_max = theta
          if(theta < 0) theta_min = theta;
          else theta_max = theta; 
          // 9 // theta ~ U(theta_min,theta_max) ///////??? There is probably a faster implementation of this loop and logic. 
          theta = R::runif(theta_min,theta_max);
          // Rcout << "theta=" << theta << std::endl;
          while_iters += 1;
          
          // if(while_iters > 10000) numerator = denominator - 1;
        }
        while (proposal_likelihood <= loglike_threshold ) ;
        
        if(while_iters > 9000){
          Rcout << "used more than 9000 iterations in while loop instance" << std::endl;
        }
        // Rcout << "after numerator=" << numerator << std::endl;
        // Rcout << "after beta_proposal=" << beta_proposal << std::endl;
        // Rcout << "after while_iters=" << while_iters << std::endl;
        
        // 6 // Accept: return f'
      }
      
      
      
    }
    
 
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
    
    
    
    /**
     * Write out 
     */
    
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