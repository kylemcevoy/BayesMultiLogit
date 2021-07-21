/**
 * ESS is a technique from the Murray, Adams, MacKay (2010) paper 
 * Elliptical Slice Sampling as found in the Journal of Machine Learning Research.
 * 
 * Pseudo-code from that paper was written into C++ by Jared Fisher and Kyle McEvoy.
 */

#include "dmvnrm_arma.h"
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

//' Multinomial Logistic Regression using Data Augmentation (Elliptical Slice Sampler)
//' 
//' @description This function implements a data augmentation method for
//' multinomial logisitic regression using MCMC. The data augmentation method is
//' outlined in a paper by Jared D. Fisher and Kyle R. McEvoy titled
//' "Bayesian Multinomial Logistic Regression for Numerous Categories",
//' currently a work in progress.
//' 
//' The sampler uses the elliptical slice sampling algorithm
//' from the Murray, Adams, Mackay (2010) paper Elliptical Slice Sampling in the Journal of
//' Machine Learning Research. This sampler requires a multivariate normal prior on the betas.
//'
//' @param Y An N by C numeric matrix where the ith row is a set of
//' indicators for observation i of N total observations giving which
//' of the C categories the observation is classified into.
//' @param X An N by P numeric matrix where N is the total number of
//' observations and P is the total number of predictor variables (including the intercept) the ith row
//' gives the values of the predictor variables for the ith outcome observation.
//' The first column of X should be an intercept column of 1s.
//' Non-intercept X columns should be centered and scaled by their standard deviations for best results. 
//' @param n_sample positive integer giving the number of samples to draw as
//'  output after burn-in.
//' @param n_burn non-negative integer giving the number of samples of burn-in
//'  before the chain output is saved.
//' @param prior_mean a vector of length equal to P the number of predictors.
//' Giving the mean for the normal prior on the coefficient vector.
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
//' out <- multilogit_C_ESS(Y, X, n_sample = 3000, n_burn = 1500, reference_cat = 1)
//'  
//' out_2 <- multilogit_C_ESS(Y, X, n_sample = 2000, n_burn = 500,
//'  prior_var = diag(x = 1, nrow = 5))
//' 
// [[Rcpp::export]]
List multilogit_C_ESS(
    arma::mat const &Y,
    arma::mat const &X, 
    size_t n_sample = 1000,
    size_t n_burn = 200,
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
    
    // beta_var = beta_var ;
    // Rcout << "Default prior variance is 10 * I_P\n";
    
      Rcout << "Default prior variance is I_P\n";
    
    
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
  arma::cube beta_out(nPred, nCat, n_sample);
  arma::mat beta(nPred, nCat); beta.zeros();
  // arma::mat candidate_sigmas(nPred,nCat); candidate_sigmas.fill(step_size);
  arma::mat acceptance(nPred, nCat); acceptance.zeros();
  arma::cube prob_out(nSub, nCat,n_sample);
  arma::mat prob(nSub, nCat);
  arma::vec phi(nSub, fill::ones);
  arma::mat phi_out(nSub, n_sample);
  double rate;
  size_t counter=0;
  
  // needed for ESS 

  // For debugging
  int while_iters; 
  // Setting initial values of beta to the linear probability model least squares...
  // @Jared do these lines still need to be in here? If so, should it be all non-reference beta columns?
  //if(reference_beta != 0) beta.col(0) = arma::inv(arma::trans(X)*X) * arma::trans(X) * Y.col(0);
  //if(reference_beta != 1) beta.col(1) = arma::inv(arma::trans(X)*X) * arma::trans(X) * Y.col(1);
  //if(reference_beta != 2) beta.col(2) = arma::inv(arma::trans(X)*X) * arma::trans(X) * Y.col(2);
  
  
  /**
   * MCMC
   */
  for(size_t iter=0; iter < (n_burn + n_sample); iter++){
    
    if( progress == true && ((iter % 1000) == 0 || (iter + 1) == n_burn + n_sample) ){
      Rcout << "iteration " << iter + 1 << " of " << n_sample + n_burn << "\n";
    }
    
    // Rcout << "iteration " << iter << std::endl;
    
    /**
     * Data augmentation
     */
    for(size_t i = 0; i < nSub;i++){ //can parallelize this in the future???
      rate = sum(exp(X.row(i) * beta)); //colptr for future speed up???
      // Rcout << "Xbeta"<< i<< " = "<< X.row(i)*beta << std::endl;
      // Rcout << "rate=" << rate << std::endl;
      
      phi(i) = (R::rgamma(nObs[i], 1)) / rate;
      
      // while denominator of probabilities is calculated, I figured we'd calculate the probs because it's not done elsewhere right now
      if(probs == true){
        prob.row(i) = exp(X.row(i) * beta) / rate;
        }
      
    }
    
    /** 
     * Sample betas - Elliptical Slice sampling. Assumes MVN prior on beta's.  
     */
      
      for(size_t j=0; j < nCat; j++){ //can parallelize this in the future???
        
        if(reference_beta == j){
          
          continue;
          
        }
        
        
        arma::vec var_weight = Rcpp::rnorm(nPred);
        
        // 1 // Choose ellipse
        nu = chol_var * var_weight;
        // Rcout << "nu=" << nu << std::endl;
        
        // 2 // Log-likelihood threshold
        // u = arma::randu(1);
        loglike_threshold =  dot(Y.col(j), X * beta.col(j)) - dot(phi, exp(X * beta.col(j))) + log(R::runif(0,1));
        // Rcout << "denominator=" << denominator << std::endl;
        
        // 3 // Draw an initial proposal, also defining a bracket
        // u = arma::randu(1);
        theta = R::runif(0, 2 * datum::pi); 
        theta_min = theta - 2 * datum::pi;
        theta_max = theta;
        // Rcout << "theta=" << theta << std::endl;
        
        // "slicing"
        proposal_likelihood = loglike_threshold +1;
        while_iters = 0;
        do{
          // Rcout << "while_iters=" << while_iters << std::endl;
          // Rcout << "iter,j = " << iter << " " << j << std::endl;
          
          // 4 // f' <- f*cos(theta) + nu*sin(theta)
          beta_proposal = beta.col(j) * cos(theta) + nu * sin(theta);
          // Rcout << "beta_proposal=" << beta_proposal << std::endl;
          
          // 5 //
          proposal_likelihood =  dot(Y.col(j), X * beta_proposal) - dot(phi, exp(X * beta_proposal));
          
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
      }
    
 
    if(probs == true){
      
      if(iter >= n_burn){
        beta_out.slice(counter) = beta;
        prob_out.slice(counter) = prob;
        phi_out.col(counter) = phi;
        counter += 1;
      }
      
      }
    
    else{
      
      if(iter >= n_burn){
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
