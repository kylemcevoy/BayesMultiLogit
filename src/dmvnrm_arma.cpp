// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//Code borrowed from Rcpp Gallery https://gallery.rcpp.org/articles/dmvnorm_arma/. Slightly modified so only takes single x sample.
//This code falls under the Gnu Public License version 2 license. Originally written in 2013 by Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen

static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
double dmvnrm_arma(arma::rowvec const &x,  
                   arma::rowvec const &mean,  
                   arma::mat const &sigma, 
                   bool const logd = false) { 
  using arma::uword;
  uword const xdim = x.n_elem;
  double out;
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  z      = (x - mean) * rooti;    
  
  out = other_terms - 0.5 * arma::dot(z, z);     
  
  if (logd)
    return out;
  return exp(out);
}