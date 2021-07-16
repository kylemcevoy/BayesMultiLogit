// Code adapted by Kyle McEvoy from Rcpp Gallery code found at https://gallery.rcpp.org/articles/dmvnorm_arma/
// downloaded in 2020.
// 
// Originally written in 2013 by Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen.
// Slight modifications made by Kyle McEvoy 2021.
// This code falls under the Gnu GPL v2 license. 
// Code modified to only take a single vector of means.

#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

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