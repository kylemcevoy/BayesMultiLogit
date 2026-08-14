// Code adapted by Kyle McEvoy from Rcpp Gallery code found at https://gallery.rcpp.org/articles/dmvnorm_arma/
// downloaded in 2020.
// 
// Originally written in 2013 by Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen.
// Slight modifications made by Kyle McEvoy 2021.
// This code falls under the Gnu GPL v3 license. 
// Code modified to only take a single vector of means.

#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

static double const log2pi = std::log(2.0 * M_PI);

//'Multivariate Normal Density Function
//'
//' @description This function calculates densities for a mutlivariate normal distribution.
//' 
//' Code adapted by Kyle McEvoy from Rcpp Gallery code found at https://gallery.rcpp.org/articles/dmvnorm_arma/
//' downloaded in 2020.
//' 
//' Originally written in 2013 by Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen.
//' Slight modifications made by Kyle McEvoy 2021.
//' This code falls under the Gnu GPL v3 license. 
//' Code modified to only take a single vector of means.
//' 
//' @param x numeric vector The vector for which you want the density to be calculated.
//' @param mean numeric vector The vector of means for the Multivariate Normal Distribution.
//' should be the same length as x.
//' @param sigma numeric matrix Should be a square positive-definite covariance matrix for the desired
//' multivariate normal distribution. Each dimension should be equal to the length of x. 
//' @param logd logical; if TRUE, return the log density.
//' @return numeric The value of the pdf of the given Multivariate Normal distribution at the specified x value.
//' 
// [[Rcpp::export]]
double dmvnrm_arma(arma::rowvec const &x,  
                   arma::rowvec const &mean,  
                   arma::mat const &sigma, 
                   bool const logd = false) { 
  
  using arma::uword;
  
  uword const xdim = x.n_elem;
  
  double out;
  
  arma::mat const root = arma::chol(sigma);

  double const rootisum = -arma::sum(log(root.diag())),
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::vec z = arma::solve(arma::trimatl(root.t()), (x - mean).t(),
                            arma::solve_opts::fast);
  
  out = other_terms - 0.5 * arma::dot(z, z);     
  
  if (logd)
    return out;
  return exp(out);
}
