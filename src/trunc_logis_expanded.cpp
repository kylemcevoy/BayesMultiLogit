// This code was written working off the pseudocode included in a 2006 paper
// by Chris C. Holmes and Leonhard Held and with corrections suggested by
// by Ralf van der Lans as found in the 2011 response:

// 1. Holmes,  C.  C.  and  Held,  L.  (2006). Bayesian  auxiliary  variable  models
// for  binary  and  multinomialregression. Bayesian Analysis, 1(1):145–168.

// 2. Holmes, C. C. and Held, L. (2011). Response to van der Lans.
// Bayesian Analysis, 6(2):357-358.

#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec random_truncated_logistic(
    size_t n,
    double location,
    double scale, 
    double left,
    double right) {
  
  if(n < 1){
    stop("Error in random_truncated_logistic: n cannot be less than 1");
  }
  
  // currently outputting as arma::vec, hopefully it works! If not, will try numeric vector 
  arma::vec out;
  double A = R::plogis(left,location,scale,true,false); //thanks to "crch" the R library for this truncated sampling trick.  
  double B = R::plogis(right,location,scale,true,false);
  NumericVector p{ 0 }; 
  
  if(n > 1){
    p = Rcpp::runif(n);
    p = (1-p)*A + p*B; // = A + p(B-A)
    out = Rcpp::qlogis(p,location,scale);
  }
  if(n == 1){ //R:: is supposedly faster than Rcpp:: but can only return scalars.
    p = R::runif(0,1);
    p = (1-p)*A + p*B; // = A + p(B-A)
    out = Rcpp::qlogis(p,location,scale);
  }

  
  return out;
}
