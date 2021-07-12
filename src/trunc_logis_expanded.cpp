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
