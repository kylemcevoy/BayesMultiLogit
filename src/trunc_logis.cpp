// This code was written working off the pseudocode included in a 2006 paper
// by Chris C. Holmes and Leonhard Held and with corrections suggested by
// by Ralf van der Lans as found in the 2011 response:

// 1. Holmes,  C.  C.  and  Held,  L.  (2006). Bayesian  auxiliary  variable  models
// for  binary  and  multinomialregression. Bayesian Analysis, 1(1):145–168.

// 2. Holmes, C. C. and Held, L. (2011). Response to van der Lans.
// Bayesian Analysis, 6(2):357-358.

//truncated logistic sampler

#include <Rcpp.h>
using namespace Rcpp;

//bool right decides whether or not we are sampling from (0, inf) if TRUE, or (-inf, 0) if FALSE.

// [[Rcpp::export]]
double trunc_logis(double location, double scale, bool right) {
  double a{ 0 };
  double b{ 0 };
  
  if (right)
  {
    //R::plogis and qlogis require all arguments specified including lower tail and log format
    a = R::plogis(0, location, scale, 1, 0);
    b = 1;
    
  }
  else
  {
    a = 0;
    b = R::plogis(0, location, scale, 1, 0);
    
  }
  
  double p = R::runif(0,1);
  p = (1-p)*a + p*b; // = a + p(b-a)
  double logisDraw = R::qlogis(p, location, scale, 1, 0);

  
  
  return logisDraw;
}