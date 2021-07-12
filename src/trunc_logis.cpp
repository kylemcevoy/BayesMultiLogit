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