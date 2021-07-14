// This code was written working off the pseudocode included in a 2006 paper
// by Chris C. Holmes and Leonhard Held and with corrections suggested by
// by Ralf van der Lans as found in the 2011 response:

// 1. Holmes,  C.  C.  and  Held,  L.  (2006). Bayesian  auxiliary  variable  models
// for  binary  and  multinomialregression. Bayesian Analysis, 1(1):145–168.

// 2. Holmes, C. C. and Held, L. (2011). Response to van der Lans.
// Bayesian Analysis, 6(2):357-358.

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
int right_interval(double U, double lambda)
{
  int OK{ 0 };
  
  double Z{ 1 };
  
  double X = exp(-0.5 * lambda);
  
  long j = 0;
  
  int exit_flag{ 0 };
  
  
  do
  {
    ++j;
    
    double d_j{ static_cast<double>(j) };
    
    double expon{ (d_j + 1) * (d_j + 1) - 1 };
    
    Z = Z - (d_j + 1) * (d_j + 1) * pow(X, expon);
    
    if (Z > U)
    {
      OK = 1;
      
      exit_flag = 1;
      
    }
    
    else
    {
      ++j;
      
      d_j = static_cast<double>(j);
      
      expon = (d_j + 1) * (d_j + 1) - 1;
      
      Z = Z + (d_j + 1) * (d_j + 1) * pow(X, expon);
      
      if (Z < U)
      {
        OK = 0;
        
        exit_flag = 1;
        
      }
      
      
    }
    
    
  } while (exit_flag == 0);
  
  return(OK);
}

// [[Rcpp::export]]
int left_interval(double U, double lambda)
{
  int OK{ 0 };
  
  int exit_flag{ 0 };
  
  double Z{ 1 };
  
  double X{ exp(-(M_PI * M_PI)) / (2 * lambda) };
  
  double H{ 0 };
  
  double K{ lambda / (M_PI * M_PI) };
  
  long j{ 0 };
  
  double log_U{ log(U) };
  
  H = (0.5 * log(2)) + (2.5 * log(M_PI)) - (2.5 * log(lambda)) - ((M_PI * M_PI) / (2 * lambda)) + (0.5 * lambda);
  
  do
  {
    
    ++j;
    
    double d_j{ static_cast<double>(j) };
    
    double expon{ d_j * d_j - 1};
    
    Z = Z - K * pow(X, expon);
    
    if ((H + log(Z)) > log_U)
    {
      
      OK = 1;
      
      exit_flag = 1;
      
    }
    
    else
    {
      ++j;
      
      d_j = static_cast<double>(j);
      
      expon = (d_j + 1) * (d_j + 1) - 1;
      
      Z = Z + (d_j + 1) * (d_j + 1) * pow(X, expon);
      
      if ((H + log(Z)) < log_U)
      {
        OK = 0;
        
        exit_flag = 1;
      }
      
    }
    
    
    
  } while (exit_flag == 0);
  
  return(OK);
  
}

// [[Rcpp::export]]
double lambda_sampler(double r)
{
  double W{ 0 };
  double U{ 0 };
  double U2{ 0 }; 
  double lambdaDraw{ 0 };
  int OK{ 0 };
  double R{ std::abs(r) };

  
  do
  {
    // Rcpp handles RNGState calls from R on export, no need to set a true random seed in function.
    // if code needs to be vectorized switch to Rcpp::rnorm/runif.
    
    W = R::rnorm(0, 1);
    
    U = R::runif(0, 1);
    
    W = W * W;
    
    W = 1 + ((W - sqrt(W * (4 * R + W))) / (2 * R));
    
    
    if (U <= (1 / (1 + W)))
    {
       
      lambdaDraw = R / W;
      
    }
    
    else
    {
      lambdaDraw = R * W;
      
    }
    
    U2 = R::runif(0, 1);
    
    //beginning the squeezing process on lambdaDraw
    
    if (lambdaDraw > (4 / 3))
    {
      
      OK =  right_interval(U2, lambdaDraw);
    }
    
    else
    {
      OK = left_interval(U2, lambdaDraw);
      
    }
    
    
  } while (OK == 0);
  
  return lambdaDraw;
  
  
}
