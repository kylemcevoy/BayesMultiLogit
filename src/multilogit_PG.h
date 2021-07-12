#ifndef MULTILOGIT_PG_H
#define MULTILOGIT_PG_H

#include "RcppArmadillo.h"
#include "helloPG.h"


using namespace Rcpp;
using namespace arma;

List multilogit_PG_C(NumericMatrix Y_,
               NumericMatrix X_, 
               size_t n_sample = 1000,
               size_t n_burn = 200);



#endif