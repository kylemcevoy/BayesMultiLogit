#ifndef MULTILOGIT_HH_INV_C_H
#define MULTILOGIT_HH_INV_C_H

#include "RcppArmadillo.h"
#include "lambda_sampler.h"
#include "trunc_logis.h"

using namespace Rcpp;
using namespace arma;

List multilogit_hh_inv_C(
    NumericMatrix Y_, 
    NumericMatrix X_, 
    NumericMatrix v_, 
    size_t n_sample = 1000,
    size_t n_burn = 200,
    bool probs = true,
    bool progress = true
);

#endif