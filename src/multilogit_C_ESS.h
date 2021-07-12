#ifndef MULTILOGIT_C_ESS_H
#define MULTILOGIT_C_ESS_H

#include "RcppArmadillo.h"

using namespace Rcpp;

List multilogit_C_ESS(
    NumericMatrix Y_,
    NumericMatrix X_, 
    size_t n_sample = 1000,
    size_t n_burn = 100,
    String prior = "normal",
    Nullable<NumericVector> prior_mean = R_NilValue,
    Nullable<NumericMatrix> prior_var = R_NilValue,
    Nullable<IntegerVector> reference_cat = R_NilValue,
    bool probs = true,
    bool progress = true
);


#endif