#ifndef MNLOGIT_H
#define MNLOGIT_H

#include "RcppArmadillo.h"

using namespace Rcpp;

List multilogit_C(
        arma::mat const &Y,
        arma::mat const &X, 
        size_t n_sample = 1000,
        size_t n_burn = 200,
        size_t n_sigma_check = 20,
        String prior = "flat",
        double step_size = 0.1,
        Nullable<NumericVector> prior_mean = R_NilValue,
        Nullable<NumericMatrix> prior_var = R_NilValue,
        Nullable<IntegerVector> reference_cat = R_NilValue,
        bool probs = true,
        bool progress = true
);


#endif