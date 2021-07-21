#ifndef MULTILOGIT_C_ESS_H
#define MULTILOGIT_C_ESS_H

#include "RcppArmadillo.h"
#include "dmvnrm_arma.h"

using namespace Rcpp;

List multilogit_C_ESS(
        arma::mat const &Y,
        arma::mat const &X, 
        size_t n_sample = 1000,
        size_t n_burn = 200,
        Nullable<NumericVector> prior_mean = R_NilValue,
        Nullable<NumericMatrix> prior_var = R_NilValue,
        Nullable<IntegerVector> reference_cat = R_NilValue,
        bool probs = true,
        bool progress = true
);


#endif