#ifndef DMVNRM_ARMA_H
#define DMVNRM_ARMA_H

#include "RcppArmadillo.h"

double dmvnrm_arma(arma::rowvec const &x,  
                   arma::rowvec const &mean,  
                   arma::mat const &sigma, 
                   bool const logd = false);


#endif