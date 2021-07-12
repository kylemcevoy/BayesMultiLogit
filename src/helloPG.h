#ifndef HELLOPG_H
#define HELLOPG_H

#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;


colvec rpg(colvec shape, colvec scale);

SEXP helloPG(int n, double z);

#endif