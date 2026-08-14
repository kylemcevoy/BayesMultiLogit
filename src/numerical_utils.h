#ifndef BAYESMULTILOGIT_NUMERICAL_UTILS_H
#define BAYESMULTILOGIT_NUMERICAL_UTILS_H

#include "RcppArmadillo.h"

inline double row_log_sum_exp(const arma::rowvec& x) {
  const double maximum = x.max();
  return maximum + std::log(arma::accu(arma::exp(x - maximum)));
}

inline double row_log_sum_exp_excluding(const arma::rowvec& x,
                                        const arma::uword excluded) {
  double maximum = -std::numeric_limits<double>::infinity();
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    if (j != excluded && x(j) > maximum) maximum = x(j);
  }

  double total = 0.0;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    if (j != excluded) total += std::exp(x(j) - maximum);
  }
  return maximum + std::log(total);
}

inline arma::rowvec softmax(const arma::rowvec& x) {
  const double maximum = x.max();
  arma::rowvec result = arma::exp(x - maximum);
  return result / arma::accu(result);
}

#endif
