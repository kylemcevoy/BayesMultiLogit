// This file is part of the R package helloPG as found at https://github.com/jgscott/helloPG, 
// written by James Scott adapting code originally written by Jesse Windle in the package
// BayesLogit https://github.com/jwindle/BayesLogit. It is distributed under the GNU General
// Public License version 3 or later and without ANY warranty, implied or otherwise.
// Downloaded in 2021 by Kyle McEvoy.

// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see
// <http://www.gnu.org/licenses/>.

#include "RcppArmadillo.h"
#include <R_ext/Utils.h>
#include <iostream>
#include <exception>
#include "RNG.h"
#include "PolyaGamma.h"
// Rcpp::depends(RcppArmadillo)

using namespace Rcpp;
using namespace arma;

colvec rpg(colvec shape, colvec scale) {
// C++-only interface to PolyaGamma class
// draws random PG variates from arma::vectors of n's and psi's
  RNG r;
  PolyaGamma pg;
#ifdef USE_R
  GetRNGstate();
#endif
  int d = shape.n_elem;
  colvec result(d);
  for(int i=0; i<d; i++) {
    result[i] = pg.draw(shape(i), scale(i), r);
  }
#ifdef USE_R
  PutRNGstate();
#endif
  return result;
}

SEXP helloPG(int n, double z) {
  // returns n draws from PG(1,z)
  colvec pgscale(n, fill::ones);
  colvec pgshape(n);
  for(int i=0; i<n; i++) {
    pgshape[i] = z;
  }
  colvec out = rpg(pgscale, pgshape);

  return Rcpp::List::create(Rcpp::Named("draws")=out
			    );
}

