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

#ifndef HELLOPG_H
#define HELLOPG_H

#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;


colvec rpg(colvec shape, colvec scale);

SEXP helloPG(int n, double z);

#endif