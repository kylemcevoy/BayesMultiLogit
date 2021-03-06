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

#include "RRNG.h"

//////////////////////////////////////////////////////////////////////
		      // R Random Variates //
//////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------
// Distributions with one parameter.

#define ONEP(NAME, CALL, P1)			\
  double BasicRNG::NAME(double P1)	\
  {						\
    return CALL (P1);				\
  }						\

ONEP(expon_mean, rexp  , mean)
ONEP(chisq     , rchisq, df  )

#undef ONEP

//--------------------------------------------------------------------
// Distributions with two parameters.

#define TWOP(NAME, CALL, P1, P2)			\
  double BasicRNG::NAME(double P1, double P2)	\
  {							\
    return CALL (P1, P2);				\
  }							\

TWOP(gamma_scale, rgamma, shape, scale)
TWOP(norm      , rnorm , mean , sd  )
TWOP(flat      , runif , a    , b   )
TWOP(beta      , rbeta , a    , b   )

// x ~ Gamma(shape=a, scale=b)
// x ~ x^{a-1} exp(x / b).

#undef TWOP

//--------------------------------------------------------------------
			    // Uniform //

double BasicRNG::unif()
{
  return unif_rand();
} // unif

//--------------------------------------------------------------------
			  // Exponential //
double BasicRNG::expon_rate(double rate)
{
  return expon_mean(1.0 / rate);
}

//--------------------------------------------------------------------
			    // Normal //

double BasicRNG::norm(double sd)
{
  return rnorm(0, sd);
} // norm

//--------------------------------------------------------------------
			   // gamma_rate //

double BasicRNG::gamma_rate(double shape, double rate)
{
  return gamma_scale(shape, 1.0 / rate);
}

//--------------------------------------------------------------------
			   // Inv-Gamma //

// a = shape, b = scale
// x ~ IG(shape, scale) ~ x^{-a-1} exp(b / x).
// => 1/x ~ Ga(shape, scale*=1/scale).

double BasicRNG::igamma(double shape, double scale)
{
  return 1.0/rgamma(shape, 1.0 / scale);
} // igamma

////////////////////////////////////////////////////////////////////////////////

double BasicRNG::p_norm(double x, int use_log)
{
  return pnorm(x, 0.0, 1.0, 1, use_log);
}

double BasicRNG::p_gamma_rate(double x, double shape, double rate, int use_log)
{
  double scale = 1.0 / rate;
  return pgamma(x, shape, scale, 1, use_log);
}

////////////////////////////////////////////////////////////////////////////////

double BasicRNG::Gamma (double x, int use_log)
{
  double y = lgammafn(x);
  if (!use_log) y = exp(y);
  return y;
}

////////////////////////////////////////////////////////////////////////////////

double BasicRNG::d_beta(double x, double a, double b)
{
  return dbeta(x, a, b, false);
}

////////////////////////////////////////////////////////////////////////////////


