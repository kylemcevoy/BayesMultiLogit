# BayesMultiLogit

BayesMultiLogit implements five MCMC approaches to Bayesian multinomial
logistic regression: direct adaptive Metropolis-Hastings, random-walk
Metropolis-Hastings with gamma augmentation, elliptical slice sampling with
gamma augmentation, Polya-Gamma augmentation, and Holmes-Held augmentation. See
"Bayesian Multinomial Logistic Regression for Numerous Categories" by Jared D. Fisher and Kyle 
R. McEvoy for further detail.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kylemcevoy/BayesMultiLogit")
```

## Input format

All wrapper functions accept an `N` by `C` response matrix `Y` and an `N` by
`P` design matrix `X`. `Y` may contain non-negative category counts for
`multilogit_AMH()`, `multilogit()`, `multilogit_ESS()`, and `multilogit_PG()`.
The Holmes-Held
sampler requires a one-hot response with exactly one 1 per row. The first
column of `X` should normally be an intercept column of ones.

```r
X <- cbind(1, scale(iris[, 1:4]))
Y <- model.matrix(~ Species - 1, data = iris)

fit <- multilogit_ESS(
  Y, X,
  n_sample = 1000,
  n_burn = 500,
  reference_cat = 1,
  probs = TRUE
)

dim(fit$posterior_coef) # P x C x n_sample
dim(fit$posterior_prob) # N x C x n_sample
```

For a direct posterior sampler without augmentation, use
`multilogit_AMH()`. It returns final proposal scales and retained-sampling
acceptance rates in addition to posterior draws.

The normal-prior samplers require positive-definite covariance matrices.
Using a reference category is recommended when identifiable category-specific
coefficients are required.

## Attribution

The data augmentation method underlying `multilogit()` and
`multilogit_ESS()` was developed by Jared D. Fisher and Kyle R. McEvoy.

The code for all functions using the Holmes-Held methods were written by us following the pseudo-code template provided in the paper:
"Bayesian Auxiliary Variable Models for Binary and Multinomial Regression" by Chris C. Holmes and Leonhard Held, Bayesian Analysis (2006).

The code for the Polya-Gamma functions was copied from the public repository https://github.com/jgscott/helloPG by James G. Scott who was adapting code
from the package https://github.com/jwindle/BayesLogit. The BayesLogit repository includes code written by Jesse Windle, Nicholas Polson, and James G. Scott.

The multivariate normal density function was copied from work by Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen on the Rcpp Gallery
located at https://gallery.rcpp.org/articles/dmvnorm_arma/.

## License

All of the work in this package falls under the Gnu General Public License Version 3 found at https://www.gnu.org/licenses/gpl-3.0.en.html.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
