% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dmvnrm_arma}
\alias{dmvnrm_arma}
\title{Multivariate Normal Density Function}
\usage{
dmvnrm_arma(x, mean, sigma, logd = FALSE)
}
\arguments{
\item{x}{numeric vector The vector for which you want the density to be calculated.}

\item{mean}{numeric vector The vector of means for the Multivariate Normal Distribution.
should be the same length as x.}

\item{sigma}{numeric matrix Should be a square positive semi-definite Covariance matrix for the desired
multivariate normal distribution. Each dimension should be equal to the length of x.}
}
\value{
numeric The value of the pdf of the given Multivariate Normal distribution at the specified x value.
}
\description{
This function calculates densities for a mutlivariate normal distribution.

Code adapted by Kyle McEvoy from Rcpp Gallery code found at https://gallery.rcpp.org/articles/dmvnorm_arma/
downloaded in 2020.

Originally written in 2013 by Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen.
Slight modifications made by Kyle McEvoy 2021.
This code falls under the Gnu GPL v3 license. 
Code modified to only take a single vector of means.
}
