\name{BayesMultiLogit-package}
\alias{BayesMultiLogit-package}
\alias{BayesMultiLogit}
\docType{package}
\title{
\packageTitle{BayesMultiLogit}
}
\description{
\packageDescription{BayesMultiLogit}
}
\details{
This package implements specific methods for performing Bayesian Multinomial Logistic Regression. The main functions included are wrappers multilogit, multilogit_PG, multilogit_holmesheld for the functions written in C++. The multilogit function implements a method explored in our paper "Bayesian Multinomial Logistic Regression for Numerous Categories". While the function multilogit_PG implements the Polya-Gamma method from the 2013 paper by Nicholas G. Polson and  James G. Scott and Jesse Windle, and the function multilogit_holmesheld implements the method proposed by Chris C. Holmes and Leonhard Held in their 2006 paper.

multilogit_holmesheld was programmed following the pseudo-code from the Holmes-Held paper referenced below.
multilogit_PG is adapted from code included in the BayesLogit thesis version by Jesse Windle: https://github.com/jwindle/BayesLogit-Thesis/.

Mistakes and bugs are our own.
}
\author{
\packageAuthor{BayesMultiLogit}

Maintainer: \packageMaintainer{BayesMultiLogit}
}
\references{
1. Holmes, C. C. and Held, L. (2006). Bayesian auxiliary variable models for binary and multinomial regression. Bayesian Analysis, 1(1):145 to 168.

2. Polson, N. G., Scott, J. G., and Windle, J. (2013). Bayesian inference for logistic models using Polya–gamma latent variables.Journal of the American Statistical Association, 108(504):1339 to 1349.


}
\keyword{ package }
\seealso{
\code{\link{multilogit}}, \code{\link{multilogit_ESS}},
\code{\link{multilogit_PG}}, and \code{\link{multilogit_holmesheld}}.
}
\examples{
Y <- matrix(0, nrow = 150, ncol = 3)
Y[1:50, 1] <- 1
Y[51:100, 2] <- 1
Y[101:150, 3] <- 1

X <- cbind(1, iris[ , 1:4])
X <- as.matrix(X)

out1 <- multilogit(Y, X, n_sample = 500, n_burn = 200, prior = "normal", progress = FALSE)

out2 <- multilogit_PG(Y, X, n_sample = 500, n_burn = 200, progress = FALSE)

out3 <- multilogit_holmesheld(Y, X, v = diag(ncol(X)), n_sample = 500, n_burn = 200, progress = FALSE)
}
