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
This package implements several MCMC methods for Bayesian multinomial logistic regression. The main R wrappers are multilogit_AMH, multilogit, multilogit_ESS, multilogit_PG, and multilogit_holmesheld. multilogit_AMH targets the multinomial likelihood directly with adaptive coordinate-wise Metropolis-Hastings. multilogit and multilogit_ESS implement methods explored in "Bayesian Multinomial Logistic Regression for Numerous Categories". multilogit_PG implements the Polya-Gamma method, and multilogit_holmesheld implements the method proposed by Chris C. Holmes and Leonhard Held.

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
\code{\link{multilogit_AMH}}, \code{\link{multilogit}}, \code{\link{multilogit_ESS}},
\code{\link{multilogit_PG}}, and \code{\link{multilogit_holmesheld}}.
}
\examples{
Y <- matrix(0, nrow = 150, ncol = 3)
Y[1:50, 1] <- 1
Y[51:100, 2] <- 1
Y[101:150, 3] <- 1

X <- cbind(1, iris[ , 1:4])
X <- as.matrix(X)

out0 <- multilogit_AMH(Y, X, n_sample = 20, n_burn = 10, progress = FALSE)

out1 <- multilogit(Y, X, n_sample = 500, n_burn = 200, prior = "normal", progress = FALSE)

out2 <- multilogit_PG(Y, X, n_sample = 500, n_burn = 200, progress = FALSE)

out3 <- multilogit_holmesheld(Y, X, v = diag(ncol(X)), n_sample = 500, n_burn = 200, progress = FALSE)
}
