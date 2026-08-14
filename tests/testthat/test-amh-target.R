test_that("adaptive MH uses the multinomial likelihood and requested prior", {
  Y <- rbind(c(2, 0), c(0, 3))
  X <- matrix(1, nrow = 2, ncol = 1)
  prior_mean <- 0.3
  prior_var <- matrix(2, 1, 1)

  set.seed(918)
  proposal <- rnorm(1, 0, 0.2)
  eta_proposal <- rep(proposal, 2)
  log_denom_current <- rep(log(2), 2)
  log_denom_proposal <- log1p(exp(eta_proposal))
  log_likelihood_ratio <- sum(Y[, 1]) * proposal +
    sum(rowSums(Y) * (log_denom_current - log_denom_proposal))
  log_prior_ratio <- -0.5 * (
    (proposal - prior_mean)^2 / prior_var[1, 1] -
      (0 - prior_mean)^2 / prior_var[1, 1]
  )
  expected <- if (log(runif(1)) < log_likelihood_ratio + log_prior_ratio) {
    proposal
  } else {
    0
  }

  set.seed(918)
  out <- multilogit_AMH(
    Y, X, n_sample = 1, n_burn = 0, step_size = 0.2,
    prior_mean = prior_mean, prior_var = prior_var,
    reference_cat = 2, progress = FALSE
  )

  expect_equal(out$posterior_coef[1, 1, 1], expected, tolerance = 1e-14)
  expect_equal(out$posterior_coef[1, 2, 1], 0)
})
