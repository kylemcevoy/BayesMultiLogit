make_smoke_data <- function() {
  list(
    Y = rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 0, 0),
              c(0, 1, 0), c(0, 0, 1)),
    X = cbind(1, seq(-1, 1, length.out = 6))
  )
}

check_output <- function(out, X, categories, draws, probabilities) {
  n <- nrow(X)
  p <- ncol(X)
  expect_equal(dim(out$posterior_coef), c(p, categories, draws))
  expect_true(all(is.finite(out$posterior_coef)))
  if (probabilities) {
    expect_equal(dim(out$posterior_prob), c(n, categories, draws))
    expect_true(all(is.finite(out$posterior_prob)))
    expect_equal(apply(out$posterior_prob, c(1, 3), sum),
                 matrix(1, n, draws), tolerance = 1e-12)
    for (draw in seq_len(draws)) {
      eta <- X %*% out$posterior_coef[, , draw]
      expected <- exp(eta - apply(eta, 1, max))
      expected <- expected / rowSums(expected)
      expect_equal(out$posterior_prob[, , draw], expected, tolerance = 1e-12)
    }
  }
}

test_that("samplers return finite arrays with normalized probabilities", {
  dat <- make_smoke_data()
  set.seed(42)

  outputs <- list(
    multilogit(dat$Y, dat$X, n_sample = 3, n_burn = 1, n_sigma_check = 1,
               reference_cat = 3, probs = TRUE, progress = FALSE),
    multilogit_ESS(dat$Y, dat$X, n_sample = 3, n_burn = 1,
                   prior_mean = c(0.5, -0.25), reference_cat = 3,
                   probs = TRUE, progress = FALSE),
    multilogit_PG(dat$Y, dat$X, n_sample = 3, n_burn = 1,
                  probs = TRUE, progress = FALSE),
    multilogit_holmesheld(dat$Y, dat$X, n_sample = 3, n_burn = 1,
                          probs = TRUE, progress = FALSE, fast = FALSE),
    multilogit_holmesheld(dat$Y, dat$X, n_sample = 3, n_burn = 1,
                          probs = TRUE, progress = FALSE, fast = TRUE)
  )

  for (out in outputs)
    check_output(out, dat$X, ncol(dat$Y), 3, TRUE)
})
