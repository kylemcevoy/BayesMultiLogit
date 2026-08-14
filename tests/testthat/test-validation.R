make_data <- function() {
  list(
    Y = rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 0, 0)),
    X = cbind(1, c(-1, -0.5, 0.5, 1))
  )
}

test_that("all wrappers reject malformed common inputs", {
  dat <- make_data()
  wrappers <- list(
    function(Y, X, n_sample = 1, n_burn = 0, probs = FALSE, ...)
      multilogit(Y, X, n_sample, n_burn, probs = probs, progress = FALSE, ...),
    function(Y, X, n_sample = 1, n_burn = 0, probs = FALSE, ...)
      multilogit_ESS(Y, X, n_sample, n_burn, probs = probs, progress = FALSE, ...),
    function(Y, X, n_sample = 1, n_burn = 0, probs = FALSE, ...)
      multilogit_PG(Y, X, n_sample, n_burn, probs = probs, progress = FALSE, ...),
    function(Y, X, n_sample = 1, n_burn = 0, probs = FALSE, ...)
      multilogit_holmesheld(Y, X, n_sample = n_sample, n_burn = n_burn,
                            probs = probs, progress = FALSE, ...)
  )

  for (wrapper in wrappers) {
    expect_error(wrapper(dat$Y[, 1, drop = FALSE], dat$X), "at least two categories")
    expect_error(wrapper(dat$Y, dat$X[-1, ]), "same number of rows")
    expect_error(wrapper(dat$Y, replace(dat$X, 1, Inf)), "finite")
    expect_error(wrapper(dat$Y, dat$X, n_sample = 0), "positive integer")
    expect_error(wrapper(dat$Y, dat$X, probs = NA), "TRUE or FALSE")
  }
})

test_that("count and prior validation matches sampler requirements", {
  dat <- make_data()
  counts <- dat$Y
  counts[1, 1] <- 2

  expect_error(multilogit_holmesheld(counts, dat$X), "one-hot")
  expect_error(multilogit(dat$Y, dat$X, prior = "normal", prior_var = matrix(1, 2, 2)),
               "positive definite")
  expect_error(multilogit_ESS(dat$Y, dat$X, reference_cat = c(1, 2)),
               "single integer")
  expect_error(multilogit_PG(dat$Y, dat$X, n_burn = -1), "non-negative integer")
})
