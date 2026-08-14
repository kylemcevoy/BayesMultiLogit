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
    multilogit_AMH(dat$Y, dat$X, n_sample = 3, n_burn = 2,
                   n_sigma_check = 1, reference_cat = 2,
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

test_that("adaptive MH returns coherent diagnostics and fixes its reference", {
  dat <- make_smoke_data()
  set.seed(123)
  out <- multilogit_AMH(
    dat$Y, dat$X, n_sample = 5, n_burn = 3, n_sigma_check = 2,
    prior_mean = c(0.25, -0.1), prior_var = diag(c(2, 0.5)),
    reference_cat = 2, progress = FALSE
  )

  expect_true(all(out$posterior_coef[, 2, ] == 0))
  expect_equal(dim(out$acceptance_rate), c(ncol(dat$X), ncol(dat$Y)))
  expect_equal(dim(out$final_step_size), c(ncol(dat$X), ncol(dat$Y)))
  expect_true(all(is.na(out$acceptance_rate[, 2])))
  expect_true(all(is.na(out$final_step_size[, 2])))
  expect_true(all(out$acceptance_rate[, -2] >= 0 &
                  out$acceptance_rate[, -2] <= 1))
  expect_true(all(out$final_step_size[, -2] > 0))
})
