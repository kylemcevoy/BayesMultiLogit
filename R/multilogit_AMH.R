#' Adaptive Metropolis-Hastings Multinomial Logistic Regression
#'
#' Fits a Bayesian multinomial logistic regression using coordinate-wise
#' random-walk Metropolis-Hastings proposals. Proposal standard deviations are
#' adapted during burn-in and fixed while posterior draws are retained. Unlike
#' the package's other samplers, this method targets the multinomial likelihood
#' directly and does not introduce augmentation variables.
#'
#' @param Y An N by C numeric matrix of non-negative category counts. Every row
#' must have a positive total.
#' @param X An N by P numeric design matrix. The first column should generally
#' be an intercept column of ones.
#' @param n_sample Positive integer number of retained posterior draws.
#' @param n_burn Non-negative integer number of burn-in iterations.
#' @param n_sigma_check Positive integer adaptation-window length.
#' @param step_size Finite positive initial random-walk standard deviation.
#' @param acceptance_lower Lower target acceptance rate in `[0, 1)`.
#' @param acceptance_upper Upper target acceptance rate in `(0, 1]`; must be
#' greater than `acceptance_lower`.
#' @param step_increase Multiplicative proposal-scale increase greater than one.
#' @param step_decrease Multiplicative proposal-scale decrease strictly between
#' zero and one.
#' @param prior_mean Optional finite numeric vector with one value per predictor.
#' Defaults to zero.
#' @param prior_var Optional positive-definite P by P prior covariance matrix.
#' Defaults to the identity matrix.
#' @param reference_cat Optional integer reference category. Its coefficients
#' are held at zero. Defaults to the last category.
#' @param probs Logical; if TRUE, return fitted probabilities for every retained
#' draw.
#' @param progress Logical; if TRUE, print periodic progress updates.
#'
#' @return A list with `posterior_coef`, a P by C by `n_sample` array;
#' `acceptance_rate`, the retained-sampling acceptance rate for each
#' coefficient; and `final_step_size`, the post-adaptation proposal scales. If
#' `probs = TRUE`, the list also contains `posterior_prob`, an N by C by
#' `n_sample` array. Diagnostics for the reference category are `NA`.
#'
#' @examples
#' X <- cbind(1, scale(iris[, 1:4]))
#' Y <- model.matrix(~ Species - 1, data = iris)
#' fit <- multilogit_AMH(
#'   Y, X, n_sample = 20, n_burn = 10,
#'   probs = TRUE, progress = FALSE
#' )
#'
#' @export
multilogit_AMH <- function(Y,
                           X,
                           n_sample = 1000L,
                           n_burn = 200L,
                           n_sigma_check = 20L,
                           step_size = 0.1,
                           acceptance_lower = 0.2,
                           acceptance_upper = 0.4,
                           step_increase = 2,
                           step_decrease = 0.9,
                           prior_mean = NULL,
                           prior_var = NULL,
                           reference_cat = NULL,
                           probs = FALSE,
                           progress = TRUE) {
  data <- .validate_data(Y, X)
  Y <- data$Y
  X <- data$X
  n_pred <- ncol(X)

  n_sample <- .validate_scalar_integer(n_sample, "n_sample", 1L)
  n_burn <- .validate_scalar_integer(n_burn, "n_burn", 0L)
  .validate_iteration_total(n_sample, n_burn)
  n_sigma_check <- .validate_scalar_integer(
    n_sigma_check, "n_sigma_check", 1L
  )
  probs <- .validate_flag(probs, "probs")
  progress <- .validate_flag(progress, "progress")

  validate_number <- function(value, name) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value)) {
      stop(sprintf("%s must be a finite number.", name), call. = FALSE)
    }
    as.numeric(value)
  }

  step_size <- validate_number(step_size, "step_size")
  acceptance_lower <- validate_number(acceptance_lower, "acceptance_lower")
  acceptance_upper <- validate_number(acceptance_upper, "acceptance_upper")
  step_increase <- validate_number(step_increase, "step_increase")
  step_decrease <- validate_number(step_decrease, "step_decrease")

  if (step_size <= 0) stop("step_size must be positive.", call. = FALSE)
  if (acceptance_lower < 0 || acceptance_upper > 1 ||
      acceptance_lower >= acceptance_upper) {
    stop("Acceptance bounds must satisfy 0 <= acceptance_lower < acceptance_upper <= 1.",
         call. = FALSE)
  }
  if (step_increase <= 1) {
    stop("step_increase must be greater than one.", call. = FALSE)
  }
  if (step_decrease <= 0 || step_decrease >= 1) {
    stop("step_decrease must be strictly between zero and one.", call. = FALSE)
  }

  if (is.null(prior_mean)) {
    prior_mean <- numeric(n_pred)
  } else if (!is.numeric(prior_mean) || length(prior_mean) != n_pred ||
             anyNA(prior_mean) || any(!is.finite(prior_mean))) {
    stop("prior_mean must be a finite numeric vector with one value per predictor.",
         call. = FALSE)
  } else {
    prior_mean <- as.numeric(prior_mean)
  }

  if (is.null(prior_var)) {
    prior_var <- diag(n_pred)
  } else {
    prior_var <- .validate_covariance(prior_var, n_pred, "prior_var")
  }

  if (is.null(reference_cat)) reference_cat <- ncol(Y)
  reference_cat <- .validate_reference(reference_cat, ncol(Y))
  .warn_missing_intercept(X)

  multilogit_AMH_C(
    Y = Y,
    X = X,
    n_sample = n_sample,
    n_burn = n_burn,
    n_sigma_check = n_sigma_check,
    step_size = step_size,
    acceptance_lower = acceptance_lower,
    acceptance_upper = acceptance_upper,
    step_increase = step_increase,
    step_decrease = step_decrease,
    prior_mean = prior_mean,
    prior_var = prior_var,
    reference_cat = reference_cat - 1L,
    probs = probs,
    progress = progress
  )
}
