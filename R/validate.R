.validate_scalar_integer <- function(x, name, minimum) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(sprintf("%s must be a %s integer.", name,
                 if (minimum == 0L) "non-negative" else "positive"),
         call. = FALSE)
  }
  as.integer(x)
}

.validate_flag <- function(x, name) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

.validate_iteration_total <- function(n_sample, n_burn) {
  if (as.double(n_sample) + as.double(n_burn) > .Machine$integer.max) {
    stop("n_sample + n_burn must not exceed the maximum R integer.", call. = FALSE)
  }
}

.validate_data <- function(Y, X, indicator = FALSE) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)

  if (!is.numeric(Y) || !is.numeric(X)) {
    stop("X and Y must be numeric matrices.", call. = FALSE)
  }
  if (length(dim(Y)) != 2L || nrow(Y) < 1L || ncol(Y) < 2L) {
    stop("Y must have at least one row and at least two categories.", call. = FALSE)
  }
  if (length(dim(X)) != 2L || nrow(X) < 1L || ncol(X) < 1L) {
    stop("X must have at least one row and at least one predictor.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows.", call. = FALSE)
  }
  if (anyNA(X) || anyNA(Y) || any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("X and Y must contain only finite values.", call. = FALSE)
  }
  if (any(Y < 0) || any(Y != floor(Y))) {
    stop("Y must contain non-negative integer values.", call. = FALSE)
  }

  totals <- rowSums(Y)
  if (indicator) {
    if (any(Y > 1) || any(totals != 1)) {
      stop("For the Holmes-Held sampler, Y must be one-hot encoded with exactly one 1 per row.",
           call. = FALSE)
    }
  } else if (any(totals <= 0)) {
    stop("Every row of Y must have a positive total count.", call. = FALSE)
  }

  list(Y = Y, X = X)
}

.validate_covariance <- function(x, size, name) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !identical(dim(x), c(size, size)) ||
      anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("%s must be a finite %d by %d numeric matrix.", name, size, size),
         call. = FALSE)
  }
  if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must be symmetric.", name), call. = FALSE)
  }
  if (inherits(try(chol(x), silent = TRUE), "try-error")) {
    stop(sprintf("%s must be positive definite.", name), call. = FALSE)
  }
  x
}

.validate_reference <- function(reference_cat, n_cat) {
  if (is.null(reference_cat)) {
    return(NULL)
  }
  if (length(reference_cat) != 1L || !is.numeric(reference_cat) ||
      is.na(reference_cat) || !is.finite(reference_cat) ||
      reference_cat != floor(reference_cat) || reference_cat < 1L ||
      reference_cat > n_cat) {
    stop("reference_cat must be NULL or a single integer between 1 and the number of categories.",
         call. = FALSE)
  }
  as.integer(reference_cat)
}

.warn_missing_intercept <- function(X) {
  if (!isTRUE(all.equal(X[, 1L], rep(1, nrow(X)), tolerance = sqrt(.Machine$double.eps)))) {
    warning("The first column of X is expected to be an intercept column of ones.",
            call. = FALSE)
  }
}
