.log_pdf_precision_cholesky <- function(Y, chol_Q) {
  sum(log(diag(chol_Q))) - 0.5 * sum((chol_Q %*% Y) ^ 2)
}

.sample_normal_precision_chol <- function(chol_Q) {
  if (is(chol_Q, 'triangularMatrix')) {
    as.vector(solve(chol_Q, rnorm(ncol(chol_Q))))
  } else if (is(chol_Q, 'CHMfactor')) {
    as.vector(solve(chol_Q, solve(
      chol_Q,
      rnorm(ncol(chol_Q)),
      system = 'Lt'
    ), system = 'Pt'))
  } else {
    as.vector(backsolve(chol_Q, rnorm(ncol(chol_Q))))
  }
}

.sample_normal_precision <- function(Q) {
  .sample_normal_precision_chol(Cholesky(Q, LDL = FALSE))
}

.sample_normal_covariance <- function(S) {
  as.vector(crossprod(chol(S), rnorm(ncol(S))))
}

.sample_normal_woodbury <- function(W) {
  .sample_normal_precision(W@A) + as.vector(W@X %*% .sample_normal_precision(W@B))
}

.dmvnorm <- function(x, mean, covariance, precision, log = FALSE) {
  if (missing(mean)) mean <- 0

  z <- x - mean
  output <- if (!missing(covariance)) {
    as.numeric(
      - 0.5 * (length(x) * log(2 * pi))
      - 0.5 * determinant(covariance, logarithm = TRUE)$modulus
      - 0.5 * crossprod(z, solve(covariance, z))
    )
  } else {
    as.numeric(
      - 0.5 * (length(x) * log(2 * pi))
      + 0.5 * determinant(precision, logarithm = TRUE)$modulus
      - 0.5 * crossprod(z, precison %*% z)
    )
  }

  if (log) output else exp(output)
}
