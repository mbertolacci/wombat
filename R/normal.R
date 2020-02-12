.log_pdf_woodbury <- function(Y, D, chol_Q, Q_small, X) {
  chol_Q_small <- chol(Q_small)
  chol_D <- chol(D)
  log_det_Sigma_Y <- (
    2 * sum(log(diag(chol_Q)))
    - 2 * sum(log(diag(chol_D)))
    - 2 * sum(log(diag(chol_Q_small)))
  )
  a1 <- as.vector(chol_D %*% Y)
  a2 <- as.vector(D %*% Y)
  rhs <- crossprod(X, a2)
  sum_squares <- (
    sum(a1 ^ 2)
    - sum(solve(
      t(chol_Q),
      rhs
    ) ^ 2)
  )

  (
    - 0.5 * log_det_Sigma_Y
    - 0.5 * sum_squares
  )
}

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
