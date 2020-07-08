.make_w_sampler <- function(model) {
  n_times <- ncol(model$H) / length(model$regions)
  n_w <- length(model$regions)

  rg <- rgamma
  if ('lower' %in% names(model$w_prior)) {
    rg <- function(...) rtgamma(
      ...,
      lower = model$w_prior[['lower']],
      upper = model$w_prior[['upper']]
    )
  }

  function(current) {
    chol_Q_alpha_t <- chol(ar1_Q(n_times, current$a))
    alpha_matrix <- matrix(current$alpha, nrow = n_times, byrow = TRUE)

    current$w <- rg(
      n_w,
      shape = model$w_prior[['shape']] + 0.5 * n_times,
      rate = (
        model$w_prior[['rate']]
        + 0.5 * colSums((chol_Q_alpha_t %*% alpha_matrix) ^ 2)
      )
    )
    if (!is.null(model[['w']])) {
      has_preset <- !is.na(model$w)
      current$w[has_preset] <- model$w[has_preset]
    }
    current
  }
}
