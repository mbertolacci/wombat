.make_a_sampler <- function(
  model,
  tuning = list(w = 0.25, max_evaluations = 100),
  Q_alpha = .make_Q_alpha(model)
) {
  if (!is.null(model[['a']])) {
    return(function(current, ...) {
      current$a <- model[['a']]
      current
    })
  }

  a_slice <- do.call(slice, tuning)

  function(current, warming_up) {
    output <- a_slice(current$a, function(a) {
      if (a <= 0 || a >= 1) return(-Inf)
      current2 <- current
      current2$a <- a

      log_prior_value <- log_prior(model, current2)
      if (!is.finite(log_prior_value)) return(log_prior_value)

      output <- log_prior_value + .log_pdf_precision_cholesky(
        current$alpha,
        chol(Q_alpha(current2))
      )
    }, learn = warming_up, include_n_evaluations = TRUE)

    log_trace(
      'a sampler took {output$n_evaluations} evaluations, w = {output$w}'
    )

    current$a <- output$sample
    current
  }
}
