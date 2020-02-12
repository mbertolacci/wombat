.make_gamma_sampler <- function(
  measurement_model,
  process_model,
  tuning = list(w = 1, max_evaluations = 100),
  variables = c('alpha', 'eta', 'beta'),
  X = .make_X_omega(
    process_model,
    measurement_model,
    variables
  ),
  X_c = .make_X_omega(
    process_model,
    measurement_model,
    setdiff(c('alpha', 'eta', 'beta'), variables)
  ),
  Q_epsilon = .make_Q_epsilon(measurement_model),
  chol_Q_omega_c_conditional = .make_chol_Q_omega_conditional(
    process_model,
    measurement_model,
    setdiff(c('alpha', 'eta', 'beta'), variables)
  ),
  Z2_tilde = calculate(measurement_model, 'Z2_tilde', process_model)
) {
  if (!is.null(measurement_model[['gamma']])) {
    return(function(current, ...) {
      current$gamma <- measurement_model[['gamma']]
      current
    })
  }

  n_gamma <- nlevels(measurement_model$attenuation_factor)
  omega <- .make_omega(process_model, measurement_model, variables)

  if (is.null(X_c)) {
    # Nothing has been marginalised out

    rg <- rgamma
    if ('lower' %in% names(measurement_model$gamma_prior)) {
      rg <- function(...) rtgamma(
        ...,
        lower = measurement_model$gamma_prior[['lower']],
        upper = measurement_model$gamma_prior[['upper']]
      )
    }

    function(current, ...) {
      Z_tilde <- as.vector(
        chol(measurement_model$measurement_precision)
        %*% (Z2_tilde - X %*% omega(current))
      )

      n <- tabulate(as.integer(measurement_model$attenuation_factor), n_gamma)
      sum_squares <- sapply(1 : n_gamma, function(i) {
        included <- as.integer(measurement_model$attenuation_factor) == i
        sum(Z_tilde[included] ^ 2)
      })

      current$gamma <- rg(
        n_gamma,
        shape = (
          measurement_model$gamma_prior[['shape']]
          + 0.5 * n
        ),
        rate = (
          measurement_model$gamma_prior[['rate']]
          + 0.5 * sum_squares
        )
      )

      current
    }
  } else {
    # Something has been marginalised out; use Woodbury
    gamma_slice <- lapply(seq_len(n_gamma), function(i) {
      do.call(slice, tuning)
    })

    function(current, warming_up) {
      Z_tilde <- as.vector(Z2_tilde - X %*% omega(current))

      current$gamma <- sapply(1 : n_gamma, function(i) {
        output <- gamma_slice[[i]](current$gamma[i], function(gamma_i) {
          current2 <- current
          current2$gamma[i] <- gamma_i

          log_prior_value <- log_prior(measurement_model, current2)
          if (!is.finite(log_prior_value)) return(log_prior_value)

          log_prior_value + .log_pdf_woodbury(
            Z_tilde,
            Q_epsilon(current2),
            chol_Q_omega_c_conditional(current2),
            process_model$eta_prior_precision,
            X_c
          )
        }, learn = warming_up, include_n_evaluations = TRUE)

        .log_trace(
          'gamma[%d] took %d evaluations, w = %f',
          i, output$n_evaluations, output$w
        )

        output$sample
      })

      current
    }
  }
}
