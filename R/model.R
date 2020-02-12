.make_omega_parts <- function(
  process_model,
  measurement_model,
  variables = c('alpha', 'beta', 'eta')
) {
  n_parts <- c(
    alpha = ncol(process_model$H),
    beta = ncol(measurement_model$A),
    eta = ncol(process_model$Psi)
  )
  do.call(c, lapply(variables, function(variable) {
    rep(variable, n_parts[variable])
  }))
}

.make_omega <- function(
  process_model,
  measurement_model,
  variables = c('alpha', 'beta', 'eta')
) {
  omega_parts <- .make_omega_parts(process_model, measurement_model, variables)

  function(parameters) {
    omega <- rep(0, length(omega_parts))
    for (variable in variables) {
      omega[omega_parts == variable] <- parameters[[variable]]
    }
    omega
  }
}

.make_omega_unpack <- function(
  process_model,
  measurement_model,
  variables = c('alpha', 'beta', 'eta')
) {
  omega_parts <- .make_omega_parts(process_model, measurement_model, variables)

  function(parameters, omega) {
    for (variable in variables) {
      parameters[[variable]] <- omega[omega_parts == variable]
    }
    parameters
  }
}

.make_mu_omega <- function(
  process_model,
  measurement_model,
  variables = c('alpha', 'beta', 'eta')
) {
  function(params) {
    do.call(c, lapply(variables, function(variable) {
      if (variable == 'alpha') {
        process_model$alpha_prior_mean
      } else if (variable == 'beta') {
        measurement_model$beta_prior_mean
      } else if (variable == 'eta') {
        process_model$eta_prior_mean
      }
    }))
  }
}

.make_Q_omega <- function(
  process_model,
  measurement_model,
  variables = c('alpha', 'beta', 'eta')
) {
  if ('alpha' %in% variables) {
    Q_alpha <- .make_Q_alpha(process_model)
  }

  function(params) {
    bdiag(lapply(variables, function(variable) {
      if (variable == 'alpha') {
        Q_alpha(params)
      } else if (variable == 'beta') {
        measurement_model$beta_prior_precision
      } else if (variable == 'eta') {
        process_model$eta_prior_precision
      }
    }))
  }
}

.make_X_omega <- function(
  process_model,
  measurement_model,
  variables = c('alpha', 'beta', 'eta')
) {
  do.call(cbind, lapply(variables, function(variable) {
    if (variable == 'alpha') {
      measurement_model$C %*% process_model$H
    } else if (variable == 'beta') {
      measurement_model$A
    } else if (variable == 'eta') {
      measurement_model$C %*% process_model$Psi
    }
  }))
}

.make_chol_Q_omega_conditional <- function(
  process_model,
  measurement_model,
  variables = c('alpha', 'beta', 'eta'),
  X = .make_X_omega(process_model, measurement_model, variables),
  Xt_Q_epsilon_X = .make_Xt_Q_epsilon_X(X, measurement_model),
  Q_omega = .make_Q_omega(process_model, measurement_model, variables)
) {
  function(parameters) {
    chol(.fast_add(Q_omega(parameters), Xt_Q_epsilon_X(parameters)))
  }
}

.alpha_beta_eta_hat <- function(current, measurement_model, process_model) {
  Z2_tilde <- calculate(measurement_model, 'Z2_tilde', process_model)

  C <- measurement_model$C
  X <- cbind(
    C %*% process_model$H,
    measurement_model$A,
    C %*% process_model$Psi
  )

  n_alpha <- ncol(process_model$H)
  n_times <- n_alpha / length(process_model$regions)
  n_w <- length(process_model$regions)
  n_beta <- ncol(measurement_model$A)
  n_eta <- ncol(process_model$Psi)

  D <- Diagonal(x = current$gamma[measurement_model$attenuation_factor])
  Q_epsilon <- D %*% measurement_model$measurement_precision
  Q_alpha_t <- ar1_Q(n_times, current$a)
  Q_alpha_s <- Diagonal(n_w, x = current$w)
  Q_alpha <- kronecker(Q_alpha_t, Q_alpha_s)
  Q_omega <- bdiag(
    Q_alpha,
    measurement_model$beta_prior_precision,
    process_model$eta_prior_precision
  )

  .log_trace('Finding chol_Q_omega')
  chol_Q_omega_post <- chol(Q_omega + crossprod(chol(Q_epsilon) %*% X))
  omega_hat <- .chol_solve(chol_Q_omega_post, as.vector(
    crossprod(X, Q_epsilon %*% Z2_tilde)
  ))

  list(
    alpha = omega_hat[1 : n_alpha],
    beta = if (n_beta > 0) {
      omega_hat[(n_alpha + 1) : (n_alpha + n_beta)]
    } else {
      NULL
    },
    eta = as.vector(tail(omega_hat, n_eta)),
    chol_Q_omega = chol_Q_omega_post
  )
}
