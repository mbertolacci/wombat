#' @export
inversion_map_em <- function(
  measurement_model,
  process_model,
  max_iterations = 100,
  objective_tolerance = sqrt(.Machine$double.eps),
  step_tolerance = sqrt(.Machine$double.eps),
  start = NULL,
  trace = 0,
  nlm_args = list(),
  ...
) {
  n_alpha <- ncol(process_model$H)
  n_beta <- ncol(measurement_model$A)
  n_regions <- length(process_model$regions)
  n_times <- n_alpha / length(process_model$regions)
  n_w <- length(process_model$regions)
  n_gamma <- nlevels(measurement_model$attenuation_factor)

  find_a <- is.null(process_model[['a']])
  find_w <- is.null(process_model[['w']])
  find_gamma <- is.null(measurement_model[['gamma']])
  get_a <- function(theta) {
    if (find_a) .expit_transform(theta[1]) else process_model$a
  }
  get_w <- function(theta) {
    if (find_w) {
      .softplus_transform(tail(theta, -1))
    } else process_model$w
  }

  .log_debug('Precomputing various quantities')
  attenuation_index <- as.integer(measurement_model$attenuation_factor)
  Z <- measurement_model$soundings$xco2 - as.vector(measurement_model$C %*% process_model$control$xco2)
  C <- measurement_model$C
  X <- cbind(
    C %*% process_model$H,
    measurement_model$A,
    C %*% process_model$Psi
  )
  Xt_Q_epsilon_X <- .make_Xt_Q_epsilon_X(X, measurement_model)
  Q_epsilon <- .make_Q_epsilon(measurement_model)
  Q_alpha <- .make_Q_alpha(process_model)
  Q_omega <- .make_Q_omega(
    process_model,
    measurement_model,
    c('alpha', 'beta', 'eta')
  )

  .log_debug('Setting start values')
  start <- .extend_list(
    c(
      generate(process_model)[c('a', 'w')],
      generate(measurement_model)[c('gamma')]
    ),
    start,
    .remove_nulls(process_model[c('a', 'w')]),
    .remove_nulls(process_model[c('gamma')])
  )[c('a', 'w', 'gamma')]

  current <- start
  previous <- start
  log_posterior <- -Inf
  code <- 3

  .log_debug('Starting EM algorithm')
  for (iteration in 1 : max_iterations) {
    .log_trace('[%d] Finding sufficient statistics', iteration, max_iterations)
    chol_Q_omega_conditional <- chol(
      Xt_Q_epsilon_X(current) + Q_omega(current)
    )
    mu_omega_conditional <- .chol_solve(
      chol_Q_omega_conditional,
      crossprod(X, Q_epsilon(current) %*% Z)
    )

    E_omega_t_omega <- (
      chol2inv(chol_Q_omega_conditional)
      + tcrossprod(mu_omega_conditional)
    )

    if (find_a || find_w) {
      .log_trace('[%d] Optimising for a and w', iteration, max_iterations)
      T <- E_omega_t_omega[1 : n_alpha, 1 : n_alpha]
      output <- suppressWarnings(do.call(nlm, c(list(function(theta) {
        a <- get_a(theta)
        w <- get_w(theta)
        Q_alpha_current <- Q_alpha(list(a = a, w = w))
        -(
          0.5 * determinant(Q_alpha_current, logarithm = TRUE)$modulus
          - 0.5 * sum(diag(Q_alpha_current %*% T))
          + log_prior(process_model, a, w)
        )
      }, c(
        if (find_a) .logit_transform(current$a) else NULL,
        if (find_w) .inv_softplus_transform(current$w) else NULL
      )), nlm_args)))
      # stopifnot(output$code <= 3)
      current$a <- get_a(output$estimate)
      current$w <- get_w(output$estimate)
    }

    if (find_gamma) {
      .log_trace('[%d] Optimising for gamma', iteration)
      Y_tilde <- as.vector(sqrt(measurement_model$measurement_precision) %*% (Z - X %*% mu_omega_conditional))

      current$gamma <- sapply(1 : n_gamma, function(gamma_index) {
        indices <- attenuation_index == gamma_index
        n <- sum(indices)

        Xt_Q_epsilon0_X_part <- Xt_Q_epsilon_X(parts = gamma_index)

        shape <- measurement_model$gamma_prior[1] + n / 2
        rate <- measurement_model$gamma_prior[2] + 0.5 * (
          sum(Y_tilde[indices] ^ 2)
          - as.vector(crossprod(mu_omega_conditional, Xt_Q_epsilon0_X_part %*% mu_omega_conditional))
          + sum(diag(E_omega_t_omega %*% Xt_Q_epsilon0_X_part))
        )
        (shape - 1) / rate
      })
    }

    .log_trace('[%d] Calculating current log posterior', iteration)
    Q_omega_current <- Q_omega(current)
    log_posterior_previous <- log_posterior
    log_posterior <- (
      .log_pdf_woodbury(
        Z,
        Q_epsilon(current),
        chol(Q_omega_current + Xt_Q_epsilon_X(current)),
        Q_omega_current,
        X
      )
      + log_prior(process_model, a_current, w_current)
      + log_prior(measurement_model, gamma_current)
    )

    step_max_abs <- max(abs(do.call(c, current) - do.call(c, previous)))
    previous <- current

    if (trace > 0 && (iteration - 1) %% trace == 0) {
      cat('=====', iteration, '\n')
      cat('a =\n', sprintf('%.06f', current$a), '\n')
      cat('w =\n', sprintf('%.06f', current$w), '\n')
      cat('gamma =\n', sprintf('%.06f', current$gamma), '\n')
      cat('nlm iterations =', output$iterations, '\n')
      cat('log_posterior =', sprintf('%.06f', log_posterior), '\n')
      cat('step_max_abs =', step_max_abs, '\n')
    }

    if (abs(log_posterior - log_posterior_previous) < objective_tolerance) {
      .log_trace('[%d] Terminating because change in log posterior is within tolerance', iteration)
      code <- 1
      break
    }

    if (step_max_abs < step_tolerance) {
      .log_trace('[%d] Terminating because step is within tolerance', iteration)
      code <- 2
      break
    }
  }

  if (iteration == max_iterations) {
    warning('Max iterations reached, convergence may not have been achieved')
  }

  c(mget(c('log_posterior', 'code')), current, .alpha_beta_eta_hat(
    current,
    measurement_model,
    process_model
  ))
}
