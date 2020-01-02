#' @export
inversion_mcmc <- function(
  n_iterations = 1000,
  measurement_model,
  process_model,
  start = NULL,
  warm_up = 100,
  tuning = list(
    slice_a_options = list(w = 0.25, max_evaluations = 100),
    slice_gamma_options = list(w = 1, max_evaluations = 100)
  ),
  show_progress = TRUE
) {
  if (!missing(tuning)) {
    tuning <- .extend_list(eval(formals(inversion_mcmc)$tuning), tuning)
  }

  n_alpha <- ncol(process_model$H)
  n_beta <- ncol(measurement_model$A)
  n_times <- n_alpha / length(process_model$regions)
  n_w <- length(process_model$regions)
  n_gamma <- nlevels(measurement_model$attenuation_factor)

  find_a <- is.null(process_model[['a']])
  find_w <- is.null(process_model[['w']])
  find_gamma <- is.null(measurement_model[['gamma']])

  .log_debug('Precomputing various quantities')
  Z <- measurement_model$soundings$xco2 - as.vector(measurement_model$C %*% process_model$control$xco2)
  C <- measurement_model$C
  X <- cbind(
    C %*% process_model$H,
    measurement_model$A
  )
  Psi_C <- C %*% process_model$Psi
  Xt_Q_epsilon_X <- .make_Xt_Q_epsilon_X(X, measurement_model)
  Psi_Ct_Q_epsilon_Psi_C <- .make_Xt_Q_epsilon_X(Psi_C, measurement_model)
  Psi_Ct_Q_epsilon_X <- .make_A_Q_epsilon_B(t(Psi_C), X, measurement_model)
  Q_epsilon <- .make_Q_epsilon(measurement_model)
  Q_alpha <- .make_Q_alpha(process_model)
  Q_omega <- .make_Q_omega(
    process_model,
    measurement_model,
    c('alpha', 'beta')
  )
  chol_Q_eta_conditional <- function(params) {
    chol(.fast_add(
      process_model$eta_prior_precision,
      Psi_Ct_Q_epsilon_Psi_C(params)
    ))
  }

  .log_debug('Setting start values')
  start <- .extend_list(
    c(
      generate(process_model)[c('alpha', 'a', 'w')],
      generate(measurement_model)[c('beta', 'gamma')]
    ),
    start,
    .remove_nulls(process_model[c('alpha', 'a', 'w')]),
    .remove_nulls(process_model[c('beta', 'gamma')])
  )[c('alpha', 'beta', 'a', 'w', 'gamma')]

  alpha_samples <- matrix(NA, nrow = n_iterations, ncol = n_alpha)
  beta_samples <- matrix(NA, nrow = n_iterations, ncol = n_beta)
  a_samples <- rep(NA, n_iterations)
  w_samples <- matrix(NA, nrow = n_iterations, ncol = n_w)
  gamma_samples <- matrix(NA, nrow = n_iterations, ncol = n_gamma)

  current <- start
  alpha_samples[1, ] <- current$alpha
  beta_samples[1, ] <- current$beta
  a_samples[1] <- current$a
  w_samples[1, ] <- current$w
  gamma_samples[1, ] <- current$gamma

  a_slice <- do.call(slice, tuning$slice_a_options)
  gamma_slice <- lapply(1 : n_gamma, function(i) {
    do.call(slice, tuning$slice_gamma_options)
  })

  if (show_progress) {
    pb <- progress::progress_bar$new(
      total = n_iterations,
      format = '[:bar] :current/:total eta: :eta'
    )
    pb$tick()
  }

  .log_debug('Starting sampler')
  for (iteration in 2 : n_iterations) {
    if (show_progress) {
      pb$tick()
    }

    .log_trace('[%d/%d] Sampling alpha and beta', iteration, n_iterations)
    chol_Q_eta_conditional_i <- chol_Q_eta_conditional(current)
    chol_Q_omega_conditional <- chol(
      Q_omega(current)
      + Xt_Q_epsilon_X(current)
      - crossprod(solve(t(chol_Q_eta_conditional_i), Psi_Ct_Q_epsilon_X(current)))
    )
    mu_omega_conditional <- as.vector(.chol_solve(
      chol_Q_omega_conditional,
      crossprod(
        X,
        Q_epsilon(current) %*% (
          Z - Psi_C %*% .chol_solve(
            chol_Q_eta_conditional_i,
            crossprod(Psi_C, Q_epsilon(current) %*% Z)
          )
        )
      )
    ))
    omega <- mu_omega_conditional + .sample_normal_precision_chol(chol_Q_omega_conditional)
    current$alpha <- head(omega, n_alpha)
    current$beta <- tail(omega, n_beta)

    if (find_a) {
      .log_trace('[%d/%d] Sampling a', iteration, n_iterations)
      output <- a_slice(current$a, function(a) {
        if (a <= 0 || a >= 1) return(-Inf)
        current2 <- current
        current2$a <- a

        output <- (
          .log_pdf_precision_cholesky(current$alpha, chol(Q_alpha(current2)))
          + log_prior(process_model, a, current$w)
        )
        output
      }, learn = (iteration <= warm_up), include_n_evaluations = TRUE)
      .log_trace('[%d/%d] a took %d evaluations, w = %f', iteration, n_iterations, output$n_evaluations, output$w)

      current$a <- output$sample
    }

    if (find_w) {
      .log_trace('[%d/%d] Sampling w', iteration, n_iterations)
      chol_Q_alpha_t <- chol(ar1_Q(n_times, current$a))
      alpha_matrix <- matrix(current$alpha, nrow = n_times, byrow = TRUE)
      current$w <- rgamma(
        n_w,
        shape = process_model$w_prior[1] + 0.5 * n_times,
        rate =  process_model$w_prior[2] + 0.5 * colSums((chol_Q_alpha_t %*% alpha_matrix) ^ 2)
      )
    }

    if (find_gamma) {
      .log_trace('[%d/%d] Sampling gamma', iteration, n_iterations)
      Z_tilde <- as.vector(Z - X %*% c(current$alpha, current$beta))
      current$gamma <- sapply(1 : n_gamma, function(i) {
        output <- gamma_slice[[i]](current$gamma[i], function(gamma_i) {
          if (gamma_i <= 0) return(-Inf)
          current2 <- current
          current2$gamma[i] <- gamma_i

          (
            .log_pdf_woodbury(
              Z_tilde,
              Q_epsilon(current2),
              chol_Q_eta_conditional(current2),
              process_model$eta_prior_precision,
              Psi_C
            )
            + log_prior(
              measurement_model,
              current2$gamma
            )
          )
        }, learn = (iteration <= warm_up), include_n_evaluations = TRUE)

        .log_trace('[%d/%d] gamma[%d] took %d evaluations, w = %f', iteration, n_iterations, i, output$n_evaluations, output$w)

        output$sample
      })
    }

    alpha_samples[iteration, ] <- current$alpha
    beta_samples[iteration, ] <- current$beta
    a_samples[iteration] <- current$a
    w_samples[iteration, ] <- current$w
    gamma_samples[iteration, ] <- current$gamma
  }

  list(
    alpha = alpha_samples,
    beta = beta_samples,
    a = a_samples,
    w = w_samples,
    gamma = gamma_samples
  )
}
