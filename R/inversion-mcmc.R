#' @export
inversion_mcmc <- function(
  n_iterations = 1000,
  measurement_model,
  process_model,
  start = NULL,
  warm_up = 100,
  tuning = list(
    a = list(w = 0.25, max_evaluations = 100),
    gamma = list(w = 1, max_evaluations = 100)
  ),
  marginalised_variables = NULL,
  show_progress = TRUE
) {
  if (!missing(tuning)) {
    tuning <- .extend_list(eval(formals(inversion_mcmc)$tuning), tuning)
  }

  if (!is.null(marginalised_variables)) {
    marginalised_variables <- match.arg(
      marginalised_variables,
      c('eta', 'beta'),
      several.ok = TRUE
    )
  }
  included_variables <- setdiff(
    c('alpha', 'eta', 'beta'),
    marginalised_variables
  )

  n_alpha <- ncol(process_model$H)
  n_eta <- ncol(process_model$Psi)
  n_beta <- ncol(measurement_model$A)
  n_times <- n_alpha / length(process_model$regions)
  n_w <- length(process_model$regions)
  n_gamma <- nlevels(measurement_model$attenuation_factor)

  .log_debug('Precomputing various quantities')
  omega_sampler <- .make_omega_sampler(
    measurement_model,
    process_model,
    variables = included_variables
  )
  a_sampler <- .make_a_sampler(process_model, tuning[['a']])
  w_sampler <- .make_w_sampler(process_model)
  gamma_sampler <- .make_gamma_sampler(
    measurement_model,
    process_model,
    variables = included_variables,
    tuning = tuning[['gamma']]
  )

  .log_debug('Setting start values')
  start <- .extend_list(
    c(
      generate(process_model)[c('alpha', 'eta', 'a', 'w')],
      generate(measurement_model)[c('beta', 'gamma')]
    ),
    start,
    .remove_nulls(process_model[c('alpha', 'eta', 'a', 'w')]),
    .remove_nulls(process_model[c('beta', 'gamma')])
  )[c('alpha', 'eta', 'beta', 'a', 'w', 'gamma')]
  for (variable in marginalised_variables) {
    start[[variable]] <- NA
  }

  alpha_samples <- matrix(NA, nrow = n_iterations, ncol = n_alpha)
  eta_samples <- matrix(NA, nrow = n_iterations, ncol = n_eta)
  a_samples <- matrix(NA, nrow = n_iterations, ncol = 1)
  w_samples <- matrix(NA, nrow = n_iterations, ncol = n_w)
  beta_samples <- matrix(NA, nrow = n_iterations, ncol = n_beta)
  gamma_samples <- matrix(NA, nrow = n_iterations, ncol = n_gamma)

  current <- start
  alpha_samples[1, ] <- current$alpha
  eta_samples[1, ] <- current$eta
  a_samples[1, ] <- current$a
  w_samples[1, ] <- current$w
  beta_samples[1, ] <- current$beta
  gamma_samples[1, ] <- current$gamma

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

    .log_trace('[%d/%d] Sampling omega', iteration, n_iterations)
    current <- omega_sampler(current)

    .log_trace('[%d/%d] Sampling a', iteration, n_iterations)
    current <- a_sampler(current, iteration <= warm_up)

    .log_trace('[%d/%d] Sampling w', iteration, n_iterations)
    current <- w_sampler(current)

    .log_trace('[%d/%d] Sampling gamma', iteration, n_iterations)
    current <- gamma_sampler(current, iteration <= warm_up)

    alpha_samples[iteration, ] <- current$alpha
    eta_samples[iteration, ] <- current$eta
    a_samples[iteration, ] <- current$a
    w_samples[iteration, ] <- current$w
    beta_samples[iteration, ] <- current$beta
    gamma_samples[iteration, ] <- current$gamma
  }

  region_month <- expand.grid(
    region = process_model$regions,
    month_index = seq_len(length(unique(process_model$emissions$month_start)))
  )
  colnames(alpha_samples) <- sprintf(
    'alpha[%d, %d]',
    region_month$region,
    region_month$month_index
  )
  colnames(eta_samples) <- sprintf('eta[%d]', seq_len(ncol(eta_samples)))
  colnames(a_samples) <- 'a'
  colnames(w_samples) <- sprintf('w[%d]', seq_len(ncol(w_samples)))
  colnames(beta_samples) <- sprintf('beta[%d]', seq_len(ncol(beta_samples)))
  colnames(gamma_samples) <- sprintf('gamma[%d]', seq_len(ncol(gamma_samples)))

  structure(
    list(
      alpha = coda::mcmc(alpha_samples),
      eta = coda::mcmc(eta_samples),
      a = coda::mcmc(a_samples),
      w = coda::mcmc(w_samples),
      beta = coda::mcmc(beta_samples),
      gamma = coda::mcmc(gamma_samples)
    ),
    class = 'flux_inversion_mcmc'
  )
}

#' @export
window.flux_inversion_mcmc <- function(object, ...) {
  for (name in names(object)) {
    object[[name]] <- window(object[[name]], ...)
  }
  object
}

#' @export
plot_traces <- function(object, n_columns = 4) {
  matrix_to_long_df <- function(x) {
    as.data.frame(x) %>%
      mutate(iteration = 1 : n()) %>%
      tidyr::pivot_longer(-iteration) %>%
      mutate(name = factor(name, levels = colnames(x)))
  }

  trace_plot <- function(x) {
    df <- matrix_to_long_df(x)

    df_summary <- df %>%
      group_by(name) %>%
      summarise(
        q10 = quantile(value, probs = 0.10),
        q50 = quantile(value, probs = 0.50),
        mean = mean(value),
        q90 = quantile(value, probs = 0.90)
      ) %>%
      ungroup() %>%
      mutate(
        label = sprintf(
          'mean = %.03g | 10%% = %.03g | 50%% = %.03g | 90%% = %.03g',
          mean,
          q10,
          q50,
          q90
        )
      )

    ggplot2::ggplot(
      df,
      ggplot2::aes(iteration, value)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_text(
        data = df_summary,
        mapping = aes(x = -Inf, y = Inf, label = label),
        hjust = 'left',
        vjust = 'top'
      ) +
      ggplot2::facet_wrap(
        ~ name,
        scales = 'free_y',
        ncol = n_columns,
        labeller = 'label_parsed'
      ) +
      ggplot2::labs(x = NULL, y = NULL)
  }

  alpha_subset <- object$alpha[
    ,
    seq(1, ncol(object$alpha), by = ceiling(ncol(object$alpha) / 8))
  ]
  eta_subset <- object$eta[
    ,
    seq(1, ncol(object$eta), by = ceiling(ncol(object$eta) / 8))
  ]

  gridExtra::grid.arrange(
    trace_plot(alpha_subset),
    trace_plot(eta_subset),
    trace_plot(object$beta),
    trace_plot(object$a),
    trace_plot(object$w),
    trace_plot(object$gamma),
    left = 'Value',
    bottom = 'Iteration',
    heights = ceiling(c(
      ncol(alpha_subset),
      ncol(eta_subset),
      ncol(object$beta),
      n_columns,
      ncol(object$w),
      ncol(object$gamma)
    ) / n_columns)
  )
}
