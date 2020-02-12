#' @export
flux_process_model <- function(
  emissions,
  control,
  sensitivities,
  lag = months(3),
  H = transport_matrix(
    emissions,
    control,
    sensitivities,
    lag
  ),
  alpha,
  alpha_prior_mean = 0,
  a,
  a_prior = c(0.5, 1),
  w,
  w_prior = gamma_quantile_prior(1 / 10 ^ 2, 1 / 0.1 ^ 2),
  Psi = latitudinal_random_effects(control),
  eta,
  eta_prior_mean = 0,
  eta_prior_variance = 4,
  eta_prior_precision = Diagonal(ncol(Psi), 1 / eta_prior_variance)
) {
  control <- control %>%
    arrange(model_id) %>%
    mutate(month_start = to_month_start(time))

  stopifnot(!is.unsorted(control$time))

  regions <- sort(unique(emissions$region))
  emissions <- emissions %>%
    arrange(month_start, region)

  if (!missing(alpha)) {
    alpha <- .recycle_vector_to(alpha, ncol(H))
  }

  alpha_prior_mean <- .recycle_vector_to(alpha_prior_mean, ncol(H))

  if (!missing(w)) {
    w <- .recycle_vector_to(w, length(regions))
  }

  .log_debug('Constructing emissions basis')
  Phi <- emissions %>%
    group_by(month_start, region) %>%
    group_map(~ .x$flux_density) %>%
    bdiag()

  eta_prior_mean <- .recycle_vector_to(eta_prior_mean, ncol(Psi))
  if (!missing(eta)) {
    eta <- .recycle_vector_to(eta, ncol(Psi))
  }

  stopifnot(ncol(H) == ncol(Phi))
  stopifnot(length(a_prior) == 2)
  stopifnot(is.list(w_prior))
  stopifnot(length(w_prior) == 2 || length(w_prior) == 4)
  stopifnot(nrow(eta_prior_precision) == ncol(Psi))
  stopifnot(ncol(eta_prior_precision) == ncol(Psi))

  structure(.remove_nulls_and_missing(mget(c(
    'emissions',
    'control',
    'sensitivities',
    'lag',
    'H',
    'Phi',
    'regions',
    'alpha',
    'alpha_prior_mean',
    'a',
    'a_prior',
    'w',
    'w_prior',
    'Psi',
    'eta',
    'eta_prior_mean',
    'eta_prior_precision'
  ))), class = 'flux_process_model')
}

#' @export
update.flux_process_model <- function(model, ...) {
  current_arguments <- .remove_nulls_and_missing(model[c(
    'emissions',
    'control',
    'sensitivities',
    'lag',
    'H',
    'alpha',
    'alpha_prior_mean',
    'a',
    'a_prior',
    'w',
    'w_prior',
    'Psi',
    'eta',
    'eta_prior_mean',
    'eta_prior_precision'
  )])
  update_arguments <- list(...)

  # Ensure that current_arguments doesn't override anything
  terminal_arguments <- list(
    H = c('emissions', 'control', 'sensitivities', 'lag'),
    Psi = c('control'),
    eta_prior_precision = c('control', 'Psi')
  )
  for (name in names(terminal_arguments)) {
    if (any(terminal_arguments[[name]] %in% names(update_arguments))) {
      current_arguments[[name]] <- NULL
    }
  }

  do.call(flux_process_model, .extend_list(
    current_arguments,
    update_arguments
  ))
}

#' @export
transport_matrix <- function(
  emissions,
  control,
  sensitivities,
  lag = months(3)
) {
  .log_debug('Truncating sensitivities')
  # NOTE(mgnb): do lag computation on the month_start from control because it's
  # much shorter than doing it after it's joined to sensitivities
  control_month_start_lag <- control$month_start - lag
  truncated_sensitivities <- sensitivities %>%
    mutate(
      month_start_lag = control_month_start_lag[
        match(model_id, control$model_id)
      ]
    ) %>%
    filter(
      month_start_lag <= from_month_start
    )

  .log_debug('Adding row and column indices')
  column_indices <- expand.grid(
    region = sort(unique(emissions$region)),
    from_month_start = sort(unique(emissions$month_start))
  ) %>%
    mutate(column_index = 1 : n())

  row_indices <- data.frame(
    model_id = sort(unique(control$model_id))
  ) %>%
    mutate(row_index = 1 : n())

  indexed_sensitivities <- truncated_sensitivities %>%
    left_join(column_indices, by = c('region', 'from_month_start')) %>%
    left_join(row_indices, by = c('model_id'))

  .log_debug('Constructing transport matrix')
  sparseMatrix(
    i = indexed_sensitivities$row_index,
    j = indexed_sensitivities$column_index,
    x = indexed_sensitivities$xco2_sensitivity,
    dims = c(nrow(control), nrow(column_indices))
  )
}

#' @export
latitudinal_random_effects <- function(
  df,
  n_latitude_bands = 5,
  latitude_basis_scale = 50,
  time_varying = TRUE,
  intercept = TRUE
) {
  .log_debug('Contructing basis and design matrix')
  basis <- FRK::local_basis(
    FRK::real_line(),
    loc = matrix(seq(-90, 90, length = n_latitude_bands)),
    scale = rep(latitude_basis_scale, n_latitude_bands)
  )
  df <- df %>% arrange(model_id)
  if (time_varying) {
    output <- df %>%
      select(time, latitude) %>%
      mutate(month_start = to_month_start(time)) %>%
      group_by(month_start) %>%
      group_map(
        ~ FRK::eval_basis(basis, as.matrix(.x$latitude))
      ) %>%
      bdiag()
  } else {
    output <- FRK::eval_basis(basis, as.matrix(df$latitude))
  }

  if (intercept) {
    output <- cbind(1, output)
  }
  attr(output, 'is_latitudinal') <- TRUE
  attr(output, 'n_latitude_bands') <- n_latitude_bands
  attr(output, 'time_varying') <- time_varying
  attr(output, 'intercept') <- intercept
  attr(output, 'basis') <- basis
  output
}

#' @export
generate.flux_process_model <- function(model) {
  if (is.null(model[['a']])) {
    model$a <- .rnorm_truncated(1, model$a_prior[1], model$a_prior[2], 0, 1)
  }

  if (is.null(model[['w']])) {
    # Not yet supported
    stopifnot(!('lower' %in% model$w_prior))
    model$w <- rgamma(
      length(model$regions),
      shape = model$w_prior[['shape']],
      rate = model$w_prior[['rate']]
    )
  }

  if (is.null(model[['alpha']])) {
    model$alpha <- .sample_normal_precision(.make_Q_alpha(model)(model))
  }

  if (is.null(model[['eta']])) {
    model$eta <- .sample_normal_precision(model$eta_prior_precision)
  }

  model
}

#' @export
calculate.flux_process_model <- function(
  x,
  name = c('Y2_tilde', 'Y2', 'Y1_tilde', 'Y1', 'H_alpha'),
  parameters = x
) {
  name <- match.arg(name)

  if (name == 'Y2_tilde') {
    as.vector(x$H %*% parameters$alpha + x$Psi %*% parameters$eta)
  } else if (name == 'Y2') {
    x$control$xco2 + calculate(x, 'Y2_tilde', parameters)
  } else if (name == 'Y1_tilde') {
    as.vector(x$Phi %*% parameters$alpha)
  } else if (name == 'Y1') {
    x$emissions$flux_density + calculate(x, 'Y1_tilde', parameters)
  } else if (name == 'H_alpha') {
    as.vector(x$H %*% parameters$alpha)
  }
}

#' @export
log_prior.flux_process_model <- function(model, parameters = model) {
  output <- 0

  if (is.null(model[['a']])) {
    if (parameters$a <= 0 || parameters$a >= 1) return(-Inf)

    output <- output + dnorm(
      parameters$a,
      mean = model$a_prior[1],
      sd = model$a_prior[2],
      log = TRUE
    )
  }

  if (is.null(model[['w']])) {
    if (any(parameters$w <= 0)) return(-Inf)

    output <- output + sum(dgamma(
      parameters$w,
      shape = model$w_prior[['shape']],
      rate = model$w_prior[['rate']],
      log = TRUE
    ))
  }

  output
}

#' Aggregate monthly fluxes according to a constraint
#' @export
aggregate_flux <- function(model, filter_expr, parameters = model) {
  filter_expr <- enquo(filter_expr)

  Phi_aggregate <- model$emissions %>%
    mutate(
      # kgCO2 / year / m ^ 2 => PgC / year
      flux = if_else(
        !! filter_expr,
        flux_density * area * 10 ^ (6 + 3 - 15) / 44.01 * 12.01,
        0
      )
    ) %>%
    group_by(month_start, region) %>%
    summarise(
      total_flux = sum(flux)
    ) %>%
    group_by(month_start) %>%
    group_map(
      ~ .x$total_flux
    ) %>%
    bdiag() %>%
    t()

  area_df <- model$emissions %>%
    mutate(
      area = if_else(!! filter_expr, area, 0)
    ) %>%
    group_by(month_start) %>%
    summarise(
      area = sum(area)
    )

  tibble::tibble(
    month_start = sort(unique(model$emissions$month_start)),
    flux = if (is.null(dim(parameters$alpha))) {
      as.vector(Phi_aggregate %*% (parameters$alpha + 1))
    } else {
      as.matrix(Phi_aggregate %*% t(parameters$alpha + 1))
    }
  ) %>%
    left_join(area_df, by = 'month_start')
}

.make_Q_alpha <- function(model) {
  n_alpha <- ncol(model$H)
  n_w <- length(model$regions)
  n_times <- n_alpha / n_w

  function(params) {
    Q_alpha_t <- ar1_Q(n_times, params$a)
    Q_alpha_s <- t(sparseMatrix(
      i = seq_len(n_w),
      j = seq_len(n_w),
      x = params$w,
      symmetric = TRUE
    ))
    fast_kronecker(Q_alpha_t, Q_alpha_s)
  }
}
