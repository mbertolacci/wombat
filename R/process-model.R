#' @export
flux_process_model <- function(
  control_emissions,
  control_mole_fraction,
  perturbations,
  sensitivities,
  lag = months(3),
  H = transport_matrix(
    perturbations,
    control_mole_fraction,
    sensitivities,
    lag
  ),
  alpha,
  alpha_prior_mean = 0,
  a,
  a_prior = c(0.5, 1),
  w,
  w_prior = gamma_quantile_prior(1 / 10 ^ 2, 1 / 0.1 ^ 2),
  Psi = latitudinal_random_effects(control_mole_fraction),
  eta,
  eta_prior_mean = 0,
  eta_prior_variance = 4,
  eta_prior_precision = Diagonal(ncol(Psi), 1 / eta_prior_variance)
) {
  stopifnot_within <- function(x, y) {
    stopifnot(
      all(x %in% y) && all(y %in% x)
    )
  }

  stopifnot_within(control_emissions$model_id, perturbations$model_id)
  stopifnot_within(perturbations$from_month_start, sensitivities$from_month_start)
  stopifnot_within(perturbations$region, sensitivities$region)
  stopifnot_within(sensitivities$model_id, control_mole_fraction$model_id)

  control_mole_fraction <- control_mole_fraction %>%
    arrange(model_id)

  stopifnot(!is.unsorted(control_mole_fraction$time))

  regions <- sort(unique(perturbations$region))
  perturbations <- perturbations %>%
    arrange(from_month_start, region, model_id)

  if (!missing(alpha)) {
    alpha <- .recycle_vector_to(alpha, ncol(H))
  }

  alpha_prior_mean <- .recycle_vector_to(alpha_prior_mean, ncol(H))

  if (!missing(w)) {
    w <- .recycle_vector_to(w, length(regions))
  }

  log_debug('Constructing perturbations basis')
  row_indices <- control_emissions %>%
    mutate(i = 1 : n()) %>%
    select(model_id, i)
  column_indices <- expand.grid(
    region = regions,
    from_month_start = sort(unique(perturbations$from_month_start))
  ) %>%
    mutate(j = 1 : n())
  Phi_df <- perturbations %>%
    left_join(row_indices, by = 'model_id') %>%
    left_join(column_indices, by = c('from_month_start', 'region'))
  Phi <- sparseMatrix(
    i = Phi_df$i,
    j = Phi_df$j,
    x = Phi_df$flux_density,
    dims = c(nrow(control_emissions), nrow(column_indices))
  )

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
    'control_emissions',
    'control_mole_fraction',
    'perturbations',
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
    'control_emissions',
    'control_mole_fraction',
    'perturbations',
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
    H = c('perturbations', 'control_mole_fraction', 'sensitivities', 'lag'),
    Psi = c('control_mole_fraction'),
    eta_prior_precision = c('control_mole_fraction', 'Psi')
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
  perturbations,
  control_mole_fraction,
  sensitivities,
  lag = months(3)
) {
  log_debug('Truncating sensitivities')
  # NOTE(mgnb): do lag computation on the month_start from control because it's
  # much shorter than doing it after it's joined to sensitivities
  control_month_start_lag <- to_month_start(control_mole_fraction$time) - lag
  truncated_sensitivities <- sensitivities %>%
    mutate(
      month_start_lag = control_month_start_lag[
        match(model_id, control_mole_fraction$model_id)
      ]
    ) %>%
    filter(
      month_start_lag <= from_month_start
    )

  log_debug('Adding row and column indices')
  column_indices <- expand.grid(
    region = sort(unique(perturbations$region)),
    from_month_start = sort(unique(perturbations$from_month_start))
  ) %>%
    mutate(column_index = 1 : n())

  row_indices <- data.frame(
    model_id = sort(unique(control_mole_fraction$model_id))
  ) %>%
    mutate(row_index = 1 : n())

  indexed_sensitivities <- truncated_sensitivities %>%
    left_join(column_indices, by = c('region', 'from_month_start')) %>%
    left_join(row_indices, by = c('model_id'))

  log_debug('Constructing transport matrix')
  sparseMatrix(
    i = indexed_sensitivities$row_index,
    j = indexed_sensitivities$column_index,
    x = indexed_sensitivities$co2_sensitivity,
    dims = c(nrow(control_mole_fraction), nrow(column_indices))
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
  log_debug('Contructing basis and design matrix')
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
generate.flux_process_model <- function(model, n_samples = 1) {
  if (is.null(model[['a']])) {
    model$a <- replicate(
      n_samples,
      .rnorm_truncated(1, model$a_prior[1], model$a_prior[2], 0, 1)
    )
  }

  # Not yet supported
  stopifnot(!('lower' %in% model$w_prior))
  w_sample <- t(replicate(n_samples, rgamma(
    length(model$regions),
    shape = model$w_prior[['shape']],
    rate = model$w_prior[['rate']]
  )))
  if (!is.null(model[['w']])) {
    is_fixed <- which(!is.na(model$w))
    fixed_w <- model$w[is_fixed]
    model$w <- w_sample
    for (i in seq_along(is_fixed)) {
      model$w[, is_fixed[i]] <- fixed_w[i]
    }
  } else {
    model$w <- w_sample
  }

  .recycle_index <- function(i, n) {
    1 + ((i - 1) %% n)
  }
  .as_matrix <- function(x) {
    if (length(dim(x)) == 2) x else t(x)
  }

  if (is.null(model[['alpha']])) {
    w_matrix <- .as_matrix(model$w)
    model$alpha <- t(sapply(seq_len(n_samples), function(index) {
      .sample_normal_precision(.make_Q_alpha(model)(list(
        a = model$a[.recycle_index(index, length(model$a))],
        w = w_matrix[.recycle_index(index, nrow(w_matrix)), ]
      )))
    }))
  }

  if (is.null(model[['eta']])) {
    model$eta <- t(replicate(
      n_samples,
      .sample_normal_precision(model$eta_prior_precision)
    ))
  }

  simplify <- function(x) {
    if (is.matrix(x) && nrow(x) == 1) x[1, ] else x
  }

  for (x in c('w', 'alpha', 'eta')) {
    model[[x]] <- simplify(model[[x]])
  }

  model
}

#' @export
calculate.flux_process_model <- function(
  x,
  name = c('Y2_tilde', 'Y2', 'Y2_control', 'Y1_tilde', 'Y1', 'H_alpha'),
  parameters = x
) {
  name <- match.arg(name)

  parameters <- lapply(parameters[c('alpha', 'eta')], function(x) {
    if (is.null(x)) return(x)
    if (is.vector(x)) t(x) else as.matrix(x)
  })

  get_samples <- function(l, r) {
    as.matrix(t(tcrossprod(l, r)))
  }

  add_rowwise <- function(x, y) {
    if (is.vector(x) && is.matrix(y)) {
      # HACK(mgnb): this is the easiest way to add to each row
      t(x + t(y))
    } else {
      x + y
    }
  }

  output <- if (name == 'Y2_tilde') {
    get_samples(x$H, parameters$alpha) + get_samples(x$Psi, parameters$eta)
  } else if (name == 'Y2') {
    add_rowwise(x$control_mole_fraction$co2, calculate(x, 'Y2_tilde', parameters))
  } else if (name == 'Y2_control') {
    x$control_mole_fraction$co2
  } else if (name == 'Y1_tilde') {
    get_samples(x$Phi, parameters$alpha)
  } else if (name == 'Y1') {
    add_rowwise(x$control_emissions$flux_density, calculate(x, 'Y1_tilde', parameters))
  } else if (name == 'H_alpha') {
    get_samples(x$H, parameters$alpha)
  }

  if (is.matrix(output) && nrow(output) == 1) {
    output[1, ]
  } else {
    output
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
aggregate_flux <- function(model, filter_expr = TRUE, parameters = model) {
  if (missing(parameters) && !('alpha' %in% names(model))) {
    parameters <- list(alpha = model$alpha_prior_mean)
  }

  filter_expr <- enquo(filter_expr)

  control_aggregate <- model$control_emissions %>%
    mutate(
      test_condition = !! filter_expr,
      flux = if_else(
        test_condition,
        # kgCO2 / s / m ^ 2 => PgC
        flux_density * area * 10 ^ (3 - 15) / 44.01 * 12.01 * lubridate::days_in_month(month_start) * 24 * 60 * 60,
        0
      )
    ) %>%
    group_by(month_start) %>%
    summarise(
      total_flux = sum(flux)
    )

  matching_control <- model$control_emissions %>%
    filter(!! filter_expr)
  row_indices <- data.frame(
    month_start = sort(unique(model$control_emissions$month_start))
  ) %>%
    mutate(i = 1 : n())
  column_indices <- expand.grid(
    region = model$regions,
    from_month_start = sort(unique(model$perturbations$from_month_start))
  ) %>%
    mutate(j = 1 : n())
  Phi_aggregate_df <- model$perturbations %>%
    filter(model_id %in% matching_control$model_id) %>%
    left_join(
      matching_control %>% select(model_id, month_start, area),
      by = 'model_id'
    ) %>%
    mutate(
      # kgCO2 / s / m ^ 2 => PgC
      flux = flux_density * area * 10 ^ (3 - 15) / 44.01 * 12.01 * lubridate::days_in_month(month_start) * 24 * 60 * 60
    ) %>%
    group_by(from_month_start, region, month_start) %>%
    summarise(
      total_flux = sum(flux)
    ) %>%
    left_join(row_indices, by = 'month_start') %>%
    left_join(column_indices, by = c('from_month_start', 'region'))

  Phi_aggregate <- sparseMatrix(
    i = Phi_aggregate_df$i,
    j = Phi_aggregate_df$j,
    x = Phi_aggregate_df$total_flux,
    dims = c(
      nrow(row_indices),
      nrow(column_indices)
    )
   )

  area_df <- model$control_emissions %>%
    mutate(
      filter_condition = !! filter_expr,
      area = if_else(filter_condition, area, 0)
    ) %>%
    group_by(month_start) %>%
    summarise(
      area = sum(area)
    )

  tibble::tibble(
    month_start = sort(unique(model$control_emissions$month_start)),
    flux = if (is.vector(parameters[['alpha']])) {
      control_aggregate$total_flux + as.vector(Phi_aggregate %*% parameters[['alpha']])
    } else {
      control_aggregate$total_flux + as.matrix(Phi_aggregate %*% t(parameters[['alpha']]))
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
    kronecker(Q_alpha_t, Q_alpha_s)
  }
}
