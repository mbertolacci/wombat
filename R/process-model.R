#' @export
flux_process_model <- function(
  emissions,
  control,
  sensitivities,
  lag = months(3),
  a = NULL,
  a_prior = c(0.5, 1),
  w = NULL,
  w_prior = gamma_quantile_prior(1 / 0.5 ^ 2, 1 / 0.1 ^ 2),
  Psi = latitudinal_random_effects(control),
  eta = NULL,
  eta_prior_mean = rep(0, ncol(Psi)),
  eta_prior_variance = 1,
  eta_prior_precision = Diagonal(ncol(Psi), 1 / eta_prior_variance)
) {
  stopifnot(ncol(eta_prior_precision) == nrow(eta_prior_precision))

  control <- control %>%
    arrange(model_id) %>%
    mutate(month_start = to_month_start(time))

  stopifnot(!is.unsorted(control$time))

  regions <- sort(unique(emissions$region))
  emissions <- emissions %>%
    arrange(month_start, region)

  .log_debug('Constructing emissions basis')
  Phi <- emissions %>%
    group_by(month_start, region) %>%
    group_map(~ .x$flux) %>%
    bdiag()

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
    region = regions,
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
  H <- sparseMatrix(
    i = indexed_sensitivities$row_index,
    j = indexed_sensitivities$column_index,
    x = indexed_sensitivities$xco2_sensitivity,
    dims = c(nrow(control), nrow(column_indices))
  )

  stopifnot(ncol(H) == ncol(Phi))
  stopifnot(ncol(Psi) == ncol(eta_prior_precision))

  structure(.remove_nulls(mget(c(
    'emissions',
    'Phi',
    'H',
    'control',
    'sensitivities',
    'regions',
    'a',
    'a_prior',
    'w',
    'w_prior',
    'Psi',
    'eta',
    'eta_prior_precision'
  ))), class = 'flux_process_model')
}

#' @export
latitudinal_random_effects <- function(
  df,
  n_latitude_bands = 5,
  latitude_basis_scale = 50
) {
  .log_debug('Contructing basis and design matrix')
  basis <- FRK::local_basis(
    FRK::real_line(),
    loc = matrix(seq(-90, 90, length = n_latitude_bands)),
    scale = rep(latitude_basis_scale, n_latitude_bands)
  )
  Psi <- df %>%
    select(time, latitude) %>%
    mutate(month_start = to_month_start(time)) %>%
    group_by(month_start) %>%
    group_map(
      ~ FRK::eval_basis(basis, as.matrix(.x$latitude))
    ) %>%
    bdiag()

  output <- cbind(1, Psi)
  attr(output, 'basis') <- basis
  output
}

#' @export
generate.flux_process_model <- function(model) {
  if (is.null(model[['a']])) {
    model$a <- .rnorm_truncated(1, model$a_prior[1], model$a_prior[2], 0, 1)
  }

  if (is.null(model[['w']])) {
    model$w <- rgamma(
      length(model$regions),
      shape = model$w_prior[1],
      rate = model$w_prior[2]
    )
  }

  n_alpha <- ncol(model$H)
  n_times <- n_alpha / length(model$regions)

  # Sample alpha
  Q_alpha_t <- ar1_Q(n_times, model$a)
  Q_alpha_s <- Diagonal(length(model$regions), x = model$w)
  Q_alpha <- kronecker(Q_alpha_t, Q_alpha_s)
  alpha <- .sample_normal_precision(Q_alpha)

  # Sample eta
  if (is.null(model[['eta']])) {
    model$eta <- .sample_normal_precision(model$eta_prior_precision)
  }

  # Calculate values
  Y1_tilde <- as.vector(model$Phi %*% alpha)
  Y2_tilde <- as.vector(
    model$H %*% alpha
    + model$Psi %*% model$eta
  )

  c(
    model[c('a', 'eta', 'w')],
    mget(c('alpha', 'Q_alpha', 'Y1_tilde', 'Y2_tilde'))
  )
}

#' @export
log_prior.flux_process_model <- function(model, a, w) {
  output <- 0

  if (is.null(model$w)) {
    output <- output + dnorm(
      a,
      mean = model$a_prior[1],
      sd = model$a_prior[2],
      log = TRUE
    )
  }

  if (is.null(model$w)) {
    output <- output + sum(dgamma(
      w,
      shape = model$w_prior[1],
      rate = model$w_prior[2],
      log = TRUE
    ))
  }

  output
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
