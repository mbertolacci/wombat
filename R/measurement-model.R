#' @export
flux_measurement_model <- function(
  soundings,
  biases,
  matching,
  process_model,
  C = sparseMatrix(
    i = seq_len(nrow(soundings)),
    j = matching,
    dims = c(nrow(soundings), nrow(process_model$control))
  ),
  attenuation_variables,
  attenuation_factor = if (!missing(attenuation_variables)) {
    interaction(
      as.list(soundings[attenuation_variables]),
      sep = ':'
    )
  } else {
    factor(rep(1, nrow(soundings)))
  },
  measurement_variance = soundings$xco2_error ^ 2,
  measurement_precision = Diagonal(
    nrow(soundings),
    1 / measurement_variance
  ),
  A = model.matrix(biases, soundings)[, -1, drop = FALSE],
  beta,
  beta_prior_mean = 0,
  beta_prior_variance = 4,
  beta_prior_precision = Diagonal(
    ncol(A),
    1 / beta_prior_variance
  ),
  gamma,
  gamma_prior = gamma_quantile_prior(0.1, 1.9)
) {
  if (!missing(beta)) {
    beta <- .recycle_vector_to(beta, ncol(A))
  }
  beta_prior_mean <- .recycle_vector_to(beta_prior_mean, ncol(A))
  if (!missing(gamma)) {
    gamma <- .recycle_vector_to(gamma, nlevels(attenuation_factor))
  }

  stopifnot(nrow(C) == nrow(soundings))
  stopifnot(nrow(attenuation_factor) == nrow(soundings))
  stopifnot(nrow(measurement_precision) == nrow(soundings))
  stopifnot(ncol(measurement_precision) == nrow(soundings))
  stopifnot(nrow(A) == nrow(soundings))
  if (!missing(beta)) {
    stopifnot(length(beta) == ncol(A))
  }
  stopifnot(ncol(beta_prior_precision) == ncol(A))
  stopifnot(nrow(beta_prior_precision) == ncol(A))
  stopifnot(is.list(gamma_prior))
  stopifnot(length(gamma_prior) %in% c(2, 4))

  structure(.remove_nulls_and_missing(mget(c(
    'soundings',
    'C',
    'measurement_precision',
    'biases',
    'A',
    'beta',
    'beta_prior_mean',
    'beta_prior_precision',
    'attenuation_variables',
    'attenuation_factor',
    'gamma',
    'gamma_prior'
  ))), class = 'flux_measurement_model')
}

#' @export
update.flux_measurement_model <- function(model, ...) {
  current_arguments <- .remove_nulls_and_missing(model[c(
    'soundings',
    'C',
    'measurement_precision',
    'biases',
    'A',
    'beta',
    'beta_prior_mean',
    'beta_prior_precision',
    'attenuation_variables',
    'attenuation_factor',
    'gamma',
    'gamma_prior'
  )])
  update_arguments <- list(...)

  # Ensure that current_arguments doesn't override anything
  terminal_arguments <- list(
    A = 'biases',
    C = c('matching', 'process_model'),
    attenuation_factor = 'attenuation_variables',
    measurement_precision = 'measurement_variance',
    beta_prior_precision = c('biases', 'beta_prior_variance')
  )
  for (name in names(terminal_arguments)) {
    if (any(terminal_arguments[[name]] %in% names(update_arguments))) {
      current_arguments[[name]] <- NULL
    }
  }

  do.call(flux_measurement_model, .extend_list(
    current_arguments,
    update_arguments
  ))
}

#' @export
generate.flux_measurement_model <- function(model, process_model) {
  if (is.null(model[['beta']])) {
    model$beta <- (
      model$beta_prior_mean
      + .sample_normal_precision(model$beta_prior_precision)
    )
  }

  if (is.null(model[['gamma']])) {
    model$gamma <- rgamma(
      nlevels(model$attenuation_factor),
      shape = model$gamma_prior[['shape']],
      rate = model$gamma_prior[['rate']]
    )
  }

  if (!missing(process_model)) {
    epsilon <- .sample_normal_precision(.make_Q_epsilon(model)(model))
    model$soundings$xco2 <- as.vector(
      model$C %*% calculate(process_model, 'Y2')
      + model$A %*% model$beta
      + epsilon
    )
  }

  model
}

#' @export
filter.flux_measurement_model <- function(model, expr) {
  expr_s <- substitute(expr)
  indices <- eval(expr_s, model$soundings, parent.frame())

  model$soundings <- model$soundings[indices, , drop = FALSE]
  model$A <- model$A[indices, , drop = FALSE]
  model$C <- model$C[indices, , drop = FALSE]
  model$measurement_precision <- model$measurement_precision[
    indices,
    indices,
    drop = FALSE
  ]
  model$attenuation_factor <- model$attenuation_factor[indices]

  model
}

#' @export
calculate.flux_measurement_model <- function(
  x,
  name = c('Z2_tilde', 'Z2_tilde_debiased'),
  process_model,
  parameters = x
) {
  name <- match.arg(name)

  if (name == 'Z2_tilde') {
    x$soundings$xco2 - as.vector(x$C %*% process_model$control$xco2)
  } else if (name == 'Z2_tilde_debiased') {
    calculate(x, 'Z2_tilde', process_model) - as.vector(x$A %*% parameters$beta)
  }
}

#' @export
log_prior.flux_measurement_model <- function(model, parameters = model) {
  if (is.null(model[['gamma']])) {
    gamma <- parameters$gamma

    if ('lower' %in% names(model$gamma_prior)) {
      if (
        any(gamma <= model$gamma_prior[['lower']])
        || any(gamma >= model$gamma_prior[['upper']])
      ) {
        return(-Inf)
      }
    } else {
      if (any(gamma <= 0)) {
        return(-Inf)
      }
    }

    sum(dgamma(
      gamma,
      shape = model$gamma_prior[['shape']],
      rate = model$gamma_prior[['rate']],
      log = TRUE
    ))
  } else {
    0
  }
}

.make_Xt_Q_epsilon_X <- function(X, model) {
  attenuation_index <- as.integer(model$attenuation_factor)
  n_gamma <- nlevels(model$attenuation_factor)

  Xt_Q_epsilon0_X_parts <- vector('list', n_gamma)
  for (gamma_index in 1 : n_gamma) {
    indices <- attenuation_index == gamma_index
    X_part <- X[indices, ]
    Q_epsilon0_part <- model$measurement_precision[indices, indices]
    Xt_Q_epsilon0_X_parts[[gamma_index]] <- crossprod(
      chol(Q_epsilon0_part) %*% X_part
    )
  }

  function(params = list(gamma = rep(1, n_gamma)), parts = 1 : n_gamma) {
    Reduce(function(l, i) {
      if (is.null(l)) {
        params$gamma[i] * Xt_Q_epsilon0_X_parts[[i]]
      } else {
        .fast_add(l, params$gamma[i] * Xt_Q_epsilon0_X_parts[[i]])
      }
    }, parts, NULL)
  }
}

.make_A_Q_epsilon_B <- function(A, B, model) {
  attenuation_index <- as.integer(model$attenuation_factor)
  n_gamma <- nlevels(model$attenuation_factor)

  A_Q_epsilon0_B_parts <- vector('list', n_gamma)
  for (gamma_index in 1 : n_gamma) {
    indices <- attenuation_index == gamma_index
    A_part <- A[, indices]
    B_part <- B[indices, ]
    Q_epsilon0_part <- model$measurement_precision[indices, indices]
    A_Q_epsilon0_B_parts[[gamma_index]] <- (
      A_part %*% Q_epsilon0_part %*% B_part
    )
  }

  function(params = list(gamma = rep(1, n_gamma)), parts = 1 : n_gamma) {
    Reduce(function(l, i) {
      if (is.null(l)) {
        params$gamma[i] * A_Q_epsilon0_B_parts[[i]]
      } else {
        .fast_add(l, params$gamma[i] * A_Q_epsilon0_B_parts[[i]])
      }
    }, parts, NULL)
  }
}

.make_Q_epsilon <- function(model) {
  n_gamma <- nlevels(model$attenuation_factor)
  function(params = list(gamma = rep(1, n_gamma))) {
    output <- model$measurement_precision
    output@x <- output@x * params$gamma[model$attenuation_factor]
    output
  }
}
