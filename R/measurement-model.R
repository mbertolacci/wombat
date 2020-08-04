#' @export
flux_measurement_model <- function(
  observations,
  biases,
  matching,
  process_model,
  C = sparseMatrix(
    i = seq_len(nrow(observations)),
    j = matching,
    dims = c(nrow(observations), nrow(process_model$control_mole_fraction))
  ),
  attenuation_variables,
  attenuation_factor = if (!missing(attenuation_variables)) {
    interaction(
      as.list(observations[attenuation_variables]),
      sep = ':'
    )
  } else {
    factor(rep(1, nrow(observations)))
  },
  measurement_variance = observations$co2_error ^ 2,
  measurement_covariance = Diagonal(
    nrow(observations),
    measurement_variance
  ),
  A = model.matrix(biases, observations)[, -1, drop = FALSE],
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
  stopifnot(is.factor(attenuation_factor))
  stopifnot(nrow(C) == nrow(observations))
  stopifnot(length(attenuation_factor) == nrow(observations))
  stopifnot(nrow(measurement_covariance) == nrow(observations))
  stopifnot(ncol(measurement_covariance) == nrow(observations))
  stopifnot(nrow(A) == nrow(observations))
  if (!missing(beta)) {
    stopifnot(length(beta) == ncol(A))
  }
  stopifnot(ncol(beta_prior_precision) == ncol(A))
  stopifnot(nrow(beta_prior_precision) == ncol(A))
  stopifnot(is.list(gamma_prior))
  stopifnot(length(gamma_prior) %in% c(2, 4))

  structure(.remove_nulls_and_missing(mget(c(
    'observations',
    'C',
    'measurement_covariance',
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
    'observations',
    'C',
    'measurement_covariance',
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
    measurement_covariance = 'measurement_variance',
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
    epsilon <- .sample_normal_covariance(model$measurement_covariance) / sqrt(
      model$gamma[model$attenuation_factor]
    )
    model$observations$co2 <- as.vector(
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
  indices <- eval(expr_s, model$observations, parent.frame())

  model$observations <- model$observations[indices, , drop = FALSE]

  model$A <- model$A[indices, , drop = FALSE]
  non_empty_columns <- colSums(abs(model$A)) != 0
  model$A <- model$A[, non_empty_columns, drop = FALSE]
  model$beta_prior_mean <- model$beta_prior_mean[non_empty_columns]
  model$beta_prior_precision <- model$beta_prior_precision[
    non_empty_columns,
    non_empty_columns,
    drop = FALSE
  ]
  if (!is.null(model[['beta']])) {
    model[['beta']] <- model[['beta']][non_empty_columns]
  }

  model$C <- model$C[indices, , drop = FALSE]
  model$measurement_covariance <- model$measurement_covariance[
    indices,
    indices,
    drop = FALSE
  ]
  original_levels <- levels(model$attenuation_factor)
  model$attenuation_factor <- droplevels(model$attenuation_factor[indices])
  if (!is.null(model[['gamma']])) {
    model[['gamma']] <- model[['gamma']][
      original_levels %in% levels(model$attenuation_factor)
    ]
  }

  model
}

#' @export
calculate.flux_measurement_model <- function(
  x,
  name = c(
    'Z2',
    'Z2_hat',
    'Z2_debiased',
    'Z2_tilde',
    'Z2_tilde_debiased',
    'Z2_tilde_hat',
    'Y2',
    'Y2_control',
    'Y2_tilde'
  ),
  process_model,
  parameters = x
) {
  name <- match.arg(name)

  parameters <- lapply(parameters[c('beta', 'eta', 'alpha')], function(x) {
    if (is.null(x)) return(x)
    if (is.vector(x)) {
      t(x)
    } else {
      if (ncol(x) == 0) {
        x
      } else {
        as.matrix(x)
      }
    }
  })

  calculate_bias_correction <- function() {
    if (ncol(x$A) == 0) 0
    else as.matrix(t(tcrossprod(x$A, parameters$beta)))
  }

  add_rowwise <- function(x, y) {
    if (is.vector(y) && is.matrix(x)) {
      tmp <- x
      x <- y
      y <- tmp
    }
    if (is.vector(x) && is.matrix(y)) {
      # HACK(mgnb): this is the easiest way to add to each row
      t(x + t(y))
    } else {
      x + y
    }
  }

  output <- if (name == 'Z2') {
    x$observations$co2
  } else if (name == 'Z2_hat') {
    add_rowwise(
      calculate(x, 'Y2', process_model, parameters),
      calculate_bias_correction()
    )
  } else if (name == 'Z2_debiased') {
    add_rowwise(
      calculate(x, 'Z2', process_model, parameters),
      -calculate_bias_correction()
    )
  } else if (name == 'Z2_tilde') {
    add_rowwise(
      calculate(x, 'Z2', process_model, parameters),
      -calculate(x, 'Y2_control', process_model, parameters)
    )
  } else if (name == 'Z2_tilde_debiased') {
    add_rowwise(
      calculate(x, 'Z2_tilde', process_model, parameters),
      -calculate_bias_correction()
    )
  } else if (name == 'Z2_tilde_hat') {
    add_rowwise(
      calculate(x, 'Z2_hat', process_model, parameters),
      -calculate(x, 'Y2_control', process_model, parameters)
    )
  } else {
    # Anything left falls through to the process model
    rhs <- calculate(process_model, name, parameters)
    if (is.matrix(rhs)) {
      as.matrix(tcrossprod(rhs, x$C))
    } else {
      as.matrix(t(x$C %*% rhs))
    }
  }

  if (is.matrix(output) && nrow(output) == 1) {
    output[1, ]
  } else {
    output
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
  n_per_gamma <- table(model$attenuation_factor)

  zero_rows <- function(A, include_rows) {
    At <- as(A, 'dgTMatrix')
    include <- (At@i + 1) %in% include_rows
    At@i <- At@i[include]
    At@j <- At@j[include]
    At@x <- At@x[include]
    as(At, 'dgCMatrix')
  }

  part_ij <- rbind(
    cbind(seq_len(n_gamma), seq_len(n_gamma)),
    if (n_gamma > 1) t(combn(n_gamma, 2)) else NULL
  )
  log_trace('Forming parts of Xt_Q_epsilon_X')
  part_x <- apply(part_ij, 1, function(ij) {
    i <- ij[1]
    j <- ij[2]

    if (n_per_gamma[i] == 0 || n_per_gamma[j] == 0) {
      return(sparseMatrix(i = NULL, j = NULL, dims = c(ncol(X), ncol(X))))
    }

    Xi <- zero_rows(X, which(attenuation_index == i))
    if (i == j) {
      forceSymmetric(crossprod(Xi, solve(model$measurement_covariance, Xi)))
    } else {
      Xj <- zero_rows(X, which(attenuation_index == j))
      forceSymmetric(
        crossprod(Xi, solve(model$measurement_covariance, Xj))
        + crossprod(Xj, solve(model$measurement_covariance, Xi))
      )
    }
  })

  # NOTE(mgnb): drop parts equal to zero
  keep_part <- sapply(part_x, function(x) {
    nnzero(x) > 0
  })
  part_ij <- part_ij[keep_part, , drop = FALSE]
  part_x <- part_x[keep_part]

  function(params = list(gamma = rep(1, n_gamma))) {
    gamma_ij <- (
      sqrt(params$gamma)[part_ij[, 1]] * sqrt(params$gamma)[part_ij[, 2]]
    )

    Reduce(function(l, k) {
      if (is.null(l)) {
        gamma_ij[k] * part_x[[k]]
      } else {
        .fast_add(l, gamma_ij[k] * part_x[[k]])
      }
    }, seq_len(nrow(part_ij)), NULL)
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

.make_solve_Q_epsilon <- function(model) {
  function(x, params = list(gamma = rep(1, n_gamma))) {
    D <- Diagonal(x = sqrt(params$gamma[model$attenuation_factor]))
    D %*% solve(model$measurement_covariance, D %*% x)
  }
}

# .make_Q_epsilon <- function(model) {
#   n_gamma <- nlevels(model$attenuation_factor)
#   function(params = list(gamma = rep(1, n_gamma))) {
#     output <- model$measurement_precision
#     output@x <- output@x * params$gamma[model$attenuation_factor]
#     output
#   }
# }
