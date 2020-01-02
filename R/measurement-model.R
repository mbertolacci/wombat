#' @export
flux_measurement_model <- function(
  soundings,
  biases,
  matching,
  process_model,
  C = sparseMatrix(
    i = 1 : nrow(soundings),
    j = matching,
    dims = c(nrow(soundings), nrow(process_model$control))
  ),
  attenuation_variables = NULL,
  measurement_variance = soundings$xco2_error ^ 2,
  measurement_precision = Diagonal(
    x = 1 / measurement_variance
  ),
  A = model.matrix(update(biases, ~ . - 1), soundings),
  beta = NULL,
  beta_prior_mean = rep(0, ncol(A)),
  beta_prior_variance = 4,
  beta_prior_precision = Diagonal(ncol(A), 1 / beta_prior_variance),
  gamma = NULL,
  gamma_prior = gamma_quantile_prior(0.1, 1.9)
) {
  soundings <- soundings %>%
    arrange(sounding_id)

  if (!is.null(beta)) {
    beta_prior_mean <- NULL
    beta_prior_precision <- NULL
  }

  if (!is.null(attenuation_variables)) {
    attenuation_factor <- interaction(
      as.list(soundings[attenuation_variables]),
      sep = ':'
    )
  } else {
    attenuation_factor <- factor(rep(1, nrow(soundings)))
    gamma <- 1
    gamma_prior <- NULL
  }

  structure(.remove_nulls(mget(c(
    'soundings',
    'C',
    'measurement_precision',
    'A',
    'beta',
    'beta_prior_mean',
    'beta_prior_precision',
    'attenuation_factor',
    'gamma',
    'gamma_prior'
  ))), class = 'flux_measurement_model')
}

#' @export
generate.flux_measurement_model <- function(model, process_sample) {
  if (is.null(model[['beta']])) {
    model$beta <- (
      model$beta_prior_mean
      + .sample_normal_precision(model$beta_prior_precision)
    )
  }

  if (is.null(model[['gamma']])) {
    model$gamma <- rgamma(
      nlevels(model$attenuation_factor),
      shape = model$gamma_prior[1],
      rate = model$gamma_prior[2]
    )
  }

  D <- Diagonal(x = model$gamma[as.integer(model$attenuation_factor)])
  epsilon_precision <- D %*% model$measurement_precision
  epsilon <- .sample_normal_precision(epsilon_precision)

  output <- c(
    model[c('beta', 'gamma')],
    mget('epsilon')
  )

  if (!missing(process_sample)) {
    output$Z2_tilde <- as.vector(
      model$C %*% process_sample$Y2_tilde
      + model$A %*% model$beta
      + epsilon
    )
  }

  output
}

#' @export
log_prior.flux_measurement_model <- function(model, gamma) {
  if (is.null(model$gamma)) {
    sum(dgamma(
      gamma,
      shape = model$gamma_prior[1],
      rate = model$gamma_prior[2],
      log = TRUE
    ))
  } else {
    0
  }
}

# .make_Xt_Q_epsilon_X <- function(X, model) {
#   attenuation_index <- as.integer(model$attenuation_factor)
#   n_gamma <- nlevels(model$attenuation_factor)

#   Xt_Q_epsilon0_X_parts <- array(NA, dim = c(n_gamma, ncol(X), ncol(X)))
#   for (gamma_index in 1 : n_gamma) {
#     indices <- attenuation_index == gamma_index
#     X_part <- X[indices, ]
#     Q_epsilon0_part <- model$measurement_precision[indices, indices]
#     Xt_Q_epsilon0_X_parts[gamma_index, , ] <- as.matrix(crossprod(
#       sqrt(Q_epsilon0_part) %*% X_part
#     ))
#   }

#   function(params = list(gamma = rep(1, n_gamma)), parts = 1 : n_gamma) {
#     colSums(params$gamma[parts] * Xt_Q_epsilon0_X_parts[parts, , , drop = FALSE])
#   }
# }


.make_Xt_Q_epsilon_X <- function(X, model) {
  attenuation_index <- as.integer(model$attenuation_factor)
  n_gamma <- nlevels(model$attenuation_factor)

  Xt_Q_epsilon0_X_parts <- vector('list', n_gamma)
  for (gamma_index in 1 : n_gamma) {
    indices <- attenuation_index == gamma_index
    X_part <- X[indices, ]
    Q_epsilon0_part <- model$measurement_precision[indices, indices]
    Xt_Q_epsilon0_X_parts[[gamma_index]] <- crossprod(
      sqrt(Q_epsilon0_part) %*% X_part
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

# .make_A_Q_epsilon_B <- function(A, B, model) {
#   attenuation_index <- as.integer(model$attenuation_factor)
#   n_gamma <- nlevels(model$attenuation_factor)

#   A_Q_epsilon0_B_parts <- array(NA, dim = c(n_gamma, nrow(A), ncol(B)))
#   for (gamma_index in 1 : n_gamma) {
#     indices <- attenuation_index == gamma_index
#     A_part <- A[, indices]
#     B_part <- B[indices, ]
#     Q_epsilon0_part <- model$measurement_precision[indices, indices]
#     A_Q_epsilon0_B_parts[gamma_index, , ] <- as.matrix(
#       A_part %*% Q_epsilon0_part %*% B_part
#     )
#   }

#   function(params = list(gamma = rep(1, n_gamma)), parts = 1 : n_gamma) {
#     colSums(params$gamma[parts] * A_Q_epsilon0_B_parts[parts, , , drop = FALSE])
#   }
# }

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
