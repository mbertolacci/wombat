.make_gamma_sampler <- function(
  measurement_model,
  process_model,
  tuning = list(w = 1, max_evaluations = 100),
  X = .make_X_omega(process_model, measurement_model),
  Q_epsilon = .make_Q_epsilon(measurement_model),
  Z2_tilde = calculate(measurement_model, 'Z2_tilde', process_model)
) {
  if (!is.null(measurement_model[['gamma']])) {
    return(function(current, ...) {
      current$gamma <- measurement_model[['gamma']]
      current
    })
  }

  rg <- rgamma
  if ('lower' %in% names(measurement_model$gamma_prior)) {
    rg <- function(...) rtgamma(
      ...,
      lower = measurement_model$gamma_prior[['lower']],
      upper = measurement_model$gamma_prior[['upper']]
    )
  }

  attenuation_index <- as.integer(measurement_model$attenuation_factor)
  n_gamma <- nlevels(measurement_model$attenuation_factor)
  omega <- .make_omega(process_model, measurement_model)

  gamma_slice <- lapply(seq_len(n_gamma), function(i) {
    do.call(slice, tuning)
  })

  n_per_gamma <- table(measurement_model$attenuation_factor)

  part_ij <- rbind(
    cbind(seq_len(n_gamma), seq_len(n_gamma)),
    if (n_gamma > 1) t(combn(n_gamma, 2)) else NULL
  )
  included <- apply(part_ij, 1, function(ij) {
    n_per_gamma[ij[1]] > 0 && n_per_gamma[ij[2]] > 0
  })
  part_ij <- part_ij[included, , drop = FALSE]

  function(current, warming_up) {
    z <- as.vector(Z2_tilde - X %*% omega(current))

    z_part <- lapply(seq_len(n_gamma), function(i) {
      z_i <- z
      z_i[attenuation_index != i] <- 0
      z_i
    })
    s_z_part <- lapply(seq_len(n_gamma), function(i) {
      as.vector(solve(measurement_model$measurement_covariance, z_part[[i]]))
    })

    part_x <- apply(part_ij, 1, function(ij) {
      i <- ij[1]
      j <- ij[2]

      if (i == j) {
        as.numeric(crossprod(z_part[[i]], s_z_part[[i]]))
      } else {
        z_j <- z
        z_j[attenuation_index != j] <- 0
        as.numeric(
          crossprod(z_part[[i]], s_z_part[[j]])
          + crossprod(z_part[[j]], s_z_part[[i]])
        )
      }
    })

    for (i in seq_len(n_gamma)) {
      if (n_per_gamma[i] == 0) {
        current$gamma[i] <- rg(
          1,
          shape = measurement_model$gamma_prior[['shape']],
          rate = measurement_model$gamma_prior[['rate']]
        )
      } else {
        output <- gamma_slice[[i]](current$gamma[i], function(gamma_i) {
          current2 <- current
          current2$gamma[i] <- gamma_i

          log_prior_value <- log_prior(measurement_model, current2)
          if (!is.finite(log_prior_value)) return(log_prior_value)

          sqrt_gamma <- sqrt(current2$gamma)
          gamma_ij <- (
            sqrt_gamma[part_ij[, 1]] * sqrt_gamma[part_ij[, 2]]
          )
          (
            log_prior_value
            + 0.5 * sum(n_per_gamma * log(current2$gamma))
            - 0.5 * sum(gamma_ij * part_x)
          )
        }, learn = warming_up, include_n_evaluations = TRUE)

        log_trace(
          'gamma[{i}] took {output$n_evaluations} evaluations, w = {output$w}'
        )
        current$gamma[i] <- output$sample
      }
    }

    current
  }
}
