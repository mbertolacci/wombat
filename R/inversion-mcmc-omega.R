.make_omega_sampler <- function(
  measurement_model,
  process_model,
  variables = c('alpha', 'beta'),
  X = .make_X_omega(process_model, measurement_model, variables),
  X_c = .make_X_omega(
    process_model,
    measurement_model,
    setdiff(c('alpha', 'eta', 'beta'), variables)
  ),
  Q_epsilon = .make_Q_epsilon(measurement_model),
  Xt_Q_epsilon_X = .make_Xt_Q_epsilon_X(X, measurement_model),
  X_ct_Q_epsilon_X = .make_A_Q_epsilon_B(t(X_c), X, measurement_model),
  mu_omega = .make_mu_omega(
    process_model,
    measurement_model,
    variables
  ),
  Q_omega = .make_Q_omega(
    process_model,
    measurement_model,
    variables
  ),
  chol_Q_omega_conditional = .make_chol_Q_omega_conditional(
    process_model,
    measurement_model,
    variables,
    X = X,
    Xt_Q_epsilon_X = Xt_Q_epsilon_X,
    Q_omega = Q_omega
  ),
  chol_Q_omega_c_conditional = .make_chol_Q_omega_conditional(
    process_model,
    measurement_model,
    setdiff(c('alpha', 'eta', 'beta'), variables),
    X = X_c
  ),
  Z2_tilde = calculate(measurement_model, 'Z2_tilde', process_model)
) {
  omega_unpack <- .make_omega_unpack(
    process_model,
    measurement_model,
    variables
  )

  if (is.null(X_c)) {
    # Nothing has been marginalised out
    function(current) {
      chol_Q_omega_conditional_i <- chol_Q_omega_conditional(current)
      mu_omega_conditional <- as.vector(.chol_solve(
        chol_Q_omega_conditional_i,
        crossprod(X, Q_epsilon(current) %*% Z2_tilde)
        + Q_omega(current) %*% mu_omega(current)
      ))
      omega <- (
        mu_omega_conditional
        + .sample_normal_precision_chol(chol_Q_omega_conditional_i)
      )
      omega_unpack(current, omega)
    }
  } else {
    # Something has been marginalised out; use Woodbury
    function(current) {
      chol_Q_omega_c_conditional_i <- chol_Q_omega_c_conditional(current)
      Q_epsilon_i <- Q_epsilon(current)
      chol_Q_omega_conditional <- chol(
        Q_omega(current)
        + Xt_Q_epsilon_X(current)
        - crossprod(solve(
          t(chol_Q_omega_c_conditional_i),
          X_ct_Q_epsilon_X(current)
        ))
      )
      mu_omega_conditional <- as.vector(.chol_solve(
        chol_Q_omega_conditional,
        crossprod(X, Q_epsilon_i %*% (
          Z2_tilde - X_c %*% .chol_solve(
            chol_Q_omega_c_conditional_i,
            crossprod(X_c, Q_epsilon_i %*% Z2_tilde)
          )
        ))
        + Q_omega(current) %*% mu_omega(current)
      ))
      omega <- (
        mu_omega_conditional
        + .sample_normal_precision_chol(chol_Q_omega_conditional)
      )
      omega_unpack(current, omega)
    }
  }
}
