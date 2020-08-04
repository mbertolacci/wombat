.make_omega_sampler <- function(
  measurement_model,
  process_model,
  X = .make_X_omega(process_model, measurement_model),
  solve_Q_epsilon = .make_solve_Q_epsilon(measurement_model),
  Xt_Q_epsilon_X = .make_Xt_Q_epsilon_X(X, measurement_model),
  mu_omega = .make_mu_omega(process_model, measurement_model),
  Q_omega = .make_Q_omega(process_model, measurement_model),
  chol_Q_omega_conditional = .make_chol_Q_omega_conditional(
    process_model,
    measurement_model,
    X = X,
    Xt_Q_epsilon_X = Xt_Q_epsilon_X,
    Q_omega = Q_omega
  ),
  Z2_tilde = calculate(measurement_model, 'Z2_tilde', process_model)
) {
  omega_unpack <- .make_omega_unpack(process_model, measurement_model)

  function(current) {
    chol_Q_omega_conditional_i <- chol_Q_omega_conditional(current)
    mu_omega_conditional <- as.vector(.chol_solve(
      chol_Q_omega_conditional_i,
      crossprod(X, solve_Q_epsilon(Z2_tilde, current))
      + Q_omega(current) %*% mu_omega(current)
    ))
    omega <- (
      mu_omega_conditional
      + .sample_normal_precision_chol(chol_Q_omega_conditional_i)
    )
    omega_unpack(current, omega)
  }
}
