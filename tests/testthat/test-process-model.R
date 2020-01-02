context('process-model')

test_that('process model outputs have correct dimensions', {
  model <- flux_process_model(emissions, control, sensitivities)
  for (name in c('emissions', 'control', 'a_prior', 'w_prior')) {
    expect_false(is.null(model[[name]]))
  }
  for (name in c('a', 'w', 'eta')) {
    expect_true(is.null(model[[name]]))
  }
  n_alpha <- length(regions) * length(month_starts)
  expect_equal(dim(model$Phi), c(nrow(emissions), n_alpha))
  expect_equal(dim(model$H), c(nrow(control), n_alpha))
  expect_equal(nrow(model$Psi), nrow(control))
  expect_equal(dim(model$eta_prior_precision), rep(ncol(model$Psi), 2))
})

test_that('flux process model samples have correct dimensions', {
  model <- flux_process_model(emissions, control, sensitivities)
  process_sample <- generate(model)
  expect_length(process_sample$a, 1)
  expect_length(process_sample$alpha, ncol(model$H))
  expect_equal(
    dim(process_sample$Q_alpha),
    rep(ncol(model$H), 2)
  )
  expect_length(process_sample$eta, ncol(model$Psi))
  expect_length(process_sample$Y1_tilde, nrow(model$emissions))
  expect_length(process_sample$Y2_tilde, nrow(model$control))
})
