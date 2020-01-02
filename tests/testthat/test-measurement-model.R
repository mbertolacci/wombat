context('measurement-model')

test_that('measurement model outputs have correct dimensions', {
  process_model <- flux_process_model(emissions, control, sensitivities)
  model <- flux_measurement_model(
    soundings,
    ~ instrument_mode,
    soundings$sounding_id,
    process_model,
    attenuation_variables = 'instrument_mode'
  )
  for (name in c('soundings', 'gamma_prior')) {
    expect_false(is.null(model[[name]]))
  }
  for (name in c('beta')) {
    expect_true(is.null(model[[name]]))
  }
  n_beta <- 2
  expect_equal(dim(model$C), c(nrow(soundings), nrow(control)))
  expect_equal(dim(model$measurement_precision), rep(nrow(soundings), 2))
  expect_equal(dim(model$A), c(nrow(soundings), n_beta))
  expect_length(model$beta_prior_mean, n_beta)
  expect_equal(dim(model$beta_prior_precision), rep(n_beta, 2))
  expect_length(model$attenuation_factor, nrow(soundings))
})

test_that('flux measurement model samples have correct dimensions', {
  process_model <- flux_process_model(emissions, control, sensitivities)
  model <- flux_measurement_model(
    soundings,
    ~ instrument_mode,
    soundings$sounding_id,
    process_model,
    attenuation_variables = 'instrument_mode'
  )

  measurement_sample <- generate(model, generate(process_model))

  expect_length(measurement_sample$beta, 2)
  expect_length(measurement_sample$gamma, 2)
  expect_length(measurement_sample$Z2_tilde, nrow(soundings))
})
