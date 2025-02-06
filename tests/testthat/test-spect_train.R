flog.appender(appender.file("unit_tests.log"))

test_that("spect_train return data dimensions are correct.", {
 
  rng_seed <- 42
  set.seed(rng_seed)
  event_indicator_var <- "cancel_event_detected"
  survival_time_var <- "total_months"
  obs_window <- 48
  alg="glm"
  
  syn_data <- create_synthetic_data(sample_size=2500,
                                    censor_percentage = 0.1,
                                    perturbation_shift = 6)
  
  
  source_data <- select(syn_data, -c(baseline_time_to_cancel, perturbed_baseline))
  
  predict_data <- source_data[1:10,]
  modeling_data <- source_data[11:nrow(source_data),]
  
  
  result <- spect_train(model_algorithm=alg, modeling_data=modeling_data,
                        event_indicator_var=event_indicator_var,
                        survival_time_var=survival_time_var,
                        obs_window=obs_window, use_parallel=FALSE)


  expect_equal(dim(result$data.pp), c(19920, 8))
  expect_equal(dim(result$data.train), c(1992, 4))
  expect_equal(dim(result$data.modeling), c(12063, 4))
  expect_equal(dim(result$data.test_predictions), c(4980, 12))
  expect_equal(dim(result$data.test), c(498, 4))
  
  expect_equal(length(result$bounds$upper_bounds), 10)
  expect_equal(length(result$bounds$lower_bounds), 10)
  
  
})
