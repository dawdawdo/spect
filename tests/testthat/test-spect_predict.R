flog.threshold("ERROR")

test_that("spect_predict data has the correct dimensions.", {

  rng_seed <- 42
  set.seed(rng_seed)
  
  syn_data <- create_synthetic_data(sample_size=2500,
                                    censor_percentage = 0.1,
                                    perturbation_shift = 6)
  
  
  source_data <- select(syn_data, -c(baseline_time_to_cancel, perturbed_baseline))
  
  predict_data <- source_data[1:10,]
  modeling_data <- source_data[11:nrow(source_data),]
  
  event_indicator_var <- "cancel_event_detected"
  survival_time_var <- "total_months"
  obs_window <- 48
  alg="glm"
  
  result <- spect_train(model_algorithm=alg, modeling_data=modeling_data,
                        event_indicator_var=event_indicator_var,
                        survival_time_var=survival_time_var,
                        obs_window=obs_window, use_parallel=FALSE)
  

  predictions <- spect_predict(result, predict_data)
  
  expect_equal(dim(predictions), c(100,12))
  
})
