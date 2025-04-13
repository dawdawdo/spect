flog.threshold("ERROR")

test_that("plot_survivial_curve does not throw an error or warning", {

  rng_seed <- 42
  set.seed(rng_seed)
  
  syn_data <- create_synthetic_data(sample_size=500,
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

  id <- 1
  curve <- "both"
  
  expect_no_error(plot_survival_curve(result, individual_id=id, curve_type=curve))

})

test_that("throws an error when an invalid curve_type is passed", {
  
  rng_seed <- 42
  set.seed(rng_seed)
  
  syn_data <- create_synthetic_data(sample_size=500,
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
  
  id <- 1
  curve <- "survival"
  
  expect_error(plot_survival_curve(result, individual_id=id, curve_type=curve),
               regexp="is an invalid curve type")
  
})

test_that("throws an error when an invalid individual_id is passed", {
  
  rng_seed <- 42
  set.seed(rng_seed)
  
  syn_data <- create_synthetic_data(sample_size=500,
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
  
  id <- -1
  curve <- "both"
  
  expect_error(plot_survival_curve(result, individual_id=id, curve_type=curve),
               regexp="is an invalid id for plotting")
  
})