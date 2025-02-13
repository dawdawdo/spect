flog.threshold("ERROR")

test_that("plot_synthetic_data runs with no errors", {
  
  syn_data <- create_synthetic_data(sample_size=5000, censor_percentage=0.1
                                    , max_censor_amount=2, perturbation_shift=6)
  
  expect_error(plot_synthetic_data(syn_data), NA)
    
})
