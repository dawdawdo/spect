flog.threshold("ERROR") 

test_that("plot_km does not return an error", {
  
  spect_train_data <- readRDS(test_path("testdata", "helper_train_data.rds"))
  
  spect_train_data <- readRDS(test_path("testdata", "helper_train_data.rds"))
  
  expect_error(plot_km(spect_train_data$result, 
                       prediction_threshold_search_granularity=0.2),
               NA)
  
})
