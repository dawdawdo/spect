flog.appender(appender.file("unit_tests.log"))

test_that("multiplication works", {
  
  spect_train_data <- readRDS(test_path("testdata", "helper_train_data.rds"))
  
  prediction_times <- c(12, 24, 40)
  
  expect_error(eval <- evaluate_model(train_result=spect_train_data$result, 
                              prediction_times, 
                              plot_roc=FALSE)
               , NA)
  
  # Note, here we should probe eval for some known values
})
