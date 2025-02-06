flog.appender(appender.file("unit_tests.log"))

test_that("spect_predict data has the correct dimensions.", {

  spect_train_data <- readRDS(test_path("testdata", "helper_train_data.rds"))
  
  result <- spect_train_data$result
  predict_data <- spect_train_data$predict_data
    
  predictions <- spect_predict(result, predict_data)
  
  expect_equal(dim(predictions), c(100,12))
  
})
