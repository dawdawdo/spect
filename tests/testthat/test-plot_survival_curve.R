flog.threshold("ERROR")

test_that("plot_survivial_curve does not throw an error or warning", {

  spect_train_data <- readRDS(test_path("testdata", "helper_train_data.rds"))

  result <- spect_train_data$result

  id <- 1
  curve <- "both"
  
  expect_no_error(plot_survival_curve(result, individual_id=id, curve_type=curve))

})

test_that("throws an error when an invalid curve_type is passed", {
  
  spect_train_data <- readRDS(test_path("testdata", "helper_train_data.rds"))
  
  result <- spect_train_data$result

  id <- 1
  curve <- "survival"
  
  expect_error(plot_survival_curve(result, individual_id=id, curve_type=curve),
               regexp="is an invalid curve type")
  
})

test_that("throws an error when an invalid individual_id is passed", {
  
  spect_train_data <- readRDS(test_path("testdata", "helper_train_data.rds"))
  
  result <- spect_train_data$result
  
  id <- -1
  curve <- "both"
  
  expect_error(plot_survival_curve(result, individual_id=id, curve_type=curve),
               regexp="is an invalid id for plotting")
  
})