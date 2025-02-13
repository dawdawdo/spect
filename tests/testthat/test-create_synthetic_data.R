
# Note: There is currently no test that the plot_synthetic_data function actually generates
# the expected plot. Because it is a bit of a burden to add on multiple platforms and is 
# relatively simple to implement as a user if absolutely necessary, it is only minimally
# probed here.

flog.threshold("ERROR")

test_that("create synthetic data creates a data with the right columns", {
  data <- create_synthetic_data()
  
  expected_columns <- c("incomes", "watchtimes", "total_months", "cancel_event_detected"
                        , "baseline_time_to_cancel", "perturbed_baseline")
  expect_named(data, expected_columns)
})

test_that("syntehtic data has the right sample size", {
  rowcount <- 100
  data <- create_synthetic_data(sample_size = rowcount)
  expect_vector(data, size=100)
})




