flog.threshold("ERROR") 

test_that("generate_bounds doesn't create more intervals than events", {

  # Note: There are only two events detected, but five intervals suggested.
  # We expect only 2 intervals to be created in this case.
  train_data <- data.frame(a=c(1,2,3,4), b=c(5,6,7,8), survival_time=c(1,4,5,6), 
                           event_detected=c(1,0,0,1))
  bins <- 5
  window <- 8
  
  bounds <- generate_bounds(train_data, time_col="survival_time", event_col="event_detected", suggested_intervals=bins, obs_window=window)
  
  events <- nrow(train_data[train_data$event_detected==1,])
  intervals <- length(bounds$upper_bounds)
  
  expect_equal(events, intervals)

})

test_that("generate_bounds uses suggested_bins appropriately", {
  
  train_data <- data.frame(a=c(1,2,3,4,5,6,7,8,9,10), 
                           b=c(10,9,8,7,6,5,4,3,2,1), 
                           survival_time=c(1,2,3,8,9,10,11,12,13,20), 
                           event_detected=c(1,1,0,1,0,1,1,1,0,1))
  
  bins <- 4
  window <- 20
  
  bounds <- generate_bounds(train_data, time_col="survival_time", event_col="event_detected", suggested_intervals=bins, obs_window=window)
  
  intervals <- length(bounds$upper_bounds)
  
  expect_equal(intervals, bins)
  
})


test_that("generate_bounds truncates the last interval at the observation window", {
  
  train_data <- data.frame(a=c(1,2,3,4,5,6,7,8,9,10), 
                           b=c(10,9,8,7,6,5,4,3,2,1), 
                           survival_time=c(1,2,3,8,9,10,11,12,13,20), 
                           event_detected=c(1,1,0,1,0,1,1,1,0,1))
  bins <- 5
  window <- 13
  
  bounds <- generate_bounds(train_data, time_col="survival_time", event_col="event_detected", suggested_intervals=bins, obs_window=window)
  
  last_bound <- bounds$upper_bounds[length(bounds$upper_bounds)]

  expect_equal(last_bound, window)
  
})

test_that("generate_bounds fails when the survival time column isn't in the data set", {
  
  train_data <- data.frame(a=c(1,2,3,4,5,6,7,8,9,10), 
                           b=c(10,9,8,7,6,5,4,3,2,1), 
                           survival_time=c(1,2,3,8,9,10,11,12,13,20), 
                           event_detected=c(1,1,0,1,0,1,1,1,0,1))
  
  bins <- 4
  window <- 20
  
  t_col <- "arrival_time" # Note: This column is not in the data
  e_col <- "event_detected"
  
  expect_error(generate_bounds(train_data, time_col=t_col, event_col=e_col,suggested_intervals=bins, obs_window=window), 
               regexp = "Could not find survivial time column", 
               fixed=TRUE)  
})

test_that("generate_bounds fails when the event indicator column isn't in the data set", {
  
  train_data <- data.frame(a=c(1,2,3,4,5,6,7,8,9,10), 
                           b=c(10,9,8,7,6,5,4,3,2,1), 
                           survival_time=c(1,2,3,8,9,10,11,12,13,20), 
                           event_detected=c(1,1,0,1,0,1,1,1,0,1))
  
  bins <- 4
  window <- 20
  
  t_col <- "survival_time"
  e_col <- "event_occurred" # Note: This column is not in the data

  expect_error(generate_bounds(train_data, time_col=t_col, event_col=e_col,
                                       suggested_intervals=bins, obs_window=window), 
               regexp = "Could not find event indicator column", 
               fixed=TRUE)
})



  