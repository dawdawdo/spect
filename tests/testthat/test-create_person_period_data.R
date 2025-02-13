flog.threshold("ERROR")

test_that("create_person_period_data generates the correct number of rows", {
  
  train_data <- data.frame(a=c(1,2,3,4,5,6,7,8,9,10), 
                           b=c(10,9,8,7,6,5,4,3,2,1), 
                           survival_time=c(1,2,3,8,9,10,11,12,13,20), 
                           event_detected=c(1,1,0,1,0,1,1,1,0,1))
  bins <- 5
  window <- 13
  
  t_col <- "survival_time"
  e_col <- "event_detected"
  
  bounds <- generate_bounds(train_data, time_col=t_col, event_col=e_col, 
                            suggested_intervals=bins, obs_window=window)
  
  pp_data <- create_person_period_data(train_data, bounds)
  
  expected_rows <- nrow(train_data) * length(bounds$upper_bounds)
  actual_rows <- nrow(pp_data)
    
  expect_equal(actual_rows, expected_rows)
    
})
