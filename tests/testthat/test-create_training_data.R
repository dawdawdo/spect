flog.appender(appender.file("unit_tests.log"))

test_that("create_training_data captures data at the end of the last interval", {
  
  cens <- "prev"
  
  # Individual did not experience the event (event_detected = 0) and was last observed at the
  # ver end of interval 3 (survival_time = 9) so it should contribute 3 rows of data 
  # (0, 0, 0).
  person_period_data <- data.frame(
    individual_id=c(1, 1, 1),
    interval=c(1, 2, 3),
    lower_bound=c(0.0, 1.5, 4.0),
    upper_bound=c(1.5, 4.0, 9.0),
    covar=c(5, 10, 15),
    survival_time=c(9, 9, 9),
    event_detected=c(0, 0, 0))
  
  t_col <- "survival_time"
  e_col <- "event_detected"
  
  training_data <- create_training_data(person_period_data, t_col, e_col, cens)  
  
  events <- training_data$event_class
  expected_events <- as.character(c("X0", "X0", "X0"))
  
  expect_equal(events, expected_events)  
  
})

test_that("create_training_data raises an error when the censorship type is not in the list of options.", {

  cens <- "self" # Not on the list of available types
  
  # Individual experienced the event (event_detected = 1) in interval 3 (survival_time = 6)
  # so it should contribute 3 rows of data (0, 0, 1).
  person_period_data <- data.frame(
    individual_id=c(1, 1, 1),
    interval=c(1, 2, 3),
    lower_bound=c(0.0, 1.5, 4.0),
    upper_bound=c(1.5, 4.0, 9.0),
    covar=c(5, 10, 15),
    survival_time=c(6, 6, 6),
    event_detected=c(1, 1, 1))
  
  t_col <- "survival_time"
  e_col <- "event_detected"
  
  expect_error(create_training_data(person_period_data, t_col, e_col, cens),
               regexp=" Please choose from 'prev', 'same', and 'half'",
               fixed=TRUE)
})

test_that("create_training_data handles the case when the event was detected", {
  
  cens <- "same"
  
  # Individual experienced the event (event_detected = 1) in interval 3 (survival_time = 6)
  # so it should contribute 3 rows of data (0, 0, 1).
  person_period_data <- data.frame(
    individual_id=c(1, 1, 1),
    interval=c(1, 2, 3),
    lower_bound=c(0.0, 1.5, 4.0),
    upper_bound=c(1.5, 4.0, 9.0),
    covar=c(5, 10, 15),
    survival_time=c(6, 6, 6),
    event_detected=c(1, 1, 1))
  
  t_col <- "survival_time"
  e_col <- "event_detected"
  
  # Test censoring in an early period
  training_data <- create_training_data(person_period_data, t_col, e_col, cens)
  
  events <- training_data$event_class
  expected_events <- as.character(c("X0", "X0", "X1"))
  
  expect_equal(events, expected_events)  
  
})


test_that("create_training_data applies same-interval logic correctly", {

  cens <- "same"
  
  # Individual did not experience the event (event_detected = 0) and was censored at 
  # interval 1 (survival_time = 1) so it should contribute 1 row of data (0).
  person_period_data <- data.frame(
    individual_id=c(1, 1, 1),
    interval=c(1, 2, 3),
    lower_bound=c(0.0, 1.5, 4.0),
    upper_bound=c(1.5, 4.0, 9.0),
    covar=c(5, 10, 15),
    survival_time=c(1, 1, 1),
    event_detected=c(0, 0, 0))
  
  t_col <- "survival_time"
  e_col <- "event_detected"
  
  # Test censoring in an early period
  training_data <- create_training_data(person_period_data, t_col, e_col, cens)
  
  events <- training_data$event_class
  expected_events <- as.character(c("X0"))
  
  
  expect_equal(events, expected_events)
  
  
})

test_that("create_training_data applies previous-interval logic correctly", {
  
  cens <- "prev"
  
  # Individual did not experience the event (event_detected = 0) and was censored at interval 
  # 3 (survival_time = 5) so it should contribute 2 rows of data (0, 0).
  person_period_data <- data.frame(
    individual_id=c(1, 1, 1),
    interval=c(1, 2, 3),
    lower_bound=c(0.0, 1.5, 4.0),
    upper_bound=c(1.5, 4.0, 9.0),
    covar=c(5, 10, 15),
    survival_time=c(5, 5, 5),
    event_detected=c(0, 0, 0))
  
  t_col <- "survival_time"
  e_col <- "event_detected"
  
  # Test censoring in an early period
  training_data <- create_training_data(person_period_data, t_col, e_col, cens)
  
  events <- training_data$event_class
  expected_events <- as.character(c("X0", "X0"))
  
  expect_equal(events, expected_events)
  
})

test_that("create_training_data applies previous-interval logic correctly when no rows contribute.", {
  
  cens <- "prev"
  
  # Individual did not experience the event (event_detected = 0) and was censored at interval 
  # 1 (survival_time = 1) so it should contribute no data to the output.
  person_period_data <- data.frame(
    individual_id=c(1, 1, 1),
    interval=c(1, 2, 3),
    lower_bound=c(0.0, 1.5, 4.0),
    upper_bound=c(1.5, 4.0, 9.0),
    covar=c(5, 10, 15),
    survival_time=c(1, 1, 1),
    event_detected=c(0, 0, 0))
  
  t_col <- "survival_time"
  e_col <- "event_detected"
  
  # Test censoring in an early period
  training_data <- create_training_data(person_period_data, t_col, e_col, cens)
  
  expect_equal(nrow(training_data), 0)
  
})

test_that("create_training_data applies half-interval logic when the individual did not survive for half of the censored interval", {
  
  cens <- "half"
  
  # Individual did not experience the event (event_detected = 0) and was censored at 
  # interval 3 (survival_time = 6). Because it did not survive half of that interval, it
  # should contribute 2 rows of data (0, 0).
  person_period_data <- data.frame(
    individual_id=c(1, 1, 1),
    interval=c(1, 2, 3),
    lower_bound=c(0.0, 1.5, 4.0),
    upper_bound=c(1.5, 4.0, 9.0),
    covar=c(5, 10, 15),
    survival_time=c(6, 6, 6),
    event_detected=c(0, 0, 0))
  
  t_col <- "survival_time"
  e_col <- "event_detected"
  
  # Test censoring in an early period
  training_data <- create_training_data(person_period_data, t_col, e_col, cens)
  
  events <- training_data$event_class
  expected_events <- as.character(c("X0", "X0"))
  
  expect_equal(events, expected_events)
  
  
})

test_that("create_training_data applies half-interval logic when the individual survived at least half of the censored interval", {
  
  cens <- "half"
  
  # Individual did not experience the event (event_detected = 0) and was censored at 
  # interval 3 (survival_time = 8). Because it survived at least half of that interval, it
  # should contribute 3 rows of data (0, 0, 0).
  person_period_data <- data.frame(
    individual_id=c(1, 1, 1),
    interval=c(1, 2, 3),
    lower_bound=c(0.0, 1.5, 4.0),
    upper_bound=c(1.5, 4.0, 9.0),
    covar=c(5, 10, 15),
    survival_time=c(8, 8, 8),
    event_detected=c(0, 0, 0))
  
  t_col <- "survival_time"
  e_col <- "event_detected"
  
  # Test censoring in an early period
  training_data <- create_training_data(person_period_data, t_col, e_col, cens)
  
  events <- training_data$event_class
  expected_events <- as.character(c("X0", "X0", "X0"))
  
  expect_equal(events, expected_events)
  
  
})
