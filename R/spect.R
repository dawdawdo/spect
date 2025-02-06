# spect - survival predictive ensemble classification tool
#
# Purpose: A tool for ensemble classification modeling of survival data which 
# provides:
# - Simplicity - A dead-simple interface, requiring minimal parameter specification
# and tooling for visualization.
# - Flexibilty - users may provide advanced parameters, if desired
# - Robustness - tested against synthetic and standard survival benchmark data
# - Transparency - provides all intermediate data sets and parameters to the user
# for further analysis

# Not explicitly loaded with "library()", but required for the spect_train fucntion to run.
#' @import caret
#' @import caretEnsemble

##############################################################################
# fucntions generate_bounds, create_person_period_data, create_training_data, 
# and spect_predict adapted from https://github.com/ksuresh17/autoSurv
##############################################################################


#' @import futile.logger
library(futile.logger)

#' @import dplyr
# Note: This is only required to pass the R CMD check.
#' @importFrom rlang :=
library(dplyr)

#' Generates a survival data set for synthetic streaming service subscription data. 
#' The survival event in this case is a cancellation of the subscription. It is 
#' given as a function of household income and average number of hours watched 
#' in the prior month. Users can adjust the level of censoring and variance in 
#' the data with the supplied parameters or simply call with no parameters for a
#' default distribution of data.
#'
#' @param sample_size optional - size of the sample population to generate
#' @param minimum_income optional - minimum household income used to generate the distribution
#' @param median_income optional - median household income used to generate the distribution
#' @param income_variance optional - variance to use when generating the household income distribution
#' @param min_watchhours optional - minimum average number of hours watched used to generate the distribution
#' @param max_watchhours optional - minimum average number of hours watched used to generate the distribution
#' @param censor_percentage optional - percentage of population to artificially censor
#' @param min_censor_amount optional - Minimum number of months of censoring to apply to the censored population 
#' @param max_censor_amount optional - maximum number of months of censoring to apply to the censored population
#' @param study_time_in_months optional - observation horizon in months
#' @param perturbation_shift optional - defines a boundary for the amount to randomly perturb the formulaic result. Zero for no perturbation
#'
#' @return A survival data set suitable for modeling using spect_train.
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords utilities
#'
#' @examples
#' data <- create_synthetic_data()
#'
#' @export
create_synthetic_data <- function(sample_size = 250
                          , minimum_income = 5000
                          , median_income = 50000
                          , income_variance = 10000
                          , min_watchhours = 0.0
                          , max_watchhours = 6.0
                          , censor_percentage = 0.0
                          , min_censor_amount = 0.0
                          , max_censor_amount = 0.0
                          , study_time_in_months = 48
                          , perturbation_shift = 0
                          )
{
  
  flog.info("Creating %i income samples from normal distribution of median %i, variance %i 
            and watchtimes samples from uniform distribution with min: %i and max: %i", 
             sample_size, median_income, income_variance, min_watchhours, max_watchhours)

  income_candidates <- stats::rnorm(sample_size, median_income, income_variance)
  incomes <- ifelse(income_candidates < minimum_income, minimum_income, income_candidates)
  watchtimes <- stats::runif(sample_size, min_watchhours, max_watchhours)
  
  
  # Apply the formula: 26 + watch_time^2 - (Income/10000)
  # By inspection, this seems reasonable
  flog.debug("Generating time to cancel with random perturbation of %i...", perturbation_shift)

  baseline_time_to_cancel <- 26 + watchtimes * watchtimes - (incomes/10000)
  perturbed_baseline <- baseline_time_to_cancel + stats::runif(sample_size, -perturbation_shift, perturbation_shift)
  
  flog.debug("Applying random censorship to totral_months...")
  censor_selector <- ifelse(stats::runif(sample_size, 0, 1) < censor_percentage, 1, 0)
  censor_amount <- stats::runif(sample_size, min_censor_amount, max_censor_amount)
  censored_baseline <- perturbed_baseline - censor_selector * censor_amount
  
  total_months <- ifelse(censored_baseline > study_time_in_months, study_time_in_months, censored_baseline)
  
  flog.debug("Generating event indicators...")
  event_within_study_time <- ifelse(perturbed_baseline <= study_time_in_months, 1, 0)
  cancel_event_detected <- event_within_study_time * (1 - censor_selector)
  
  modeling_data <- data.frame(incomes, watchtimes, total_months, cancel_event_detected, 
                              baseline_time_to_cancel, perturbed_baseline)
  
  return(modeling_data)
  
}

#' Simple visualization of synthetic subscription data.
#'
#' @param data a data object generated by create_synthetic_data
#' 
#' @return None - prints synthetic data generated by create_synthetic_data
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords utilities
#'
#' @examples
#' data <- create_synthetic_data()
#' plot_synthetic_data(data)
#'
#' @export
plot_synthetic_data <- function(data)
{
  
  # Note - the use of .data in the aes function here avoids the "no visible binding for 
  # global variable" note in R CMD check. As the time of this writing, this appears to be
  # the only solution that doesn't raise the check NOTE and doesn't trigger a warning from
  # the plot function itself.
  flog.debug("Generating incomes plot...")
  print(ggplot2::ggplot(data, ggplot2::aes(x=.data$incomes)) + ggplot2::geom_density())

  flog.debug("Generating watchtimes plot...")
  print(ggplot2::ggplot(data, ggplot2::aes(x=.data$watchtimes)) + ggplot2::geom_density())

  flog.debug("Generating cancellation times plot...")
  print(ggplot2::ggplot(data, ggplot2::aes(x=.data$baseline_time_to_cancel)) + ggplot2::geom_density())

  flog.debug("Generating perturbed cancellation times plot...")
  print(ggplot2::ggplot(data, ggplot2::aes(x=.data$perturbed_baseline)) + ggplot2::geom_density())

  flog.debug("Generating final months plot...")
  print(ggplot2::ggplot(data, ggplot2::aes(x=.data$total_months)) + ggplot2::geom_density())

}

#' Generates the intervals based on the survival times in the supplied data set 
#' using the quantile function.
#'
#' @param train_data A survival data set containing at least three columns - one which
#' matches the string in the `time_col` parameter, one which matches the string in the 
#' `event_col` parameter, and at least one covariate column for modeling.
#' @param time_col The name of the column in `train_data` containing survivial time 
#' @param event_col The name of the column in `train_data` contaiing the event indicator. 
#' Values in this column must be either zero (0) or one (1)
#' @param suggested_intervals The number of intervals to create. If the number of events
#' in the data is less than `suggested_intervals`, it is ignored.
#' @param obs_window An artificial censoring time. Any observations in `train_data` beyond
#' this time will be administratively censored.
#' 
#' @return A list of upper an lower bounds for each generated interval.
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords utilities
#' 
#' @examples
#' df <- data.frame(a=c(1,2,3,4,5,6), surv_time=c(1,4,5,6,8,9), event=c(1,1,1,1,0,1))
#' bounds <- generate_bounds(df, time_col="surv_time", event_col="event", 
#'                           suggested_intervals=3, obs_window=8)
#'                           
#' @seealso [create_person_period_data()]
#'
#' @export
generate_bounds <- function(train_data, time_col, event_col, suggested_intervals, obs_window){
  
  flog.debug("Checking data for survivial time column, %s...", time_col)
  if(!(time_col %in% colnames(train_data))){
    msg <- sprintf("Could not find survivial time column %s in the training data. Are you sure there's a column with that name?", time_col)
    stop(msg)
  }
  
  flog.debug("Checking data for event indicator column, %s...", event_col)
  if(!(event_col %in% colnames(train_data))){
    msg <- sprintf("Could not find event indicator column %s in the training data. Are you sure there's a column with that name?", event_col)
    stop(msg)
  }
  
  time_col_index <- which(eval(time_col) == colnames(train_data))
  event_col_index <- which(eval(event_col) == colnames(train_data))
  
  flog.debug("administratively censoring beyond time %i...", obs_window)
  train_data[train_data[[time_col_index]] > obs_window, event_col_index] <- 0
  train_data[train_data[[time_col_index]] > obs_window, time_col_index] <- obs_window
  
  
  time <- train_data[,time_col_index]
  event <- train_data[,event_col_index]
  calculated_intervals <- min(suggested_intervals, length(unique(time[event==1])))

  flog.debug("Creating %i intervals...", calculated_intervals)  
  slice_indexes <- seq(from=0, to=1, length.out=(calculated_intervals + 1))
  raw_upper_bounds <- stats::quantile(time[event==1], probs=slice_indexes, names=FALSE)[-1]
  
  flog.debug("Forcing bounds between zero and %i...", obs_window)
  upper_bounds <- c(raw_upper_bounds[-length(raw_upper_bounds)], obs_window)
  lower_bounds <- c(0, upper_bounds[-length(upper_bounds)])
  
  return (mget(c("upper_bounds", "lower_bounds")))
}

#' Generates person-period data for any data set, given the bounds defined by the 
#' training set.
#'
#' @param individual_data A survival data set.
#' @param bounds Output from the `generate_bounds` function of this package.
#' 
#' @return A data set consisting of the original `individual_data` repeated once for each
#' interval defined by the `bounds` parameter. Each row will be labeled with an id and
#' an interval. The output of this function can be passed to either `create_training_data`
#' or `spect_predict` to genreate modeling data or predictions respectively.
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords utilities
#' @seealso [generate_bounds()], [spect_predict()], [create_training_data()]
#'
#' @export
create_person_period_data <- function(individual_data, bounds){
  
  flog.debug("Generating unique ID for each individual in the population...")  
  periods <- length(bounds$upper_bounds)
  individuals <- nrow(individual_data)
  individual_id <- rep(seq_len(individuals), each=periods) 
  
  flog.debug("Generating interval ids and bounds...")  
  interval_data <- cbind(interval=1:periods, lower_bound=bounds$lower_bounds, upper_bound=bounds$upper_bounds)
  interval_columns <- interval_data[rep(seq_len(nrow(interval_data)), 
                                        times=individuals), ]
  individual_period_data <- individual_data[rep(seq_len(nrow(individual_data)), 
                                                each=periods),]
  
  return (cbind(individual_id, interval_columns, individual_period_data))
}

#' Generates modeling data from a person-period data set.
#'
#' @param person_period_data A discrete-time data set. Generally, this will be output
#' from the `create_person_period_data` function.
#' @param time_col A string specifying the name of the column which contains the survival 
#' time.
#' @param event_col A string specifying the name of the column which contains the event
#' indicator.
#' @param cens Specifies how to apply censored data. Valid values are "same" - considers
#' censorship to occur in the same interval as the survival time, "prev" - considers
#' censorship to occur in the prior interval, and "half" - considers censorship to occur
#' in the same interval as survival time if the individual survived for at least half of
#' that interval.
#' 
#' @return A discrete-time data set suitable for training using any binary classifer.
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords utilities
#' @seealso [create_person_period_data()]
#'
#' @export
create_training_data <- function(person_period_data, time_col, event_col, cens) {
  
  flog.debug("Testing censorship type...")
  censorship_types <- c("prev", "same", "half")
  if(!(cens %in% censorship_types)){
    msg <- sprintf("Censorship type %s is invalid. Please choose from 'prev', 'same', and 'half'", cens)
    stop(msg)
  }

  flog.debug("applying '%s' censorship algorithm to person-period data...", cens)
  top_bound <- max(person_period_data$upper_bound)
  
  training_data_candidate <- person_period_data %>% mutate(
    event_interval_indicator = case_when(
      !!sym(time_col) > upper_bound ~ 0,
      !!sym(time_col) <= lower_bound ~ NA,
      !!sym(event_col) == 1 ~ 1,
      cens == "same" ~ 0,
      cens == "prev" ~ NA,
      cens == "half" & !!sym(time_col) >= 0.5 * (upper_bound + lower_bound) ~ 0,
      cens == "half" & !!sym(time_col) < 0.5 * (upper_bound + lower_bound) ~ NA,
      TRUE ~ -1)) %>% mutate(
        event_interval_indicator = case_when(
          !!sym(time_col) == top_bound & is.na(event_interval_indicator) ~ 0,
          TRUE ~ event_interval_indicator
        )) %>% filter(!is.na(event_interval_indicator)) %>% arrange(individual_id, interval)
  
  flog.debug("Setting event class as a factor...")
  training_data_candidate$event_class <- make.names(training_data_candidate$event_interval_indicator)
  
  flog.debug("Stripping metadata for modeling...")
  training_data <- training_data_candidate %>% dplyr::select(
    -c(event_interval_indicator, individual_id, !!sym(time_col), !!sym(event_col), 
       upper_bound, lower_bound))
  
  return(training_data)
}

#' Generates a trained caret model using the given primary binary classification. Optionally
#' generates a stacked ensemble model if a list of base learners is supplied.
#'
#' @param test_prop optional proportion of the data set to reserve for testing
#' @param censor_type optional method used to determine censorship in a given bin - may be "half", "prev" or "same". see createDiscreteDat for usage.
#' @param bin_slices optional number of intervals to use for predictions.
#' @param method optional caret parameter
#' @param resampling_number optional for repeated cv 
#' @param kfold_repeats optional number of folds
#' @param model_algorithm primary classification algorithm. Trains a stack-ensemble 
#' model if `base_learner_list` is supplied, otherwise trains a simple classifier model.
#' @param base_learner_list optional list of base learner algorithms
#' @param metric optional metric for model calibration
#' @param rng_seed optional random number generation seed for reproducibility
#' @param use_parallel optioanlly make use of the caret multicore training cluster
#' @param cores optioanl number of cores for multicore training. If zero, spect will 
#' attempt to make a good choice. Note: only relevant if `use_parallel` is set to TRUE, 
#' otherwise this parameter is ignored.
#' @param modeling_data This data set must have one column for time and one column for the 
#' event indicator. The remaining columns are treated as covariates for modeling.
#' @param event_indicator_var The name of the column containing the event indicator 
#' (values in this column must be zero or one).
#' @param survival_time_var The name of the column containing the time variable
#' @param obs_window The last time to use for generating person-period data. 
#' Any event occurring after this time will be administratively censored. In general, 
#' choosing a time at or near the end of the max observed time will include most events.
#'
#' @return A list containing all intermediate data sets created by `spect_train`, a 
#' trained caret model object, the following parameters passed to `spect_train`: 
#' `obs_window`, `survival_time_var`, `event_indicator_var`, `base_learner_list`, 
#' `bin_slices`, and the bounds of each interval generated by the training data set.
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#'
#' @export
#'
spect_train <- function(
  test_prop = 0.2
  , censor_type = 'half'
  , bin_slices = 10
  , method = 'repeatedcv' 
  , resampling_number = 10 
  , kfold_repeats = 3
  , model_algorithm
  , base_learner_list = list()
  , metric = 'Kappa'
  , rng_seed = 42
  , use_parallel = TRUE
  , cores = 0
  , modeling_data
  , event_indicator_var
  , survival_time_var
  , obs_window
)
{
  
  flog.debug("Setting RNG seed to %i", rng_seed)
  set.seed(rng_seed)  
  
  flog.info("Splitting test/train data at %f/%f...", test_prop, (1 - test_prop))
  testIndexes <- sample(x = 1:nrow(modeling_data), size = nrow(modeling_data) * test_prop)
  
  flog.debug("Creating test data...")
  data.test <- modeling_data[testIndexes, ]
  
  flog.debug("Creating training data...")
  data.train <- modeling_data[-testIndexes, ]        
  
  flog.debug("Split complete!")
  
  flog.info("Creating person-period data set...")
  bounds <- generate_bounds(data.train, survival_time_var, event_indicator_var, bin_slices, obs_window)
  
  data.pp <- create_person_period_data(data.train, bounds)
  data.modeling <- create_training_data(data.pp, survival_time_var, 
                                        event_indicator_var, censor_type)
  
  if(use_parallel){
    
    flog.info("Parallel processing specified - creating cluster for Caret modeling...")
    machine_cores <- parallel::detectCores()
    used_cores <- max(3 * machine_cores / 4, 1)
    
    flog.debug("Using %i of %i detected cores for parallel processing", used_cores, machine_cores)
    cl <- parallel::makePSOCKcluster(used_cores)
    doParallel::registerDoParallel(cl)
  }
  
  flog.info("Creating caret train control...")
  modelControl <- caret::trainControl(method=method, number=resampling_number, 
                                      repeats=kfold_repeats, savePredictions=TRUE, 
                                      classProbs=TRUE, allowParallel = use_parallel)
  if (length(base_learner_list) == 0){
    flog.info('Training %s using %s method with %i resamples and %i kfold repeats'
              , model_algorithm, method, resampling_number, kfold_repeats)
    results <- caret::train(event_class ~ ., data=data.modeling, method=model_algorithm, metric=metric, trControl=modelControl)
  } else {
    flog.debug("Pruning algorithm list...")
    algorithm_list <- base_learner_list[!is.na(base_learner_list)]
    
    flog.info("Training base models %s using %s method with %i resamples and %i kfold repeats"
              , paste(algorithm_list, collapse = ", "), method, resampling_number, kfold_repeats)

    base_models <- caretEnsemble::caretList(event_class ~ ., data=data.modeling, 
                                            trControl=modelControl, 
                                            methodList=algorithm_list)
    
    flog.info("Training the metalearner %s using %s metric...", model_algorithm, metric)
    results <- caretEnsemble::caretStack(base_models, method=model_algorithm, 
                                         metric=metric, trControl=modelControl)
  }
  
  if(use_parallel){
    flog.info("Shutting down parallel computing cluster...")
    parallel::stopCluster(cl)
    rm(cl)
  }
  
  flog.info("Calculating probabilities on %i test individuals...", nrow(data.test))

  train_result <- mget(c("results", "bounds", "base_learner_list"))
  
  # Technical note: the train_result parameter in spect_predict need only contain 
  # a trained caret binary classification model and a list of upper and lower bounds 
  # in order to generate predictions. For this reason, the function is used here to 
  # generate predictions for the test data set, but can also be called by external 
  # users, passing the return value from this function. The tryCatch clause avoids
  # destroying the unit-testability of this construction.
  
  data.test_predictions <- list()
  
  tryCatch(
    {
      data.test_predictions <- spect_predict(train_result, data.test)
    },
    error=function(e) {
      message('An Error Occurred')
      print(e)
      message("Training is complete, but data.test_predictions will be empty.")
    },
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
  )
  
  return (mget(c("data.test", "data.train", "data.pp", "data.modeling", 
                 "results", "data.test_predictions", "event_indicator_var", 
                 "survival_time_var", "bin_slices", "bounds", 
                 "obs_window", "base_learner_list")))
}

#' Plots a sample of individual survival curves from the test data set.
#'
#' @param train_result return data object from `spect_train`
#' @param individual_id identifier of the individual to plot
#' @param curve_type optional specification of the type of curve. Available options are 
#' "conditional", which plots the conditional probability of surviving each interval given
#' that the individual survived to the start of that interval, "absolute" which plots the 
#' unconditional probability of surviving each interval, and "both", the default value, 
#' which plots both curves on the same chart.
#' 
#' @return None - plots the number of requested samples
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords visualization
#'
#' @export
#'
plot_survival_curve <- function(train_result, individual_id, curve_type="both"){
  
  flog.debug("Checking that %s is a valid curve type...", curve_type)
  if(!(curve_type %in% c("both", "conditional", "absolute"))){
    msg <- sprintf("%s is an invalid curve type. Valid options are 'conditional', 'absolute', and 'both.'", curve_type)
    stop(msg)
  }

  flog.debug("Checking that %i is a valid individual_id...", individual_id)
  if(!between(individual_id, 1, nrow(train_result$data.test_predictions))){
    msg <- sprintf("%i is an invalid id for plotting. Is it larger than the size of the test data set?", individual_id)
    stop(msg)
  }
  
  event_var <- train_result$event_indicator_var
  time_var <- train_result$survival_time_var
  pred_times <- train_result$obs_window
  bin_slices <- train_result$bin_slices
  
  data <- train_result$data.test_predictions
  
  individual_predictions <- data[data$individual_id == individual_id,]
  x_max <- max(individual_predictions$upper_bound)
  
  flog.debug("Collecting actuals for the plot label...")
  actual_absolute_time <- round(unique(individual_predictions[, time_var]))
  actual_interval <- round(actual_absolute_time / (pred_times / bin_slices))
  event_detected <- unique(individual_predictions[, event_var])
  
  title_value = paste("Survival Curve: Individual", individual_id)
  subtitle_value <- paste("Actual Time:", actual_absolute_time, "Actual Interval:", actual_interval,"Event:", event_detected)
  
  flog.debug("Plotting individual %i", individual_id)    
  plot(0,0,xlim = c(0,x_max),ylim = c(0,1),type = "n", xlab="", ylab="")
  
  if (curve_type %in% c("both", "conditional")){
    graphics::lines(individual_predictions$upper_bound, individual_predictions$cond_surv_pred, type='b', col='#FFA500', lwd="4")
  }
  
  if (curve_type %in% c("both", "absolute")){
    graphics::lines(individual_predictions$upper_bound, individual_predictions$abs_surv_pred, type='b', col='#1E90FF', lwd="4")
  }
  
  graphics::title(main = title_value, sub = subtitle_value,
                  xlab = "Time Interval", ylab = "Survival Probability",
                  cex.main = 2,   font.main= 4, col.main= "#808080",
                  cex.sub = 1, font.sub = 4, col.sub = "#800000",
                  col.lab = "#808080", 
  )
}

#' Plots a series of population Kaplan-Meier curves for different thresholds for both 
#' the test predictions and the ground truth
#'
#' @param train_result return data object from `spect_train`
#' @param prediction_threshold_search_granularity optional number between zero and one 
#' which defines the granularity of searching for cumulative probability thresholds. 
#' For instance, search a value of 0.05 will search 19 thresholds (0.05, 0.10, ..., 0.95)
#' 
#' @return Data used to produce the KM curve and the passed granularity parameter. Also
#' plots the KM curves.
#' 
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords visualization
#'
#' @export
#'
plot_km <- function(train_result, prediction_threshold_search_granularity = 0.05){
  
  flog.info("Calculating plot increments...")
  survival_time_var <- train_result$survival_time_var
  event_indicator_var <- train_result$event_indicator_var
  
  start <- prediction_threshold_search_granularity
  end <- 1 - prediction_threshold_search_granularity
  increments <- (1 / prediction_threshold_search_granularity) - 1
  
  flog.debug("Store the prediction vs. groud truth for each interval...")
  comp_data <- list()
  
  flog.debug("Checking thresholds...")
  for (prediction_threshold in seq(start, end, length.out = increments)){
    
    flog.info("collecting the minimum interval which crosses the prediction threshold: %f", prediction_threshold)
    
    # Note: The minimum interval is just a choice. We could have arbitrarily chosen
    # the middle of the interval or the lower bound, or even decided dynamically 
    # based on the distribution.
    prediction_interval <- train_result$data.test_predictions[
      train_result$data.test_predictions$abs_surv_pred < prediction_threshold,]
    
    prediction_interval$interval <- as.numeric(prediction_interval$interval)
    threshold_bound <- do.call(rbind, 
                               unname(by(prediction_interval, 
                                         prediction_interval$individual_id, 
                                         function(x) x[x$interval == min(x$interval),])))
    
    flog.info("Found %i individuals...", nrow(threshold_bound))
    
    # Note: join back to the original data to pick up the lower bound of that interval.
    # This will be used to create a linear interpolation.
    # This actually removes anything for which the event occurred in the first
    # interval. It won't happen often for large enough data sets, but should probably
    # be addressed.
    under_threshold <- inner_join(train_result$data.test_predictions
                                  , threshold_bound
                                  , by=c("individual_id", 
                                         "next_interval" = "interval"))
    
    flog.debug("Interpolating to find exact threshold crossing time...")
    # Note: This is the "true" prediction time measured in fractional bins, 
    # which can be translated to the original time unit, if desired.
    under_threshold$true_interval <- with( 
      under_threshold, interval + ((prediction_threshold - abs_surv_pred.y) 
                                   / (abs_surv_pred.x - abs_surv_pred.y))) 
    
    flog.debug("Finding individuals which never cross the threshold...")
    
    removals <- unique(under_threshold$individual_id)
    over_threshold <- train_result$data.test_predictions %>% filter(
      !(.data$individual_id %in% removals))
    
    over_threshold <- over_threshold[over_threshold$interval == train_result$obs_window, ]
    
    flog.debug("Found %i individuals...", nrow(over_threshold))
    
    alt_event_indicator_var <- paste(event_indicator_var, ".x", sep="")
    alt_survival_time_var <- paste(survival_time_var, ".x", sep="")
    
    flog.debug("Combine predicted and ground truth data sets...")
    threshold_hits_under <- under_threshold %>% 
      select(predicted_crossing_time = true_interval
             , !!survival_time_var := {{alt_survival_time_var}}
             , !!event_indicator_var := {{alt_event_indicator_var}})
    
    threshold_hits_over <- over_threshold %>% 
      select(predicted_crossing_time = interval,
             all_of(event_indicator_var),
             all_of(survival_time_var))
    
    threshold_hits <- rbind(threshold_hits_under, threshold_hits_over)
    
    flog.debug("Interpolating ground truth crossing time...")
    # Note: convert the real survival time to a fractional interval for the "Ground Truth" crossing time.
    threshold_hits$GT_crossing_time <- (threshold_hits[, survival_time_var] 
                                        * (1.0 * train_result$bin_slices)
                                        / (1.0 * train_result$obs_window) )
    
    predicted_elements <- select(threshold_hits, 
                                 c(predicted_crossing_time,
                                   all_of(event_indicator_var)))
    predicted_elements$source <- "Predictions"
    colnames(predicted_elements)[1] ="TimeSlice"
    
    ground_truth_elements <- select(threshold_hits, 
                                    c(GT_crossing_time, 
                                      all_of(event_indicator_var)))
    ground_truth_elements$source <- "Ground Truth"
    colnames(ground_truth_elements)[1] ="TimeSlice"
    
    gt_comp_data <- rbind(predicted_elements, ground_truth_elements)
    colnames(gt_comp_data)[2] = "EventIndicator"
    
    gt_comp_data$TimeSlice <- as.double(gt_comp_data$TimeSlice)
    
    flog.info("Generate P-value...")    
    
    KM_diff <- survival::survdiff(survival::Surv(TimeSlice, EventIndicator) ~ source, data=gt_comp_data)
    pval<-stats::pchisq(KM_diff$chisq, length(KM_diff$n)-1, lower.tail = FALSE)
    
    display_pval <- format(round(pval, 2), nsmall = 2)
    display_prediction_threshold <- format(round(prediction_threshold, 2), 
                                           nsmall = 2)
    flog.info("Pval is %s for threhold %s",
              display_pval, 
              display_prediction_threshold)
    
    comp_info <- list(prediction_threshold = prediction_threshold
                      , pval = pval
                      , KM_diff = KM_diff
                      , gt_comp_data = gt_comp_data)
    
    comp_data <- rbind(comp_data, comp_info)
    
  }
  
  flog.info("Initializing Kaplan-Meier curve plot...")
  
  splots <- list()
  
  flog.info("Collecting data for each threshold value...")  
  
  for (i in 1:nrow(comp_data)){
    
    flog.debug(paste("Preparing ground truth fit data for plot...iteration",
                     as.character(i)))
    
    threshold_data <- comp_data[i,]
    
    fmt_pval <- format(round(threshold_data$pval, 2), nsmall = 2)
    fmt_prediction_threshold <- format(
      round(threshold_data$prediction_threshold, 2),
      nsmall = 2)
    
    threshold_params <- paste("Prediction Threhsold: ", 
                              fmt_prediction_threshold, 
                              " Pval: ",  
                              fmt_pval)
    
    flog.debug(paste("Using the following parameters - ", threshold_params))
    
    gt <- threshold_data$gt_comp_data
    
    KM_fit <- survminer::surv_fit(
      survival::Surv(TimeSlice, EventIndicator) ~ source, 
      data=gt)
    
    splots[[i]] <- survminer::ggsurvplot(KM_fit, 
                                         conf.int = TRUE, 
                                         pval = FALSE, 
                                         risk.table = FALSE,
                                         legend = "right", 
                                         censor.shape = "|",
                                         censor.size = 4,
                                         palette = c("#FF6F61", "#008080"),
                                         ylab = "Proportion Surviving", 
                                         xlab = "Time (bins)",
                                         title = threshold_params)
  }
  
  survminer::arrange_ggsurvplots(splots, print = TRUE, ncol = 1, nrow = 1, 
                                 risk.table.height = 0.2)
  
  return (mget(c("comp_data", "prediction_threshold_search_granularity")))
  
}

#' Generates evaluation metrics, include time-dependent TPR and FPR rates as well as AUC
#'
#' @param train_result return data object from spect_train
#' @param prediction_times a vecotr of times to use for generating TPR and FPR data
#' @param plot_roc optional indicator to display the time-dependent ROC curves.
#' The TPR and FPR data will be returned regardless of the value of this parameter.
#' 
#' @return Evaluation metrics. Also plots the number of requested samples
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords visualization
#'
#' @export
#'
evaluate_model <- function(train_result, prediction_times, plot_roc=TRUE){
  
  flog.info("Collecting evaluation metric data...")
  event_probs <- matrix(NA, nrow(train_result$data.test), length(prediction_times))
  rownames(event_probs) <- 1:nrow(train_result$data.test)
  
  flog.debug("Collecting probabilities for %i prediction times", length(prediction_times))
  
  upper_bounds <- train_result$bounds$upper_bounds
  
  for(i in 1:length(prediction_times)) {
    #include all intervals where t_end <= t
    prediction_bucket <- length(upper_bounds[upper_bounds<=prediction_times[i]])    
    if(prediction_bucket==0) {
      event_probs[,i] <- 0 # prob of experiencing the event at time 0 is 0
    } else {
      event_probs[,i] <- train_result$data.test_predictions$abs_event_prob[train_result$data.test_predictions$interval == prediction_bucket]
    }
  }
  
  flog.debug("Capturing metrics from riskRegression...")
  
  # It is required to load the Surv() function separately from the survival package
  # in order to dynamically pass it to the Score function for calling. I'm sure
  # it's not best practice, but it's how this function operates.
  library(survival, include.only = "Surv")
  
  eval_formula <- stats::as.formula(paste("Surv(", train_result$survival_time_var,",", train_result$event_indicator_var, ") ~ 1"))
  
  # Pass the Score function a list of times and the probability of each individual experiencing 
  # the event before that time and it will spit out the time-dependent ROC/AUC data.
  # This only works on held-out test data where the ground truth is known.
  eval_metrics <- riskRegression::Score(list("ModelResults" = event_probs),
                                        formula=eval_formula,
                                        data=train_result$data.test,
                                        times=prediction_times,
                                        summary="ibs",
                                        plots = "ROC")
  
  if(plot_roc == TRUE){
    for(i in 1:length(prediction_times)) {
      riskRegression::plotROC(eval_metrics, times=prediction_times[i])
    }
  }
  
  return(mget(c("eval_metrics", "prediction_times")))
  
}

#' Generates predictions for each individual at each interval defined by the `train_result`
#' parameter. The interval-level predictions can be combined to generate surivival curves
#' for an individual.
#'
#' @param train_result - return data object from spect_train
#' @param new_data - New data set with the same covariates as the training data set.
#' 
#' @return predictions by the trained model on a new data set
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#'
#' @export
#'
spect_predict <- function(train_result, new_data){
  
  flog.info("Transforming %i rows of data", nrow(new_data))
  predictions <- create_person_period_data(new_data, train_result$bounds)
  
  # Important! The predict function returns the probability of class zero (X0), which is the
  # probability of no event in the given interval. Therefore, the interpretation of the calculated
  # values below is as follows:
  
  # cond_surv_pred: The conditional probability of observing no event in this interval, given that
  # the event was not observed prior to this interval.
  # abs_surv_pred: The absolute probability of observing no event up to and including the 
  # given interval
  # abs_event_prob: The probability of the event occurring sometime before the end of the current
  # interval. Note that this includes the current inverval.
  
  
  # Note: I do not love sprinkling if statements throughout the train function to 
  # handle single vs. ensemble modeling. Unfortunately, caret and caretstack 
  # return very similar, but not exactly the same objects, so we have to handle 
  # that here. In a future iteration, I'd like to find a more elegant (perhaps 
  # dictionary-based?) solution for this.
  
  flog.debug("Generating predictions...")
  if (length(train_result$base_learner_list) == 0){
    predictions$cond_surv_pred <- stats::predict(train_result$results, newdata=predictions, type="prob")[["X0"]]
  } else {
    predictions$cond_surv_pred <- stats::predict(train_result$results, newdata=predictions)[["X0"]]
  }
  
  flog.debug("Calculating probabilities...")  
  predictions$abs_surv_pred <- stats::ave(predictions$cond_surv_pred, predictions$individual_id, FUN=cumprod)
  predictions$abs_event_prob <- 1 - predictions$abs_surv_pred
  
  flog.debug("Calculating inverval boundaries for interpolation...")
  predictions$interval <- as.numeric(predictions$interval)
  predictions$next_interval <- predictions$interval + 1
  
  return (predictions)
  
}
