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
# TODO: Add calculations for probability of surviving to time x and probability of failing before time X
# TODO: Add common error handling:
# 1. Passing a "times" threshold beyond the final bounds
# 2. missing a package
# 3. others as users encounter them

# Required for plotSynData and other functionality
#' @import ggplot2
library(ggplot2)

# Modeling methods
library(caret)
library(caretEnsemble)
library(survival)
library(survminer)
library(riskRegression)

# For parallelization in Caret
library(parallel)
library(doParallel)

# Logging
library(futile.logger)

library(dplyr)

##############################################################################
# Autosurv fucntions: modified from https://github.com/ksuresh17/autoSurv

# Function that returns binary process data set ---------------------------
# This function was adapted from code from the Github account of Eric Polley
# https://github.com/ecpolley/SuperLearner_Old/blob/master/R/createDiscrete.R
createDiscreteDat <- function(timeInternal, eventInternal, dataX, delta.upper, cens="same") {
  #Assumes t0=0
  delta.lower <- c(0, delta.upper[-length(delta.upper)])
  n.delta <- length(delta.upper)
  IDInternal <- 1:nrow(dataX)
  dat_i <- cbind(IDInternal, timeInternal, eventInternal, dataX)
  interval <- rep(1:length(delta.upper), times=nrow(dataX))
  
  # Create a long data set that repeats each persons data n.delta (number of intervals) times
  long.dat <- dat_i[rep(seq_len(nrow(dat_i)), each = n.delta), ]
  
  N.delta <- rep(NA, nrow(long.dat))
  long.dat <- cbind(long.dat, delta.lower, delta.upper, N.delta, interval)
  
  # Treatment of censored observations
  if(cens=="same") {
    # Include censored observations in the interval in which they were censored
    long.dat$N.delta <- ifelse(long.dat$timeInternal > long.dat$delta.upper, 0,
                               ifelse(long.dat$eventInternal==1, ifelse(long.dat$timeInternal <= long.dat$delta.lower, NA, 1),
                                      ifelse(long.dat$timeInternal>long.dat$delta.lower, 0, NA)))
  } else if(cens=="prev") {
    # Include censored observations in the previous interval
    long.dat$N.delta <- ifelse(long.dat$timeInternal > long.dat$delta.upper, 0,
                               ifelse(long.dat$eventInternal==1, ifelse(long.dat$timeInternal <= long.dat$delta.lower, NA, 1),
                                      NA))
  } else if(cens=="half") {
    # Include censored observations in interval if they have survived at least half of that interval, otherwise include in previous interval
    long.dat$N.delta <- ifelse(long.dat$timeInternal > long.dat$delta.upper, 0,
                               ifelse(long.dat$eventInternal==1, ifelse(long.dat$timeInternal <= long.dat$delta.lower, NA, 1),
                                      ifelse(long.dat$timeInternal>=0.5*(long.dat$delta.upper+long.dat$delta.lower), 0, NA)))
  }
  
  m <- delta.upper[n.delta]
  # For event time corresponding to last interval include as an event in that interval
  long.dat$N.delta <- ifelse(long.dat$timeInternal == m & long.dat$delta.upper == m,
                             ifelse(is.na(long.dat$N.delta), 0, long.dat$N.delta), long.dat$N.delta)
  
  # Drops the intervals at which the individual no longer contributes
  long.dat <- long.dat[!is.na(long.dat$N.delta), ]
  # Sort the data
  long.dat <- long.dat[order(long.dat$IDInternal, long.dat$delta.lower),]
  # Set interval and outcome as a factor
  long.dat$interval <- as.factor(long.dat$interval)
  long.dat$N.delta <- as.factor(long.dat$N.delta)

  return(long.dat)
}

# Function that identifies the endpoints for the discrete intervals --------
createDiscreteIntervals <- function(time, event, bins) {
  #Number of intervals to create (minimum of specified bins and number of events)
  n.delta <- min(bins, length(unique(time[event==1])))
  #Get percentiles for each of the bins (from 0 to 1)
  probs.delta <- seq(from=0, to=1, length.out=(n.delta+1))
  #Get event times for each of the percentiles (0: first event time in dat set, 1: last event time in data set)
  delta.upper <- quantile(time[event==1], probs=probs.delta, names=FALSE)
  #Drop the first interval
  return(delta.upper[-1])
}

# Function to create person-time data set ---------------------------------
# dat: (data.frame) data set to convert (with columns named time_var, event_var)
# bins: (integer) number of discrete intervals
# w: (integer) prediction interval
genPTdat <- function(dat, time_var = time_var, event_var = event_Var, bins = 5, w, cens) {
  time.ind <- which(eval(time_var) == colnames(dat))
  event.ind <- which(eval(event_var) == colnames(dat))
  
  dat_i <- dat
  dat_i[dat_i$survTime > w, event.ind] <- 0 # ignore events after the prediction window
  dat_i[dat_i$survTime > w, time.ind] <- w # administratively censor at w
  
  # Identify the times t1, t2,... tJ (assumes t0=0)
  # Don't need to include prediction horizon as a final interval if there are no event times in that interval (prediction close to 0)
  delta.upper <- unique(createDiscreteIntervals(time = dat_i[,time.ind],
                                                event = dat_i[event.ind],
                                                bins = bins))
  delta.upper.new <- c(delta.upper[-length(delta.upper)],w)
  delta.upper <- delta.upper.new
  
  # Creates a new data set where each person has a row corresponding to the discrete time intervals
  # 0: in intervals where alive but does not have event
  # 1: in intervals where experience the event
  # delta.upper and delta.lower columns: are the same for everyone and are the discretized time intervals based on quantiles
  dat_i.X <- createDiscreteDat(timeInternal = dat_i[,time.ind],
                               eventInternal = dat_i[,event.ind],
                               dataX = dat_i[,-c(time.ind,event.ind)],
                               delta.upper = delta.upper,
                               cens=cens)
  
  return(dat_i.X)
}

# Function that formats data we will make predictions on ------------------
# Similar to "genPTdat" function but does not create outcome variable and creates
# a row for each interval for all subjects
# Note from Stephen Abrams: Added the upper_bound column to make graphing predictions much simpler later.
createPredictData <- function(dataX, delta.upper){
  # Need to create a prediction set with the relevant data for this window
  IDInternal <- 1:nrow(dataX)
  n <- nrow(dataX)
  dataX$IDInternal <- IDInternal
  n.delta <- length(delta.upper)
  interval <- rep(1:n.delta, times=nrow(dataX))
  upper_bound <- rep(delta.upper, times=nrow(dataX))
  long.data <- dataX[rep(seq_len(nrow(dataX)), each=n.delta), ]
  long.data <- cbind(long.data, interval, upper_bound)
  # Sort the data by ID and interval
  long.data <- long.data[order(long.data$IDInternal, long.data$interval),]
  long.data$interval <- factor(long.data$interval)
  return(long.data)
}
# End autosurv functions #
#################################################################################

#################################################################################
# Synthetic Data fucntions #

#' Generates a survival data set for synthetic streaming service subscription data. The survival 
#' event in this case is a cancellation of the subscription. It is given as a function of 
#' household income and average number of hours watched in the prior month. Users can adjust the
#' level of censoring and variance in the data with teh supplied parameters.
#'
#' @param sample_size (= 250) size of the sample population to generate
#' @param minimum_income (= 5000) minimum household income used to generate the distribution
#' @param median_income (= 50000) median household income used to generate the distribution
#' @param income_variance (= 10000) variance to use when generating the household income distribution
#' @param min_watchhours (= 0.0) minimum average number of hours watched used to generate the distribution
#' @param max_watchhours (= 6.0) minimum average number of hours watched used to generate the distribution
#' @param censor_percentage (= 0.0) percentage of population to artificially censor
#' @param min_censor_amount (= 0.0) Minimum number of months of censoring to apply to the censored population 
#' @param max_censor_amount (= 0.0) maximum number of months of censoring to apply to the censored population
#' @param study_time_in_months (= 48) observation horizon in months
#' @param perturbation_shift (= 0) defines a boundary for the amount to randomly perturb the formulaic result. Zero for no perturbation
#'
#' @return nothing
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @seealso \code{\link{createSynData}}
#' @keywords utilities
#'
#' @examples
#' createSynData(param)
#'
#' @export
createSynData <- function(sample_size = 250
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
  
  income_candidates <- rnorm(sample_size, median_income, income_variance)
  incomes <- ifelse(income_candidates < minimum_income, minimum_income, income_candidates)
  watchtimes <- runif(sample_size, min_watchhours, max_watchhours)
  
  # Apply the formula: 26 + watch_time^2 - (Income/10000)
  # This seems reasonable by inspection
  baseline_time_to_cancel <- 26 + watchtimes * watchtimes - (incomes/10000)
  
  # Perturb the result by a random number of months between +perturbation_shift and -perturbation_shift
  perturbed_baseline <- baseline_time_to_cancel + runif(sample_size, -perturbation_shift, perturbation_shift)
  
  # Randomly select %5 of the population for censoring at min_censor_amount to max_censor_amount months prior to the last observation.
  censor_selector <- ifelse(runif(sample_size, 0, 1) < censor_percentage, 1, 0)
  censor_amount <- runif(sample_size, min_censor_amount, max_censor_amount)
  censored_baseline <- perturbed_baseline - censor_selector * censor_amount
  
  # With censoring applied, this is now the "time_var"
  total_months <- ifelse(censored_baseline > study_time_in_months, study_time_in_months, censored_baseline)
  
  # If the event happened within the study window AND was not artificially censored, then consider the event to have happened, otherwise consider it censored. This gives us the event_Var for the person-period data set.
  event_within_study_time <- ifelse(perturbed_baseline <= 48, 1, 0)
  cancel_event_detected <- event_within_study_time * (1 - censor_selector)
  
  modeling_data <- data.frame(incomes, watchtimes, total_months, cancel_event_detected, baseline_time_to_cancel, perturbed_baseline)
  
  return(modeling_data)
  
}

#' Simple visualization of synthetic subscription data.
#'
#' @param data a data object generated by createSynData
#' 
#' @return None - prints synthetic data generated by createSynData
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @seealso \code{\link{plotSynData}}
#' @keywords utilities
#'
#' @examples
#' data <- createSynData()
#' plotSynData(data)
#'
#' @export
plotSynData <- function(data)
{
  # Plot incomes and watchtimes
  print(ggplot(data.frame(income = data$incomes), aes(x=income)) + geom_density())
  print(ggplot(data.frame(Watch_times = data$watchtimes), aes(x=Watch_times)) + geom_density())
  print(ggplot(data.frame(cancel_months = data$baseline_time_to_cancel), aes(x=cancel_months)) + geom_density())
  print(ggplot(data.frame(perturbed_cancel_months = data$perturbed_baseline), aes(x=perturbed_cancel_months)) + geom_density())
  print(ggplot(data.frame(final_cancel_months = data$total_months), aes(x=final_cancel_months)) + geom_density())
}

# End Synthetic Data fucntions #
#################################################################################

#################################################################################
# ensemble functions

# Purpose to extend the autosurv functionality to incorporate ensemble learning methods.
# Note - this will make use of some of the existing autosurv functions, but requires additional
# parameters, notably, a stack algorithm and a list of base learner algorithms, all of which should
# be able to perform binary classification.


#' Trains an ensemble model on the supplied survival data
#'
#' Generates a stacked ensemble model with the given caret models as base learners
#'  and metalearner.
#'
#'
#' @param test_prop (= 0.2) - proportion of the data set to reserve for testing
#' @param censor_type (= 'half') method used to determine censorship in a given bin - may be "half", "prev" or "same". see createDiscreteDat for usage.
#' @param bin_slices (= 10) Number of intervals to use for predictions.
#' @param method (= 'repeatedcv') caret parameter
#' @param resampling_number (= 10) for repeated cv 
#' @param kfold_repeats (= 3) number of folds
#' @param model_algorithm Ensemble meta-learner algorithm.
#' @param base_learner_list # base learner algorithms
#' @param metric (= 'Kappa') Is there a better metric? Maybe a custom one?
#' @param rng_seed (= 42) Optionally set the RNG seed for reproducibility
#' @param use_parallel (= TRUE) make use of Caret multicore training cluster
#' @param cores (= 0) Only relevant if use_parallel is set to TRUE. If zero, spect will attempt to make a good choice
#' @param modeling_data This data set must have one column for time, one for the indicator, and the rest are treated as covariates
#' @param event_indicator_var The name of the column containing the event indicator (zero or one)
#' @param survival_time_var The name of the column containing the time variable
#' @param obs_window The last time to use for generating person-period data. 
#' Any event occurring after this time will be administratively censored. In general, 
#' choosing a time at or near the end of the max observed time will include most events.
#' @param flog_level (= INFO) Logging verbosity - options are INFO, WARN, and DEBUG. Note, this can also be set in the calling function using the fultile.logger package.
#'
#' @return Vector of character strings representing the chosen set of colors, in RGB.
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @seealso \code{\link{spect_train}}
#' @keywords utilities
#'
#' @examples
#' spect_train(param)
#'
#' @export
#'
spect_train <- function(
  # Data preparation parameters
  test_prop = 0.2 # proportion of the data set to reserve for testing
  , censor_type = 'half' # method used to determine censorship in a given bin - may be "half", "prev" or "same". see createDiscreteDat for usage.
  , bin_slices = 10 # Number of intervals to use for predictions.
  # Ensemble-specific parameters. Defaults provided for all but algorithm choices.
  , method = 'repeatedcv' 
  , resampling_number = 10 
  , kfold_repeats = 3
  , model_algorithm # Primary algorithm used for classification
  , base_learner_list = list()# base learner algorithms. When this optional parameter is supplied, the model_algorithm is interpreted as a meta-learner for the ensemble stacking approach.
  , metric = 'Kappa'
  , rng_seed = 1101 # Optionally set the RNG
  , use_parallel = TRUE # make use of Caret multicore training cluster
  , cores = 0 # Only relevant if use_parallel is set to TRUE. If zero, detect will attempt to make a good choice
  , modeling_data # This must have one column for time, one for the indicator, and the rest are covariates
  , event_indicator_var
  , survival_time_var
  , obs_window
  )
{
  flog.info("Splitting test/train data...")
  testIndexes <- sample(x = 1:nrow(modeling_data), size = nrow(modeling_data) * test_prop)
  
  flog.debug("Creating test data...")
  data.test <- modeling_data[testIndexes, ]

  flog.debug("Creating training data...")
  data.train <- modeling_data[-testIndexes, ]        

  flog.debug("Split complete!")
  
  flog.info("Creating person-period data set...")
  data.person_period <- genPTdat(dat = data.train
                                          , time_var = survival_time_var
                                          , event_var = event_indicator_var
                                          , bins = bin_slices
                                          , w = obs_window
                                          , cens = censor_type)
  
  flog.info("Prepping data for ensemble modeling...")
  
  flog.debug("Stripping identifiers for modeling only covariates...")
  data.covs <- subset(data.person_period, select= -c(IDInternal, timeInternal, eventInternal, delta.lower, delta.upper))
  
  flog.debug("Factorizing the N.delta values (X0, X1) for classification modeling...")
  data.covs$N.delta <- make.names(data.covs$N.delta)
  
  flog.debug("Setting RNG seed to %i", rng_seed)
  set.seed(rng_seed)
  
  if(use_parallel){
    
    flog.info("Parallel processing specified - creating cluster for Caret modeling...")
    machine_cores <- detectCores()
    used_cores <- max(3 * machine_cores / 4, 1)
    
    flog.debug("Using %i of %i detected cores for parallel processing", used_cores, machine_cores)
    cl <- makePSOCKcluster(used_cores)
    registerDoParallel(cl)
  }
  
  flog.info("Creating Caret train control...")
  modelControl <- trainControl(method=method, number=resampling_number, repeats=kfold_repeats,
                               savePredictions=TRUE, classProbs=TRUE, allowParallel = use_parallel)
  if (length(base_learner_list) == 0){
    flog.info('Training %s using %s method with %i resamples and %i kfold repeats'
              , model_algorithm, method, resampling_number, kfold_repeats)
    print(data.covs$N.)
    results <- train(N.delta~., data=data.covs, method=model_algorithm, metric=metric, trControl=modelControl)
  }
  else
  {
    flog.debug("Pruning algorithm list...")
    algorithm_list <- base_learner_list[!is.na(base_learner_list)]
    
    flog.info("Training base models %s using %s method with %i resamples and %i kfold repeats"
              , paste(algorithm_list, collapse = ", "), method, resampling_number, kfold_repeats)
  
  
    base_models <- caretList(N.delta ~ ., data=data.covs, trControl=modelControl, methodList=algorithm_list)
    
    flog.info("Training the metalearner %s using %s metric...", model_algorithm, metric)
    results <- caretStack(base_models, method=model_algorithm, metric=metric, trControl=modelControl)
  }
  
  if(use_parallel){
    flog.info("Stopping parallel computing cluster...")
    stopCluster(cl)
    rm(cl)
  }

  flog.info("Calculating probabilities on %i test individuals...", nrow(data.test))

  upper_bin_bounds <- unique(data.person_period$delta.upper)
  data.test_predictions <- createPredictData(dataX = data.test, delta.upper=upper_bin_bounds)

  flog.debug("calculating individual and cumulative probabilities by individual + interval...")  
  
  # Important! The predict function returns the probability of class zero (X0), which is the
  # probability of no event in the given interval. Therefore, the interpretation of the calculated
  # values below is as follows:

  # cond_surv_pred: The conditional probability of observing no event in this interval, given that
  # the event was not observed prior to this interval.
  # abs_surv_pred: The absolute probability of observing no event up to and including the 
  # given interval
  # abs_event_prob: The probability of the event occurring sometime before the end of the current
  # interval. Note that this includes the current inverval.
  data.test_predictions$cond_surv_pred <- predict(results, newdata=data.test_predictions, type="prob")
  
  # Note: I do not love sprinkling if statements throughout the train function to 
  # handle single vs. ensemble modeling. Unfortunately, caret and caretstack 
  # return very similar, but not exactly the same objects, so we have to handle 
  # that here. In a future iteration, I'd like to find a more elegant solution.
  if (length(base_learner_list) == 0){
    data.test_predictions$cond_surv_pred <- data.test_predictions$cond_surv_pred[['X0']]
  }
  
  data.test_predictions$abs_surv_pred <- ave(data.test_predictions$cond_surv_pred, data.test_predictions$IDInternal, FUN=cumprod)
  data.test_predictions$abs_event_prob <- 1 - data.test_predictions$abs_surv_pred
    
  flog.debug("Calculating inverval boundaries for interpolation...")
  data.test_predictions$interval <- as.numeric(data.test_predictions$interval)
  data.test_predictions$next_interval <- data.test_predictions$interval + 1
  
  return (mget(c("data.test", "data.train", "data.person_period", "data.covs", "results"
                 , "data.test_predictions", "event_indicator_var", "survival_time_var"
                 , "bin_slices", "upper_bin_bounds", "obs_window")))
  
}

#' Plots a sample of individual survival curves from the test data set.
#'
#'
#' @param train_result - return data object from spectTrain
#' @param samples_to_plot - number of samples to plot
#' 
#' @return None - plots the number of requested samples
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @seealso \code{\link{plot_samples}}
#' @keywords utilities
#'
#' @examples
#' plot_samples(param)
#'
#' @export
#'
plot_samples <- function(train_result, samples_to_plot){
  
  flog.info("Sampling %i test individuals for plots...", samples_to_plot)
  event_var <- train_result$event_indicator_var
  time_var <- train_result$survival_time_var
  pred_times <- train_result$obs_window
  bin_slices <- train_result$bin_slices
  
  data <- train_result$data.test_predictions
  sample <- round(runif(samples_to_plot, 1, max(data["IDInternal"])))

  flog.info("Found samples, plotting...")  
  for (individual in sample){
    
    individual_predictions <- data[data$IDInternal == individual,]
    x_max <- max(individual_predictions$upper_bound)
    
    flog.debug("Collecting actuals for the plot label...")
    actual_absolute_time <- round(unique(individual_predictions[, time_var]))
    actual_interval <- round(actual_absolute_time / (pred_times / bin_slices))
    event_detected <- unique(individual_predictions[, event_var])
    
    title_value = paste("Survival Curve: Individual", individual)
    subtitle_value <- paste("Actual Time:", actual_absolute_time, "Actual Interval:", actual_interval,"Event?", event_detected)

    flog.debug("Plotting individual %i", individual)    
    plot(0,0,xlim = c(0,x_max),ylim = c(0,1),type = "n", xlab="", ylab="")
    lines(individual_predictions$upper_bound, individual_predictions$cond_surv_pred, type='b', col='#FFA500', lwd="4")
    lines(individual_predictions$upper_bound, individual_predictions$abs_surv_pred, type='b', col='#1E90FF', lwd="4")
    
    title(main = title_value, sub = subtitle_value,
          xlab = "Time Interval", ylab = "Survival Probability",
          cex.main = 2,   font.main= 4, col.main= "#808080",
          cex.sub = 1, font.sub = 4, col.sub = "#800000",
          col.lab = "#808080", 
    )
  }
}

#' Plots a series of population Kaplan-Meier curves for different thresholds for both the test 
#' predictions and the ground truth
#'
#'
#' @param train_result - return data object from spect_train
#' @param prediction_threshold_search_granularity (= 0.05) # A number between zero and one which defines the granularity of searching for cumulative probability thresholds. 
#' For instance, search a value of 0.05 will search 19 thresholds (0.05, 0.10, ..., 0.95)
#' 
#' @return None - plots the number of requested samples
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @seealso \code{\link{plot_KM}}
#' @keywords utilities
#'
#' @examples
#' plot_KM(param)
#'
#' @export
#'
plot_KM <- function(train_result, prediction_threshold_search_granularity = 0.05){
  
  survival_time_var <- train_result$survival_time_var
  event_indicator_var <- train_result$event_indicator_var

  start <- prediction_threshold_search_granularity
  end <- 1 - prediction_threshold_search_granularity
  increments <- (1 / prediction_threshold_search_granularity) - 1
  
  # Store the prediction vs. ground truth data for each threshold
  comp_data <- list()
  
  flog.debug("Checking thresholds...")
  for (prediction_threshold in seq(start, end, length.out = increments)){
    
    flog.info("collecting the minimum interval which crosses the prediction threshold: %f", prediction_threshold)
    
    # Note: The minimum interval is just a choice. We could have arbitrarily chosen the middle of the interval or 
    # the lower bound, or even decided dynamically based on the distribution.
    prediction_interval <- train_result$data.test_predictions[train_result$data.test_predictions$abs_surv_pred < prediction_threshold,]
    prediction_interval$interval <- as.numeric(prediction_interval$interval)
    threshold_bound <- do.call(rbind, unname(by(prediction_interval, prediction_interval$IDInternal, function(x) x[x$interval == min(x$interval),])))
    
    flog.info("Found %i individuals...", nrow(threshold_bound))
    
    # Join back to the original data to pick up the lower bound of that interval.
    # This will be used to create a linear interpolation.
    # TODO: This actually removes anything for which the event occurred in the first
    # interval. It won't happen often for large bin sizes, but should probably be fixed.
    under_threshold <- inner_join(train_result$data.test_predictions
                                  , threshold_bound
                                  , by=c("IDInternal", "next_interval" = "interval"))
    
    flog.debug("Interpolating to find exact threshold crossing time...")
    # Note: This is the "true" prediction time (measured in bins, which can be translated
    # to the original time unit later)
    under_threshold$true_interval <- with( 
      under_threshold, interval + ((prediction_threshold - abs_surv_pred.y) 
                                   / (abs_surv_pred.x - abs_surv_pred.y))) 
    
    flog.debug("Finding individuals which never cross the threshold...")
    
    removals <- unique(under_threshold$IDInternal)
    over_threshold <- train_result$data.test_predictions %>% filter(!(IDInternal %in% removals))
    
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
      select(predicted_crossing_time = interval, event_indicator_var, survival_time_var)
    
    threshold_hits <- rbind(threshold_hits_under, threshold_hits_over)
    
    flog.debug("Interpolating ground truth crossing time...")
    # Note: convert the real survival time to a fractional interval for the "Ground Truth" crossing time.
    threshold_hits$GT_crossing_time <- (threshold_hits[, survival_time_var] 
                                        * (1.0 * train_result$bin_slices)
                                        / (1.0 * train_result$obs_window) )
    
    predicted_elements <- select(threshold_hits, c(predicted_crossing_time, event_indicator_var))
    predicted_elements$source <- "Predictions"
    colnames(predicted_elements)[1] ="TimeSlice"
    
    ground_truth_elements <- select(threshold_hits, c(GT_crossing_time, event_indicator_var))
    ground_truth_elements$source <- "Ground Truth"
    colnames(ground_truth_elements)[1] ="TimeSlice"
    
    gt_comp_data <- rbind(predicted_elements, ground_truth_elements)
    colnames(gt_comp_data)[2] = "EventIndicator"
    
    gt_comp_data$TimeSlice <- as.double(gt_comp_data$TimeSlice)
    
    flog.info("Generate P-value...")    
    
    KM_diff <- survdiff(Surv(TimeSlice, EventIndicator) ~ source, data=gt_comp_data)
    pval<-pchisq(KM_diff$chisq, length(KM_diff$n)-1, lower.tail = FALSE)
    
    display_pval <- format(round(pval, 2), nsmall = 2)
    display_prediction_threshold <- format(round(prediction_threshold, 2), nsmall = 2)
    flog.info("Pval is %s for threhold %s", display_pval, display_prediction_threshold)
    
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

    flog.debug(paste("Preparing ground truth fit data for plot...iteration", as.character(i)))
    threshold_data <- comp_data[i,]
    
    fmt_pval <- format(round(threshold_data$pval, 2), nsmall = 2)
    fmt_prediction_threshold <- format(round(threshold_data$prediction_threshold, 2), nsmall = 2)

    threshold_params <- paste("Prediction Threhsold: ", fmt_prediction_threshold, " Pval: ",  fmt_pval)
    flog.debug(paste("Using the following parameters - ", threshold_params))
    
    gt <- threshold_data$gt_comp_data
    
    KM_fit <- survminer::surv_fit(Surv(TimeSlice, EventIndicator) ~ source, data=gt)
    
    splots[[i]] <- ggsurvplot(KM_fit, conf.int = TRUE, pval = FALSE, risk.table = FALSE, 
                              legend = "right", censor.shape = "|", censor.size = 4,
                              palette = c("#FF6F61", "#008080"), ylab = "Proportion Surviving", 
                              xlab = "Time (bins)", title = threshold_params)
  }
  
  # Arrange multiple ggsurvplots and print the output
  arrange_ggsurvplots(splots, print = TRUE,
                      ncol = 1, nrow = 1, risk.table.height = 0.2)
  
  # TODO: It would be good to add an option to generate this data without doing the plot.
  return (mget(c("comp_data", "prediction_threshold_search_granularity")))
  
}

#' Generates evaluation metrics, include time-dependent TPR and FPR rates as well as AUC
#'
#' @param train_result - return data object from spect_train
#' @param prediction_times - a vecotr of times to use 
#' @param plot_roc (=TRUE) - indicates whether to display the time-dependent ROC curves.
#' The TPR and FPR data will be returned regardless of the value of this parameter.
#' 
#' @return None - plots the number of requested samples
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords utilities
#'
#' @examples
#' evaluate_model(param)
#'
#' @export
#'
evaluate_model <- function(train_result, prediction_times, plot_roc=TRUE){
  
  flog.info("Collecting evaluation metric data...")
  event_probs <- matrix(NA, nrow(train_result$data.test), length(prediction_times))
  rownames(event_probs) <- 1:nrow(train_result$data.test)
  
  flog.debug("Collecting probabilities for %i prediction times", length(prediction_times))
  
  for(i in 1:length(prediction_times)) {
    #include all intervals where t_end <= t
    prediction_bucket <- length(train_result$upper_bin_bounds[train_result$upper_bin_bounds<=prediction_times[i]])
    
    if(prediction_bucket==0) {
      event_probs[,i] <- 0 # prob of experiencing the event at time 0 is 0
    } else {
      event_probs[,i] <- train_result$data.test_predictions$abs_event_prob[train_result$data.test_predictions$interval == prediction_bucket]
    }
  }
  
  flog.debug("Capturing metrics from riskRegression...")
  
  eval_formula <- as.formula(paste("Surv(", train_result$survival_time_var,",", train_result$event_indicator_var, ") ~ 1"))
  
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

#' Plots a series of population Kaplan-Meier curves for different thresholds for both the test 
#' predictions and the ground truth
#'
#'
#' @param train_result - return data object from spect_train
#' @param newdata - New data set with the same covariates as the training data set.
#' 
#' @return predictions by the trained model on a new data set
#'
#' @author Stephen Abrams, \email{stephen.abrams@@louisville.edu}
#' @keywords utilities
#'
#' @examples
#' spect_predict(param)
#'
#' @export
#'
spect_predict <- function(train_result, new_data){

  flog.info("Transforming %i rows of data", length(new_data))
  predictions <- createPredictData(dataX = new_data, delta.upper=train_result$upper_bin_bounds)
  
  # Important! See notes in the spect_train fucntion for details about what exactly the 
  # predictions and calculations here mean.

  flog.debug("Generating predictions...")
  predictions$cond_surv_pred <- predict(train_result$results, newdata=predictions, type="prob")

  flog.debug("Calculating probabilities...")  
  predictions$abs_surv_pred <- ave(predictions$cond_surv_pred, predictions$IDInternal, FUN=cumprod)
  predictions$abs_event_prob <- 1 - predictions$abs_surv_pred

  return (predictions)
    
  }


# end ensemble functions
#################################################################################
