% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spect.R
\name{spect_predict}
\alias{spect_predict}
\title{Generates predictions for each individual at each interval defined by the `train_result`
parameter. The interval-level predictions can be combined to generate surivival curves
for an individual.}
\usage{
spect_predict(train_result, new_data)
}
\arguments{
\item{train_result}{- return data object from spect_train}

\item{new_data}{- New data set with the same covariates as the training data set.}
}
\value{
predictions by the trained model on a new data set
}
\description{
Generates predictions for each individual at each interval defined by the `train_result`
parameter. The interval-level predictions can be combined to generate surivival curves
for an individual.
}
\author{
Stephen Abrams, \email{stephen.abrams@louisville.edu}
}
