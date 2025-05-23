% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{acc_sim}
\alias{acc_sim}
\title{Utility function to generate accuracy metrics, for use with \code{\link[=estimate_accuracy]{estimate_accuracy()}}}
\usage{
acc_sim(
  n,
  method,
  p,
  dat,
  model,
  eta,
  nsample,
  outcome,
  power,
  effect_size,
  powersims,
  alpha,
  split,
  predictfn,
  replacement,
  ...
)
}
\arguments{
\item{n}{An integer giving the desired sample size for which the target function is to be calculated.}

\item{method}{An optional string stating the distribution from which data is to be generated. Default is i.i.d. uniform sampling. Currently also supports "Class Imbalance". Can also take a function outputting a vector of probabilities if the user wishes to specify a custom distribution.}

\item{p}{If method is 'Class Imbalance', gives the degree of weight placed on the positive class.}

\item{dat}{A rectangular \code{data.frame} or matrix-like object giving the full data from which samples are to be drawn. If left unspecified, \code{\link[=gendata]{gendata()}} is called to produce synthetic data with an appropriate structure.}

\item{model}{A function giving the model to be estimated}

\item{eta}{A real number between 0 and 1 giving the probability of misclassification error in the training data.}

\item{nsample}{A positive integer giving the number of samples to be generated for each value of $n$. Larger values give more accurate results.}

\item{outcome}{A string giving the name of the outcome variable.}

\item{power}{A logical indicating whether experimental power based on the predictions should also be reported}

\item{effect_size}{If \code{power} is \code{TRUE}, a real number indicating the scaled effect size the user would like to be able to detect.}

\item{powersims}{If \code{power} is \code{TRUE}, an integer indicating the number of simulations to be conducted at each step to calculate power.}

\item{alpha}{If \code{power} is \code{TRUE}, a real number between 0 and 1 indicating the probability of Type I error to be used for hypothesis testing. Default is 0.05.}

\item{split}{A logical indicating whether the data was passed as a single data frame or separately.}

\item{predictfn}{An optional user-defined function giving a custom predict method. If also using a user-defined model, the \code{model} should output an object of class \code{"svrclass"} to avoid errors.}

\item{replacement}{A logical flag indicating whether sampling should be performed with replacement.}

\item{...}{Additional model parameters to be specified by the user.}
}
\value{
A data frame giving performance metrics for the specified sample size.
}
\description{
Utility function to generate accuracy metrics, for use with \code{\link[=estimate_accuracy]{estimate_accuracy()}}
}
