% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_Conf-Intervals.R
\name{confidence.intervals.lassoan.test}
\alias{confidence.intervals.lassoan.test}
\title{Confidence Intervals for Lasso}
\usage{
confidence.intervals.lassoan.test(dataset, input, parameters, lasso.results,
  ndraw = 8000, burnin = 2000, alpha = 0.05)
}
\arguments{
\item{dataset}{A data.frame created with \code{OneData}}

\item{input}{A list of objects created with \code{MakeInput}}

\item{parameters}{A \code{MakeParameters} object}

\item{lasso.results}{output from the score.lasso.added.noise}

\item{ndraw}{Number of samples of optimization variables to sample.}

\item{burnin}{How many samples of optimization variable to discard.}

\item{alpha}{Significance level for confidence intervals (target is
miscoverage alpha/2 in each tail)}

\item{verbose}{logical. whether to print a message with results.}

\item{alpha}{Significance level for confidence intervals (target is
miscoverage alpha/2 in each tail)}
}
\description{
Calculates the confidence intervals for the score given a Lasso estimation
}
\examples{
input <- MakeInput()
parameters <- MakeParameters(input = input)
dataset <-OneData(input,parameters)
score.lasso(dataset, input)

}

