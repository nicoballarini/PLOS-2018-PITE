% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_scores.R
\name{score.null}
\alias{score.null}
\title{Fit linear model to create a score using lm with no interactions}
\usage{
score.null(dataset, input, parameters = NULL, verbose = TRUE,
  alpha = 0.05)
}
\arguments{
\item{dataset}{A data.frame created with \code{OneData}}

\item{input}{A list of objects created with \code{MakeInput}}

\item{parameters}{A \code{MakeParameters} object}

\item{verbose}{logical. whether to print a message with results.}

\item{alpha}{Significance level for confidence intervals (target is
miscoverage alpha/2 in each tail)}
}
\description{
Performs lasso in a dataset with main effects but no interactions. calculates
PITE and CI using the Average treatment effect (ATE)
}
\examples{
input <- MakeInput()
parameters <- MakeParameters(input = input)
dataset <- OneData(input,parameters)
score.null(dataset, input, parameters)

}

