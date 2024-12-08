% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{weight}
\alias{weight}
\title{Indicator for the doubly robust ICE estimator}
\usage{
weight(treat_model = list())
}
\arguments{
\item{treat_model}{a list of formulas specifying the treatment model for the corresponding treatment variable. The length of list must match
the number of treatment variables.}
}
\value{
treatment model specifications and treatment variable names.
}
\description{
This function identifies the doubly robust ICE estimator. The treatment models could be specified by
\code{treat_model}.
}
