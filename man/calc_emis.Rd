% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forward.R
\name{calc_emis}
\alias{calc_emis}
\title{Calculating the probability of occurance of particular values of covariates at a node given the value of target.}
\usage{
calc_emis(state, obsv, probs, pars)
}
\arguments{
\item{state}{character value of state variable at a particular node.}

\item{obsv}{character vector of values of covariates at that node.}

\item{probs}{emission probability distribution of the covariates in TAN structure.}

\item{pars}{integer vector denoting the parents of the nodes(other than root) in the TAN structure.}
}
\value{
probability of occurance of particular values of covariates at a node given the value of target.
}
\description{
Calculating the probability of occurance of particular values of covariates at a node given the value of target.
}
