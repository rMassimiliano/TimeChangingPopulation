% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apply_stopping_rules.R
\name{apply_stopping_rules}
\alias{apply_stopping_rules}
\title{Early Futility Stopping Rule}
\usage{
apply_stopping_rules(i, Active, PrPTE, stop.par = c(f = 0, gf = 1),
  Accrual, n.a.max)
}
\arguments{
\item{i}{Number of patients randomized.}

\item{Active}{A 2xA matrix which contains the activity status for each arm (1. column control).
The element \code{Active[1,a]} equals 1,0,2
if arm a is active, was dropped for futility or graduated for efficacy.}

\item{PrPTE}{\code{PrPTE} is a vector which contains, for each experimental arm,
the posterior probability of treatment efficacy
\code{PrPTE[a]} \eqn{=p(\theta_a > \theta_0 | Data)}
The vector can be generated with \code{\link{Posterior.with.control}}.}

\item{stop.par}{Vector of two elements, \code{stop.par[1]}
the treshold probability after the  maximum patient accrual is reached and
a non-negative shape parameter \code{stop.par[2]}.}

\item{Accrual}{Current number of observed outcomes for each arm.}

\item{n.a.max}{Maximum sample size for each experimental arm.}
}
\description{
Drop experimental arms for futility if the posterior probability of efficacy becomes small.
}
\details{
The function drops experimental arms for futility if the  posterior probabilities of treatment efficacy
 dropps belows a treshold. The futility treshold depends on the number of observed outcomes \eqn{N_a} for each arm a,

 \code{PrPTE[a]} \eqn{=p(\theta_a > \theta_0 | Data) \le f * (N_a/n'_E)^g}

This function is implemented in the R package AddArmToTrial; refer to \url{http://bcb.dfci.harvard.edu/∼steffen/software.html}.
}
