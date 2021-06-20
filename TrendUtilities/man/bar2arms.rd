\name{bar2arms}
\alias{bar2arms}
\title{Simulate a two-arm trial using a Bayesian adaptive design as in Thall (2015)}
\usage{
bar2arms(potential, c = 1.0,  M = 1, ret_prob = FALSE)
}
\arguments{

 \item{outcome}
 {A matrix with one potential outcome data for each arm. The first column is the enrolment time, the second a place holder for the arm the third and fourth potential outcome for the SOC and the experimental arm}
\item{M}{Randomization probabilities are updated after the enrolment of the M-th patient}
\item{c}{this is a design parameter as in Thall (2015)}
\item{ret_prob}
{if \code{TRUE} return the probability to be assigned to the experimental arm for each of the subject}
}
\value{
The function return a matrix with the same entrance as \code{potential}. The second column includes allocation of patients to arms according to a two-arm BAR design }
\description{This function implement a two-arm trial using a Bayesian adaptive design as in Thall (2015) . 
}

\references{
Thall, P., Fox, P., and Wathen, J. (2015)
\emph{Statistical controversies in clinical research: scientic and ethical problems with adaptive randomization in comparative clinical trials.}
\emph{Annals of Oncology} 26, 1621--1628.

}
