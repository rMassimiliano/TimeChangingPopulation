\name{doublyAdaptiveBiasedCoin}
\alias{doublyAdaptiveBiasedCoin}
\title{Simulate a two-arm trial using a double adaptive biased coin design, e.g, Duan and Hu (2009)}
\usage{
doublyAdaptiveBiasedCoin(potential, gamma, M =1,ret_prob =FALSE)
}
\arguments{

 \item{outcome}
 {A matrix with one potential outcome data for each arm. The first column is the enrolment time, the second a place holder for the arm the third and fourth potential outcome for the SOC and the experimental arm}
\item{M}{Randomization probabilities are updated after the enrolment of the M-th patient}
\item{gamma}{this is a design parameter as in Duan and Hu (2009)} 
\item{ret_prob}
{if \code{TRUE} return the probability to be assigned to the experimental arm for each of the subject}
}
\value{
The function return a matrix with the same entrance as \code{potential}. The second column includes allocation of patients to arms according to a two-arm BAR design }
\description{This function implement a two-arm trial using a double adaptive biased coin design, e.g, Duan and Hu (2009)}

\references{
\Duan, L. and Hu, F. (2009).
\emph{Doubly adaptive biased coin designs with heterogeneous responses.}
\emph{Journal of statistical planning and inference} 139, 3220--3230.
}
