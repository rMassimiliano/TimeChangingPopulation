\name{randomizedPlayTheWinner}
\alias{randomizedPlayTheWinner}
\title{Simulate a two-arm trial using a randomized play the winner design  as in Wei and Durham (1978)}
\usage{
randomizedPlayTheWinner(potential, beta = 1, alpha = 1,  M =1, ret_prob = FALSE)
}
\arguments{
\item{outcome}{A matrix with one potential outcome data for each arm. The first column is the enrolment time, the second a place holder for the arm the third and fourth potential outcome for the SOC and the experimental arm}
\item{M}{Randmization probabilities are updated after the enrolment of the M-th patient}
\item{alpha,beta}{these are design parameters as in Wei and Durham (1978)}
\item{ret_prob}
{if \code{TRUE} return the probability to be assigned to the experimental arm for each of the subject}
}
\value{
The function return a matrix with the same entrance as \code{potential}. The second column includes allocation of patients to arms according to a two-arm randomized play the winner design }
\description{This function implement a two-arm trial using using a randomized play the winner design  as in Wei and Durham (1978)} 


\references{
Wei, L. J. and Durham, S. (1978)
\emph{The randomized play-the-winner rule in medical trials.}
\emph{Journal of the American Statistical Association} 73, 840--843
}
