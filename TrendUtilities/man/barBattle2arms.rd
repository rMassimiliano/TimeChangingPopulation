\name{barBattle2arms}
\alias{barBattle2arms}
\title{Simulate a two-arm trial using a Bayesian adaptive design. Similar to Zhou et al (2017) we let the allocation probabilities be proportional to the posterior means computed using the accrued data.} 
	
\usage{
barBattle2arms(potential, alpha = 1.0, beta = 1.0, M = 1, minProb = 0.1, ret_prob = FALSE)
}
\arguments{

 \item{outcome}
 {A matrix with one potential outcome data for each arm. The first column is the enrolment time, the second a place holder for the arm the third and fourth potential outcome for the SOC and the experimental arm}

\item{alpha, beta}{Hyperparameters of the prior Beta(alpha,beta) used to compute the posterior mean. The same parameters are used for the two arms}
\item{M}{Randomization probabilities are updated after the enrolment of the M-th patient}
\item{minProb}{The minimal allowed probability in treatments and control arm}
\item{ret_prob}
{if \code{TRUE} return the probability to be assigned to the experimental arm for each of the subject}
}
\value{
The function return a matrix with the same entrance as \code{potential}. The second column includes allocation of patients to arms according to a two-arm BAR design }
\description{
	Simulate a two-arm trial using a Bayesian adaptive design. Similar to Zhou et al (2017) we let the allocation probabilities be proportional to the posterior means computed using the accrued data.
}

\references{
Zhou, X., Liu, S., Kim, E. S., Herbst, R. S., and Lee, J. J. (2008)
\emph{Bayesian adaptive design for targeted therapy development in lung cancer|a step toward personalized medicine.} 
\emph{Clinical Trials} 5, 181--193.
}
