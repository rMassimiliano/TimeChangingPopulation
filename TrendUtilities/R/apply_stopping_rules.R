#' Early Futility Stopping Rule
#'
#' Drop experimental arms for futility if the posterior probability of efficacy becomes small.
#'
#' @details
#' The function drops experimental arms for futility if the  posterior probabilities of treatment efficacy
#'  dropps belows a treshold. The futility treshold depends on the number of observed outcomes \eqn{N_a} for each arm a,
#'
#'  \code{PrPTE[a]} \eqn{=p(\theta_a > \theta_0 | Data) \le f * (N_a/n'_E)^g}
#'
#'
#' @param i          Number of patients randomized.
#'
#' @param Active     A 2xA matrix which contains the activity status for each arm (1. column control).
#'                   The element \code{Active[1,a]} equals 1,0,2
#'                   if arm a is active, was dropped for futility or graduated for efficacy.
#'
#' @param PrPTE      \code{PrPTE} is a vector which contains, for each experimental arm,
#'                   the posterior probability of treatment efficacy
#'                   \code{PrPTE[a]} \eqn{=p(\theta_a > \theta_0 | Data)}
#'                   The vector can be generated with \code{\link{Posterior.with.control}}.
#'
#' @param stop.par   Vector of two elements, \code{stop.par[1]}
#'                   the treshold probability after the  maximum patient accrual is reached and
#'                   a non-negative shape parameter \code{stop.par[2]}.
#'
#' @param Accrual    Current number of observed outcomes for each arm.
#' @param n.a.max    Maximum sample size for each experimental arm.
#'
#'@details This function is implemented in the R package AddArmToTrial; refer to \url{http://bcb.dfci.harvard.edu/∼steffen/software.html}.
#' @export
apply_stopping_rules = function(i, Active, PrPTE, stop.par=c(f=0, gf=1), Accrual, n.a.max)
{
  ## futility check and activity
  id                = (PrPTE <= stop.par[1]*(Accrual[-1]/n.a.max) ^ stop.par[2]) & (Active[1,-1]==1)
  Active[,-1][ ,id] = c(0,i)


  return(Active)
} 
