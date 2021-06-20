
#' Generate random numbers for a Bayesian Adaptive treatment assignment
#'
#' Generate a random vector of \code{nr.a} treatment assignment variables for the Bayesian adaptive design.
#'
#' The function generates integers between 1 and A, where A is the number of arms in the trial.
#' Each random integer \eqn{ A_i= \{1,\ldots, A\}} denotes the assignment of patient i to
#' the control (\eqn{A_i=1}) or experimental arm  (\eqn{ A_i= \{2,\ldots, A\}}).
#' The random numbers are drawn under a Bayesian outcome-adaptive assignment model.
#'
#'
#' @param Accrual    The number of patients randomized to each treatment arm.
#' @param k.c        Current number of groups of added arms (initial group counts too).
#'
#' @param batch      Group membership for each treatment arm (for control put 1).
#' @param m.k        Mean parameter of the Gompertz function for each treatment group,
#' @param n.k        Extension of the overall sample size of the trial after each group of added arms.
#'
#' @param PrPTE      \code{PrPTE} is a vector which contains, for each experimental arm,
#'                   the posterior probability of treatment efficacy
#'                   \code{PrPTE[a]} \eqn{=p(\theta_a > \theta_0 | Data)}.
#'                   The vector can be generated with \code{\link{Posterior.with.control}}.
#'
#' @param rand       A vector three tuning parameter for adaptive randomization \code{rand =c(a,b,c)}.
#'                   The first and second element correspond to the maximum
#'                   and the shap parameter of the power function \eqn{h(x) = a * (x/n)^b}
#'                   The last parameter controls the randomization to the control.
#'
#' @param q.vec      The 3 common parameter \eqn{ (r_0, r_1,r_3) } of the Gompertz function
#'                   \eqn{q_k(x) = r_0 + r_1 \exp(-\exp( -r_2 [x-m_k] ))}.
#'
#' @param n.a.max    Vector of maximum sample sizes for each experimental arm.
#' @param Active     Activity status for each experimental arm at the i-th arrival,
#'                   i.e. \code{Active[k]=1} or \code{TRUE} if  k-th arm is active.
#'
#' @param nr.a       How many random numbers should be generated
#'
#' @return
#' Returns a set of \code{nr.a} treatment assignment numbers.
#'@keywords internal
#'@details This function include minor changes  for the one in the AddArmToTrial; refer to \url{http://bcb.dfci.harvard.edu/∼steffen/software.html}.
#' @export
randomize_BAR_add = function(Accrual, k.c, batch, m.k, n.k, PrPTE, rand, q.vec, n.a.max, Active=rep(1,length(PrPTE)), nr.a=1){

  Accrual.to.batch =  sapply(1:k.c, function(k)  sum(Accrual[which(batch[-1]==k)+1]) )
  q.batch             =  q.vec[1] + q.vec[2] * exp(-exp(q.vec[3]*(Accrual.to.batch-m.k)))
  h.batch             =  rand[1] * (Accrual.to.batch/n.k)^rand[2]
  id.k                   =  Accrual.to.batch>=n.k
  id.n                   =  Accrual[-1] < n.a.max
  h.batch[id.k]     =  rand[1]

  h.a                    = h.batch[batch[-1]]
  q.a                    = q.batch[batch[-1]]

  # randomization probabilities for exparimetal arms
  id.active           = Active==1
  Pr.a                  = PrPTE^h.a * q.a * id.n * id.active

  # randomization probability for control arm
  P0 = exp(rand[3] * (max(Accrual[-1])-Accrual[1])) / sum(id.n * id.active)

  # combine all randomization probabilities
  Pr = c(P0, Pr.a/sum(Pr.a) )
  Pr = Pr/sum(Pr)

  # avoide extreme unbalancing
  id.u = (Pr[-1] < 1/(3* sum(id.active & id.n) )) & id.active & id.n
  Pr[-1][id.u]= 1/(3* sum(id.active & id.n))

  return(sample.int(length(Pr), size=1, prob=Pr))}
 
