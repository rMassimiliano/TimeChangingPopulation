#' Simulate a Trial using Bayesian response-adaptive (BAR) randomization
#'
#' Generate a multi-arm trial under BAR whit the possibility of adding arms to the trial.
#'
#' @param Potential.Outcome
#' A matrix with one potential outcome data for each arm.
#' The data can be generated with \code{\link{Generat.Outcome}}.
#'
#' @param added.at
#' A vector of adding times for each arm.
#' For all initial arms use \code{0}.
#'
#' @param Outcome.delay
#' Expected time between enrolments and avalability of the outcome.
#' I.e. if the endpoint is respons to treatment at 6 month, \code{Outcome.delay=6}.
#'
#' @param n.a.max    Vector of maximum sample sizes for each experimental arm.
#'
#' @param rand       A vector \code{rand =c(H, gamma, b)} of randomization three tuning parameter.
#'                   The first and second element correspond to the maximum
#'                   and the shap parameter of the power function \eqn{h(x) = H * (x/n)^gamma}
#'                   The last parameter controls the randomization to the control.
#'
#' @param q.vec      The 3 common parameter \eqn{ (r_0, r_1,r_3) } of the Gompertz function
#'                   \eqn{q_k(x) = r_0 + r_1 \exp(-\exp( -r_2 [x-m_k] ))}.
#'
#' @param m.k                Mean parameter of the Gompertz function for each batch
#' @param n.k                Extension of the overall sample size for each added group
#'
#' @param stop.par           Vector of stopping parameters with default c(f=0, g=1).
#' @param do.futi.stop       If \code{TRUE} do sequantial futility analysis.
#'
#' @return
#' The function returns a list \code{list(TrialData = Potential.Outcome, Outcome=Outcome, Active=Active)} of three elements.
#'
#' \code{TrialData} is identical to the impute matrix \code{Potential.Outcome}
#' except that the 2nd  column now contains treatment assignment variables, i.e.
#' \code{Potential.Outcome[10,2]} equals one if the 10th patients was Accrual to the control.
#' \code{Potential.Outcome[10,2]} equals 2 if the 10th patients was Accrual to the arm 1.
#'
#' \code{Outcome} is an array of outcome data,
#' where \code{Outcome[,,k]} corresponds to the outcome data in the kth group of added arms.
#'
#' \code{Active} corresponds to the activity status of each arm.
#'
#'@details This function includes minor changes form the one implemented in the R package AddArmToTrial; refer to \url{http://bcb.dfci.harvard.edu/∼steffen/software.html}.
#' @export
simulate_BAR_add_trial =
  function(Potential.Outcome, added.at=rep(1, ncol(Potential.Outcome)), Outcome.delay=0, n.a.max,
           rand, m.k, n.k, q.vec, stop.par=c(f=0, gf=1), do.futi.stop=TRUE){

    N.A              = dim(Potential.Outcome) - c(0,2)
    M.k              = unique(added.at)
    A.k              = sapply(M.k, function(M.k) sum(added.at[-1]==M.k))
    K                 = length(M.k)
    batch          = sapply(added.at, function(a) which(M.k==a))
    Outcome     = array(0L, dim=c(2, N.A[2], K) )
    Active         = rbind(rep(1L,N.A[2]), rep(0L,N.A[2]))
    Accrual       = rep(0, N.A[2])
    n.last          = 0

##############  active accrual period #####################
for(i in 1:N.A[1])
{
      ## update data during the active accrual period
      id              = added.at<=i
      k.c            = sum(M.k<=i)
      t.i              = Potential.Outcome[i,1]
      Update.Outcome = update_data(i, t.i, Potential.Outcome, Outcome.delay, Outcome, n.last, batch, N.target.0=N.A[1], M.k, during.accrual=TRUE)
      Outcome             = Update.Outcome$Outcome
      n.last                   = Update.Outcome$n.last
      Outcome.c           = sapply(1:sum(id), function(a) Outcome[,a, batch[a]] )
      PrPTE                   = posterior_with_control(Outcome=Outcome.c)

      ## futility testing
      if(do.futi.stop)
        Active[,id]= apply_stopping_rules(i=i, Active=Active[,id], PrPTE=PrPTE, stop.par=stop.par, Accrual=Accrual[id], n.a.max=n.a.max[id[-1]])

      ## randomize if still active patients
      if( any(Active[1,id][-1]==1 & Accrual[id][-1] <n.a.max[id[-1]])  ){
        Potential.Outcome[i,2] = a.i = randomize_BAR_add(Accrual=Accrual[id], k.c=k.c, batch=batch[id], m.k=m.k[1:k.c], n.k=n.k[1:k.c], PrPTE=PrPTE, rand=rand, q.vec=q.vec, n.a.max=n.a.max[id[-1]], Active = Active[1,id][-1] )
        Accrual[a.i]        = Accrual[a.i] + 1
      }
}

############## follow up: update data after  the active accrual period #####################
    if(n.last< N.A[1]){
      id     = rep(TRUE, N.A[2])

      for( j in (n.last+1):N.A[1]){
        if(any(Active[1,id][-1]==1)){
                                a                    = Potential.Outcome[j,2]
          if(a>1)               Outcome[,a,batch[a]] = Outcome[,a,batch[a]] + c(1L, Potential.Outcome[j,2+a])
          if(a==1){
            for(k in 1:sum(M.k <= i)) if(j>M.k[k])
                                Outcome[,1,k]        = Outcome[,1,k]        + c(1L, Potential.Outcome[j,3])  }

          # apply futility stopping rules
        if(do.futi.stop){
          Outcome.c   = sapply(1:sum(id), function(a) Outcome[,a, batch[a]] )
          PrPTE       = posterior_with_control(Outcome=Outcome.c)
          Active[,id]= apply_stopping_rules(i=i, Active=Active[,id], PrPTE=PrPTE, stop.par=stop.par, Accrual=Accrual[id], n.a.max=n.a.max[id[-1]])
        }
        }}}

    return(list(TrialData=Potential.Outcome, Outcome=Outcome, Active=Active)) } 
