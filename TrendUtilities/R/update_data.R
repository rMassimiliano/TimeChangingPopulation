#' Update outcome data during the active accrual period and during the follow up time.
#'
#' @param i                  current patients
#' @param t.i                current trial time
#' @param Potential.Outcome  matrix of potential outcome data
#' @param Outcome.delay      Waiting time until each outcome is available
#' @param Outcome            Array of outcome data
#' @param n.last             First n.last patient data where already available
#' @param batch              Vector indicating for each arm the adding-group A.k
#' @param N.target.0         Initial target sample size for the control arm.
#' @param M.k                Vector of unique adding times.
#' @param during.accrual     TRUE is update during active accrual period, FALSE if update during follow up perior.
#'
#' @keywords internal
#' @export
update_data = function(i, t.i, Potential.Outcome, Outcome.delay, Outcome, n.last, batch, N.target.0, M.k, during.accrual=TRUE)
{
  j.before               = Potential.Outcome[(n.last+1):i,1] + Outcome.delay <= t.i

  if(any(j.before))
  {
    j.before = which(j.before) + max(n.last-1,0) 
    n.last                = max(tail(j.before, n=1), i)

    for(j in j.before)
    {    
      a                    = Potential.Outcome[j,2] 
      if(a>1)
      {
          Outcome[,a,batch[a]] = Outcome[,a,batch[a]] + c(1L, Potential.Outcome[j,2+a])
      }
      if(a==1)
      { 
       for(k in 1:sum(M.k<=i))
       {
        if(j>M.k[k] & Outcome[1,1,k]<N.target.0)
        {
         Outcome[,1,k]        = Outcome[,1,k] + c(1L, Potential.Outcome[j,3])            
        } 
       } 
      }
    }
  }

return(list(Outcome=Outcome, n.last=n.last)) }


 
