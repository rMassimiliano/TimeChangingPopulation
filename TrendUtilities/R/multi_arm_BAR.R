#' Simulate a Trial using Bayesian response-adaptive (BAR) randomization
#'
#' Generate a multi-arm trial under BAR 
#' @param potential  potential outcome matrix each row as
#'                  (enrollment_time,Assigned_place_holder, Y_0, ...,Y_A  )
#' @param rand       A vector \code{rand =c(H, gamma, b)} of randomization three tuning parameter.
#'                   The first and second element correspond to the maximum
#'                   and the shap parameter of the power function \eqn{h(x) = H * (x/n)^gamma}
#'                   The last parameter controls the randomization to the control.
#'@param seed  include a seed for the simulation.
#'@param omega tuning parameter for the exponential function to guarantee approximately N_min patients in each experimental arm.
#'@param N_min Approximate minimum for each experimental arm.
#'@param outcome_delay After how much time the outcome of patients \code{i} become available for the randomization procedure. 
#' @return
#' The function return a matrix with the same entrance as \code{potential}. The second column includs allocation of patients to arms according to the BAR procedure.
multi_arm_BAR = function(potential, rand,seed = NULL, omega =0,N_min = 0,outcome_delay =0)
{
  N = NROW(potential)
  A = NCOL(potential) -2
  if(!is.null(seed)){set.seed(seed)}
  potential[1,2] = sample(1:A, size =1)
  for(i in 2:N)
  {
    ## indices of the obseved response up to i-1
    ind = which( (potential[i,1] - potential[1:(i-1),1]) >= outcome_delay)

    if(length(ind))
    {
    suffStat = rbind(N = tabulate(potential[ind,2],A),
                    Y = table( factor(sapply(ind, function(x) potential[x,potential[x,2] + 2]),0:1), factor(potential[ind,2],1:A))[2,])
    

    PrPTE = drop(posterior_with_control(suffStat))
    ## inflation for minimum number of patients
    min_cor = exp(pmax(omega*(N_min - suffStat[1,-1]),0))
    h_i   = rand[1]*((i/N)^rand[2])
    Pa    = PrPTE^h_i
    P0    = exp(rand[3]*(max(suffStat[1,-1])) - suffStat[1,1])
    Pr    = c(1,min_cor) *  c(P0, Pa/sum(Pa) )
    Pr    = Pr/sum(Pr)
    potential[i,2] =  sample(1:A, size = 1, prob =Pr)
   }
   else
   { ## no accrued info, generate from the prior, which is approx balanced 

     suffStat = rbind(N =rep(0,A),
                      Y =  rep(0,A))
     PrPTE = drop(posterior_with_control(suffStat))
     ## inflation for minimum number of patients
     min_cor = exp(pmax(omega*(N_min - suffStat[1,-1]),0))
     h_i   = rand[1]*((i/N)^rand[2])
     Pa    = PrPTE^h_i
     P0    = exp(rand[3]*(max(suffStat[1,-1])) - suffStat[1,1])
     Pr    = c(1,min_cor) *  c(P0, Pa/sum(Pa) )
     Pr    = Pr/sum(Pr)
     potential[i,2] =  sample(1:A, size = 1, prob =Pr)
   }
  }
  potential
}
