#' Compute an Adjusted Z statistic conditioning on the discrepancy measure T_1 -T_0 
#'
#' Compute the adjusted statistic \eqn{\widehat \gamma_a^{\mbox{\tiny adj}}  = \widehat \gamma_a - \widehat{\mathbb E} [\widehat \gamma_a  - \gamma_a|   \Delta_a]}, where \eqn{\Delta_a = \bar{T}_a - \bar{T_0}}, is the mean difference of T between the treatment and the control group.
#' 
#' 
#'@param Y_0 response in the control group
#'@param Y_1 response in the treatment group
#'@param T_0 enrollment time (or a patient level covariate) for the control group
#'@param T_1 enrollment time (or a patient level covariate) for the treatment group
#' @return 
#' Returns a vector including the treatment effect estimate and its variance
#' @export
Adj_Z =function(
                Y_0,   
                Y_1,   
                T_0,   ## time(or covariates) in controll group  dim  -- matrix
                T_1    ## time(covariates) in treatment group dim  -- matrix
)
{
 if(is.null(ncol(T_0)) & is.null(ncol(T_1)))
 {
   T_0 = as.matrix(T_0)
   T_1 = as.matrix(T_1)
 }
  
## sample sizes
n_0       = NROW(Y_0)
n_1       = NROW(Y_1)
## estimated quantitises
Sigma_0   = var(T_0)                   
Sigma_1   = var(T_1)                   
Cov_YT_0  = cov(Y_0,T_0)               
Cov_YT_1  = cov(Y_1,T_1)               

## means
hat_mu_0  = mean(Y_0)
hat_mu_1  = mean(Y_1)
hat_T_0   = colMeans(T_0)
hat_T_1   = colMeans(T_1)

## empirical  estimate
hat_theta = hat_mu_1 - hat_mu_0 

## discrepancy measure
hat_delta = hat_T_1 - hat_T_0

## corrections
sigma_11  = 1/(n_0 * hat_mu_0 *(1 - hat_mu_0)) + 1/(n_1 * hat_mu_1 *(1 - hat_mu_1))
sigma_12  = Cov_YT_1/n_1 + Cov_YT_0/n_0
sigma_22  = Sigma_1/n_1 +  Sigma_0/n_0

## adjusted estimate
theta_adj = drop(hat_theta - sigma_12 %*% (solve(sigma_22)) %*% hat_delta)
## estimated variance

## arm 1
A1 = hat_T_1
A2 = hat_mu_1 
A3 = mean((T_1 - hat_T_1) * (Y_1 - hat_mu_1))  
A4 = mean((T_1 - hat_T_1)**2)

B1 = hat_T_0
B2 = hat_mu_0 
B3 = mean((T_0 - hat_T_0) * (Y_0 - hat_mu_0))  
B4 = mean((T_0 - hat_T_0)**2)
   


#variables  definition 
AA1 = T_1
AA2 = Y_1
AA3 = (T_1 - hat_T_1)*(Y_1 - hat_mu_1)
AA4 = (T_1 - hat_T_1)^2
AAA = cbind(AA1, AA2, AA3, AA4)


BB1 = T_0
BB2 = Y_0
BB3 = (T_0 - hat_T_0)*(Y_0 - hat_mu_0)
BB4 = (T_0 - hat_T_0)^2

BBB = cbind(BB1, BB2, BB3, BB4)

#estimated cov
SAAA=cov(AAA)/n_1
SBBB=cov(BBB)/n_0
S=array(0,c(8,8)); S[1:4,1:4]=SAAA;S[5:8,5:8]=SBBB;

#derivatives for delta method
AD=c(-(((n_1/n_0)*B3+A3 )  /((n_1/n_0)*B4+A4)),
1,
-(A1-B1)*((1 )  /((n_1/n_0)*B4+A4)),
(A1-B1)*(((n_1/n_0)*B3+A3 )  /((n_1/n_0)*B4+A4)^2))
BD=c((((n_1/n_0)*B3+A3 )  /((n_1/n_0)*B4+A4)),-1,-(A1-B1)*(((n_1/n_0) )  /((n_1/n_0)*B4+A4)),
(n_1/n_0)*(A1-B1)*(((n_1/n_0)*B3+A3 )  /((n_1/n_0)*B4+A4)^2))

#delta method application
est_var = drop(c(AD,BD)%*%S%*%c(AD,BD))

return(c("theta" = theta_adj, "var" = est_var))
}



