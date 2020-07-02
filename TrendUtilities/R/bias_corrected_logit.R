bias_corrected_logit =
function(
Y_0,   ## response in controll group
Y_1,   ## response in treatment group
T_0,   ## time in controll group  dim  -- matrix
T_1   ## time in treatment group dim  -- matrix
)
{
## sample sizes
n_0       = NROW(Y_0)
n_1       = NROW(Y_1)
## estimated quanities
Sigma_0   = var(T_0)                   #matrix
Sigma_1   = var(T_1)                   #matrix
Cov_YT_0  = cov(Y_0,T_0)               #vector
Cov_YT_1  = cov(Y_1,T_1)               #vector

hat_mu_0  = mean(Y_0)
hat_mu_1  = mean(Y_1)

hat_theta = log(hat_mu_1) - log(hat_mu_0) + log(1 - hat_mu_0) - log(1 - hat_mu_1)

## corrections
hat_delta = colMeans(T_1) - colMeans(T_0)
sigma_11  = 1/(n_0 * hat_mu_0 *(1 - hat_mu_0)) + 1/(n_1 * hat_mu_1 *(1 - hat_mu_1))
sigma_12  = Cov_YT_1/ (n_1 * hat_mu_1 *(1 - hat_mu_1)) + Cov_YT_0/(n_0 * hat_mu_0 *(1 - hat_mu_0))
sigma_22  = Sigma_1/n_1 +  Sigma_0/n_0


theta_adj = drop(hat_theta - sigma_12 %*% (solve(sigma_22)) %*% hat_delta)
est_var   = drop(sigma_11 - sigma_12 %*%  solve(sigma_22) %*% t(sigma_12))

 return(c("theta" = theta_adj, "var" = est_var))
}
