#' Compute a standardized difference in proportion to perform Z-test
#' 
#'
#'@param Y_0 response vector for the constrol group
#'@param Y_1 response vector for the treatment group
#'@return The standardized difference in proportion
#'@export
Z_prop = function(Y_0,Y_1)
{
  n_0 = length(Y_0)
  n_1 = length(Y_1)

  p_0 = mean(Y_0)
  p_1 = mean(Y_1)

  p   = mean(c(Y_0,Y_1))

  Z   = (p_1 - p_0)/sqrt(p*(1-p) *(1/n_0 + 1/n_1))
  return(Z)
}
