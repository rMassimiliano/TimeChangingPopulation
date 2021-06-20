
# Inference in response-adaptive clinical trials when the enrolled population varies over time

This repository contains an R package including code used to generate
and analyze clinical trials in the paper *Inference in response-adaptive
clinical trials when the enrolled population varies over time*. The
manuscript is available upon request.

The package <tt>TrendUtilities</tt> can be installed with the following
code

``` r
R CMD build TrendUtilities
R CMD INSTALL TrendUtilities -l your_path
```

or alternatively using <tt>devtools</tt>.

The package contains utilities to replicate the simulations studies
presented in Section 4 of the paper.

## Two-arm trials

We consider as an example the linear trend (Trend 1) used in Section 4.1
of the paper without treatment effect (i.e, \(gamma_1 =0\)). We allocate
patient to treatment using the design BAR-1 (function <tt>
bar2arms</tt>).

As a first step we load the package <tt>TrendUtilities</tt> assuming
that it has been installed in the folder <tt>local\_libraries</tt> in
the current path. We also load the R package <tt>mgcv</tt> used to
estimate the GAM model described in Section 3.1 of the paper.

``` r
library(TrendUtilities, lib.loc = 'local_libraries/')
library(mgcv)
```

We now create a matrix of potential outcomes with the following code

``` r
N=100 ## trial size
M =20 ## allocation probabilities are updated any 20 enrollments
T_i = seq(0,1,l=N)
potential = matrix(0, N, 4)
potential[,1] = T_i

## potential outcomes for treatment and control group
Y  = rbinom(N,1, 0.5 + 0.3*T_i) 
Y1 = rbinom(N,1, 0.5 + 0.3*T_i)

potential[,3]  = Y
potential[,4]  = Y1

## allocate patients to treatment via BAR-1 algorithm
potential = bar2arms(potential, c = 1/3, M = M, ret_prob =FALSE)
A = potential[,2]

## observed outcomes
Y_obs       = sapply(1:N, function(i) potential[i, potential[i,2] + 2 ] )

##trial data
simData = data.frame(Y = Y_obs, A = A, T_i = T_i) 


## Z-test
Z_test   = Z_prop(Y_obs[A==1],Y_obs[A==2])

## A-GAM test
mod  = gam(Y~factor(A) + s(T_i),data = simData, family = binomial)
AGAM = coef(mod)[2]/sqrt(vcov(mod, freq  =TRUE)[2,2])


## Adjusted Z-test
adjZ_stat  =  Adj_Z(Y_obs[A==1],Y_obs[A==2], T_i[A==1],T_i[A==2])
adjZ_test  =  adjZ_stat['theta']/sqrt(adjZ_stat['var'])

# The p-values of the three tests are
pvalues = c(1-pnorm(Z_test),1-pnorm(AGAM),1-pnorm(adjZ_test))
names(pvalues) = c('Z test', 'A-GAM', 'Adj-Z')
round(pvalues,4)
```

    ## Z test  A-GAM  Adj-Z 
    ## 0.8133 0.6582 0.6670

## multi-arm and platform trials

The file <tt>generata\_data.R</tt> include steps to generate all the
clinical trials used for the simulations in Section 4.1 and 4.2. The
trend function used in this trails are described in the Supplementary
Material.

We consider here and example with seasonal trend and multi-arm BAR

``` r
## simulations parameters
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 250
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10


arm_probs     = c(0.378, 0.500, 0.378,0.769,0.378)
A             = length(arm_probs)
N             = 250
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10


## time trend
time_function        =  function(t) sin(t*pi/6)


theta_0 = -0.01275168
arm_param =  c(theta_0,theta_0 + gamma_a)

## generate enrollment times
enrollmentTime = cumsum(rexp(N, rate = rates))
f_t              = sapply(enrollmentTime, time_function)

#potential outcomes
potential = cbind(T = enrollmentTime,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t)))
  )

## randomize patients to treatments
potential =  multi_arm_BAR(potential,rand_par,123,omega,N_min)

resp      =  sapply(1:N, function(i) potential[i,2 + potential[i,2]])
trialData = data.frame(Y  = resp,
                       A  = potential[,2],
                       T_i  = potential[,1])
                      
## Z test
Z_test   = sapply(2:A, function(a)  with(trialData, Z_prop(Y[A==1],Y[A==a])))

## A-GAM test
mod  = gam(Y~factor(A) + s(T_i),data = trialData, family = binomial)
AGAM = coef(mod)[2:A]/ sqrt(diag(vcov(mod, freq  =TRUE)[2:A,2:A])) 


## Adjusted Z-test
adjZ_test = numeric(A-1)
for(a in 2:A)
{
 adjZ_stat  =  with(trialData, Adj_Z(Y[A==1],Y[A==a], T_i[A==1],T_i[A==a]))
 adjZ_test[a-1] =  adjZ_stat['theta']/sqrt(adjZ_stat['var'])
}
pvalues = rbind(1-pnorm(Z_test), 1-pnorm(AGAM), 1-pnorm(adjZ_test))

colnames(pvalues) = paste('Arm ', 1:4)
rownames(pvalues) = c('Z test', 'A-GAM', 'Adj-Z')

# pvalues
round(pvalues,4)
```

    ##        Arm  1 Arm  2 Arm  3 Arm  4
    ## Z test 0.0083 0.1151 0.0011 0.3993
    ## A-GAM  0.0091 0.1387 0.0017 0.4802
    ## Adj-Z  0.0073 0.1204 0.0014 0.5553
