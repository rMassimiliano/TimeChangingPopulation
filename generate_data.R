## This script contains steps to generate clinical trials with a variety of trends 
## data are saved in a folder 'data/'
## We consider
  #- 4 trend scenarios
  #- Multi-arm and platform balanced and Bayesian response-adaptive randomization

## create data folder 
dir.create('data')

## libraries
library(magrittr)
library(TrendUtilities,lib.loc = '.') ## change lib.loc to TrendUtilities location
## we use dqsample to avoid bias in sampling a discrete variable with equal probabilities
library(dqsample)

#+++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++
#++  Multi-arm trials +++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++


#------  NO TREND  ---------------------------#
## BAYESIAN ADAPTIVE RANDOMIZATION NO TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 250
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10

study_coordinates = cbind(c(0,0.25,0.50,0.75,1), qlogis(c(0.1,0.15,0.32,0.35,0.4)))
time_function     = smooth.spline(x = study_coordinates[,1], y = study_coordinates[,2],all.knots =T, df = 5 )

theta_0 = 0

arm_param =  c(theta_0,theta_0 + gamma_a)


for(i in 1:5000)
{
  set.seed(123 + i )
  ## generate data --- BAR 
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))

  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a]))))

  potential =  multi_arm_BAR(potential,rand_par,123 +i)


  resp      =  sapply(1:N, function(i) potential[i,2 + potential[i,2]])
  trialData = data.frame(Y  = resp,
                         A  = potential[,2],
                         T  = potential[,1],
                         T0 = enrollmentTime)
  
  current_name =  sprintf("data/bar_no_trend_datat%i.csv",i)
  write.csv(trialData,file = current_name)

  ###
  if(i%%100 == 0){cat(sprintf("generated datasets # %i \n", i))
}
}


## BALANCED RANDOMIZATION NO TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 250
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10

theta_0 = 0
arm_param =  c(theta_0,theta_0 + gamma_a)


for(i in 1:5000)
{
  set.seed(123 + i )
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  arms           = dqsample::sample(1:A, size = N, replace = TRUE)
  resp           = sapply(1:N,function(i) rbinom(1,1, prob = plogis(arm_param[arms[i]] )))
                          
    trialData = data.frame(Y  = resp,
                           A  = arms,
                           T  = enrollmentTime_01,
                           T0 = enrollmentTime)


     current_name =  sprintf("data/balanced_no_trend_datat%i.csv",i)
     write.csv(trialData,file = current_name)
     if(i%%100==0){cat(sprintf("generated datasets # %i \n", i))}
}


#------  NEGATIVE TREND  ---------------------#
## BAYESIAN ADAPTIVE RANDOMIZATION NEGATIVE TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 250
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10

study_coordinates = cbind(c(0,0.25,0.50,0.75,1), qlogis(c(0.4,0.35,0.32,0.15,0.1)))
time_function     = smooth.spline(x = study_coordinates[,1], y = study_coordinates[,2],all.knots =T, df = 5 )

theta_0 = 1.33913
arm_param =  c(theta_0,theta_0 + gamma_a)


for(i in 1:5000)
{
  set.seed(123 + i )
  ## generate data --- balanced so just sample arms
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t            = predict(time_function, enrollmentTime_01)$y

  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t)))
  )

  potential =  multi_arm_BAR(potential,rand_par,123 +i,omega,N_min)


  resp      =  sapply(1:N, function(i) potential[i,2 + potential[i,2]])
  trialData = data.frame(Y = resp,
                         A = potential[,2],
                         T = potential[,1],
                         T0 = enrollmentTime)
  
   
  current_name =  sprintf("data/bar_neg_trend_datat%i.csv",i)
  write.csv(trialData,file = current_name)
  if(i%%100 ==0)cat(sprintf("generated datasets # %i \n", i))
}




## BALANCED RANDOMIZATION NEGATIVE TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 250
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10



## TREND function --- via points
study_coordinates = cbind(c(0,0.25,0.50,0.75,1), qlogis(c(0.4,0.35,0.32,0.15,0.1)))
time_function     = smooth.spline(x = study_coordinates[,1], y = study_coordinates[,2],all.knots =T, df = 5 )


theta_0 = 1.33913
arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t            = predict(time_function, enrollmentTime_01)$y

  arms           = dqsample::sample(1:A, size = N, replace = TRUE)
  resp           = sapply(1:N,function(i) rbinom(1,1, prob =
                                        plogis(arm_param[arms[i]] + f_t[i])))
                          
     trialData = data.frame(Y = resp,
                        A = arms,
                        T = enrollmentTime_01,
                        T0 = enrollmentTime)

     current_name =  sprintf("data/balanced_neg_trend_datat%i.csv",i)
     write.csv(trialData,file = current_name)
    if(i%%100==0){cat(sprintf("generated datasets # %i \n", i))}
}


#------  POSITIVE TREND  ---------------------#
## BAYESIAN ADAPTIVE RANDOMIZATION POSITIVE TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 250
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10

study_coordinates = cbind(c(0,0.25,0.50,0.75,1), qlogis(c(0.1,0.15,0.32,0.35,0.4)))
time_function     = smooth.spline(x = study_coordinates[,1], y = study_coordinates[,2],all.knots =T, df = 5 )

theta_0 = 0.884

arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t            = predict(time_function, enrollmentTime_01)$y


  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t))))

  potential =  multi_arm_BAR(potential,rand_par,123 +i,omega,N_min)

  resp      =  sapply(1:N, function(i) potential[i,2 + potential[i,2]])
  trialData = data.frame(Y  = resp,
                         A  = potential[,2],
                         T  = potential[,1],
                         T0 = enrollmentTime)


  current_name =  sprintf("data/bar_pos_trend_datat%i.csv",i)
  write.csv(trialData,file = current_name)
  if(i%%100 == 0){cat(sprintf("generated datasets # %i \n", i))}
}



## BALANCED RANDOMIZATION POSITIVE TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 250
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10


study_coordinates = cbind(c(0,0.25,0.50,0.75,1), qlogis(c(0.1,0.15,0.32,0.35,0.4)))
time_function     = smooth.spline(x = study_coordinates[,1], y = study_coordinates[,2],all.knots =T, df = 5 )

theta_0 = 0.884
arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t            = predict(time_function, enrollmentTime_01)$y
                          
  arms           = dqsample::sample(1:A, size = N, replace = TRUE)
  resp           = sapply(1:N,function(i) rbinom(1,1, prob =
                                        plogis(arm_param[arms[i]] + f_t[i])))
                          
     trialData = data.frame(Y = resp,
                        A = arms,
                        T = enrollmentTime_01,
                        T0 = enrollmentTime)
     
   
     current_name =  sprintf("data/balanced_pos_trend_datat%i.csv",i)
     write.csv(trialData,file = current_name)
    if(i%%100==0){cat(sprintf("generated datasets # %i \n", i))}
}


#------  SEASONAL TREND  ---------------------#
## BAYESIAN ADAPTIVE RANDOMIZATION SEASONAL
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


time_function        =  function(t) sin(t*pi/6)
time_function_01     =  function(t,M,m) (M-m)*time_function((M-m)*t + m)

theta_0 = -0.01275168
arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t              = sapply(enrollmentTime, time_function)
  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t)))
  )

  potential =  multi_arm_BAR(potential,rand_par,123 +i,omega,N_min)

  resp      =  sapply(1:N, function(i) potential[i,2 + potential[i,2]])
  trialData = data.frame(Y  = resp,
                         A  = potential[,2],
                         T  = potential[,1],
                         T0 = enrollmentTime)

  current_name =  sprintf("data/bar_ses_trend_datat%i.csv",i)
  write.csv(trialData,file = current_name)
  if(i%%100 == 0){cat(sprintf("generated datasets # %i \n", i))}
}



## BAYESIAN BALANCED RANDOMIZATION SEASONAL
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 250
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10


time_function        =  function(t) sin(t*pi/6)
time_function_01     =  function(t,M,m) (M-m)*time_function((M-m)*t + m)

theta_0 = -0.01275168
arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))

  f_t               = sapply(enrollmentTime, time_function)
  arms           = dqsample::sample(1:A, size = N, replace = TRUE)
  resp           = sapply(1:N,function(i) rbinom(1,1, prob =
                               plogis(arm_param[arms[i]] + f_t[i])))
     
     trialData = data.frame(Y = resp,
                            A = arms,
                            T = enrollmentTime_01,
                            T0 = enrollmentTime)
   
     current_name =  sprintf("data/balanced_ses_trend_datat%i.csv",i)
     write.csv(trialData,file = current_name)
    if(i%%100 == 0) cat(sprintf("generated datasets # %i \n", i))
}

#+++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++
#++  Platform trials +++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++

#------  NO TREND  ---------------------#

## BAYESIAN PLATFORM NO TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 300
rates         = rep(2:6, each =50)
rates = c(rates, rep(6, N - length(rates)))

## design parameters
rand_par      = c(5,1.5,1)
arm_added_at  = c(0,0,0, 50,100)


theta_0 = 0

arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))

  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob =plogis(arm_param[a])))
  )

## Platform BAR 
 Trial = simulate_BAR_add_trial( 
     Potential.Outcome = potential,
     added.at          = arm_added_at,
     Outcome.delay     = 0,
     n.a.max           = rep(65, 4),
     rand              = rand_par,
     q.vec             = c(0.5,7, 15),
     m.k               = c(20, 20,20),
     n.k               = c(3*50, 50,50),
     stop.par          = c(f=.3, gf=1),
     do.futi.stop      = FALSE)

  trialData = data.frame(Y = sapply(1:N, function(i) Trial$TrialData[i,-c(1:2)][Trial$TrialData[i,2] ]),
               A = Trial$TrialData[,2],
               T = Trial$TrialData[,1],
              T0 = enrollmentTime_01)

              
  ## save data
  write.csv(trialData,sprintf("data/plat_bar_no_trend_dataset%i.csv",i), row.names = FALSE) 
  ## progress
  if(i%%100 == 0)cat(sprintf("Done iteration %i \n",i))
} 



## BAYESIAN PLATFORM BALANCED NO TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 300
rates         = rep(2:6, each =50)
rates = c(rates, rep(6, N - length(rates)))
## design parameters
rand_par      = c(0,1.5,1)
arm_added_at  = c(0,0,0, 50,100)
theta_0 = 0
arm_param =  c(theta_0,theta_0 + gamma_a)


for(i in 1:5000)
{
  set.seed(123 + i )
  ## generate data 
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))

  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a])))
  )

## Platform BALANCED
Trial = simulate_BAR_add_trial( 
    Potential.Outcome = potential,
    added.at          = arm_added_at,
    Outcome.delay     = 0,
    n.a.max           = rep(65, 4),
    rand              = rand_par,
    q.vec             = c(0.5,7, 15),
    m.k               = c(20, 20,20),
    n.k               = c(3*50, 50,50),
    stop.par          = c(f=.3, gf=1),
    do.futi.stop      = FALSE)

  trialData = data.frame(Y = sapply(1:N, function(i) Trial$TrialData[i,-c(1:2)][Trial$TrialData[i,2] ]),
               A = Trial$TrialData[,2],
               T = Trial$TrialData[,1],
               T0 = enrollmentTime)

              
  ## save data
  write.csv(trialData,sprintf("data/plat_balanced_no_trend_dataset%i.csv",i), row.names = FALSE) 
  ## progress
  if(i%%100 == 0)cat(sprintf("Done iteration %i \n",i))
} 


#------  NEGATIVE TREND  ---------------------#
## PLATFORM BAR NEGATIVE TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 300
rates         = rep(2:6, each =50)
rand_par      = c(5,1.5,1)
omega         = 7
N_min         = 10

## design parameters
rand_par      = c(5,1.5,1)
arm_added_at  = c(0,0,0, 50,100)

## time trend
study_coordinates = cbind(c(0,0.25,0.50,0.75,1), qlogis(c(0.4,0.35,0.32,0.15,0.1)))
time_function     = smooth.spline(x = study_coordinates[,1], y = study_coordinates[,2],all.knots =T, df = 5 )

theta_0 = 1.33913
arm_param =  c(theta_0,theta_0 + gamma_a)


for(i in 1:5000)
{
  set.seed(123 + i )
  ## enrollment times
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t            = predict(time_function, enrollmentTime_01)$y


  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t)))
  )

## Platform BAR  
Trial = simulate_BAR_add_trial( 
    Potential.Outcome = potential,
    added.at          = arm_added_at,
    Outcome.delay     = 0,
    n.a.max           = rep(65, 4), 
    rand              = rand_par,
    q.vec             = c(0.5,7, 15),
    m.k               = c(20, 20,20),
    n.k               = c(3*50, 50,50),
    stop.par          = c(f=.3, gf=1),
    do.futi.stop      = FALSE)

  trialData = data.frame(Y = sapply(1:N, function(i) Trial$TrialData[i,-c(1:2)][Trial$TrialData[i,2] ]),
               A = Trial$TrialData[,2],
               T = Trial$TrialData[,1],
               T0 = enrollmentTime) 
              
  ## save data
  write.csv(trialData,sprintf("data/bar_neg_trend_dataset%i.csv",i), row.names = FALSE) 
  ## progress
  if(i%%100 == 0)cat(sprintf("Done iteration %i \n",i))
} 



## PLATFORM BALANCED NEGATIVE TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 300
rates         = rep(2:6, each =50)
rates = c(rates, rep(6, N - length(rates)))

## design parameters
rand_par      = c(0,1.5,1)
arm_added_at  = c(0,0,0, 50,100)

## time trend
study_coordinates = cbind(c(0,0.25,0.50,0.75,1), qlogis(c(0.4,0.35,0.32,0.15,0.1)))
time_function     = smooth.spline(x = study_coordinates[,1], y = study_coordinates[,2],all.knots =T, df = 5 )

theta_0 = 1.33913
arm_param =  c(theta_0,theta_0 + gamma_a)


for(i in 1:5000)
{
  set.seed(123 + i )
  ## enrollment times

  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t            = predict(time_function, enrollmentTime_01)$y

  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t)))
  )

## Platform BAR  
Trial = simulate_BAR_add_trial( 
    Potential.Outcome = potential,
    added.at          = arm_added_at,
    Outcome.delay     = 0,
    n.a.max           = rep(65, 4), 
    rand              = rand_par,
    q.vec             = c(0.5,7, 15),
    m.k               = c(20, 20,20),
    n.k               = c(3*50, 50,50),
    stop.par          = c(f=.3, gf=1),
    do.futi.stop      = FALSE)


  trialData = data.frame(Y = sapply(1:N, function(i) Trial$TrialData[i,-c(1:2)][Trial$TrialData[i,2] ]),
               A = Trial$TrialData[,2],
               T = Trial$TrialData[,1],
               T0 =   enrollmentTime) 
              
  ## save data
  write.csv(trialData,sprintf("data/balanced_neg_trend_dataset%i.csv",i), row.names = FALSE) 
  ## progress
  if(i%%100 == 0)cat(sprintf("Done iteration %i \n",i))
} 



#------ POSITIVE TREND  ---------------------#
## PLATFORM BAR POSITIVE TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 300
rates         = rep(2:6, each =50)
rates = c(rates, rep(6, N - length(rates)))

## design parameters
rand_par      = c(5,1.5,1)
arm_added_at  = c(0,0,0, 50,100)

## time trend
study_coordinates = cbind(c(0,0.25,0.50,0.75,1), qlogis(c(0.1,0.15,0.32,0.35,0.4)))
time_function     = smooth.spline(x = study_coordinates[,1], y = study_coordinates[,2],all.knots =T, df = 5 )

theta_0 = 0.884
arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  ## enrollment times
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t            = predict(time_function, enrollmentTime_01)$y

  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t)))
  )

## Platform BAR  
Trial = simulate_BAR_add_trial( 
    Potential.Outcome = potential,
    added.at          = arm_added_at,
    Outcome.delay     = 0,
    n.a.max           = rep(65, 4), 
    rand              = rand_par,
    q.vec             = c(0.5,7, 15),
    m.k               = c(20, 20,20),
    n.k               = c(3*50, 50,50),
    stop.par          = c(f=.3, gf=1),
    do.futi.stop      = FALSE)


  trialData = data.frame(Y = sapply(1:N, function(i) Trial$TrialData[i,-c(1:2)][Trial$TrialData[i,2] ]),
               A = Trial$TrialData[,2],
               T = Trial$TrialData[,1],
               T0 =  enrollmentTime)

              
  ## save data
  write.csv(trialData,sprintf("data/bar_pos_trend_dataset%i.csv",i), row.names = FALSE) 
  ## progress
  if(i%%100 == 0)cat(sprintf("Done iteration %i \n",i))
} 


## PLATFORM BALANCED POSITIVE TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 300
rates         = rep(2:6, each =50)
rates = c(rates, rep(6, N - length(rates)))

## design parameters
rand_par      = c(0,1.5,1)
arm_added_at  = c(0,0,0, 50,100)

## time trend
study_coordinates = cbind(c(0,0.25,0.50,0.75,1), qlogis(c(0.1,0.15,0.32,0.35,0.4)))
time_function     = smooth.spline(x = study_coordinates[,1], y = study_coordinates[,2],all.knots =T, df = 5 )

theta_0 = 0.884
arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  ## enrollment times
  enrollmentTime = cumsum(rexp(N, rate = rates))
  f_t            = predict(time_function, enrollmentTime)$y

  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t            = predict(time_function, enrollmentTime_01)$y

  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t)))
  )

## Platform BAR  
Trial = simulate_BAR_add_trial( 
    Potential.Outcome = potential,
    added.at          = arm_added_at,
    Outcome.delay     = 0,
    n.a.max           = rep(65, 4),
    rand              = rand_par,
    q.vec             = c(0.5,7, 15),
    m.k               = c(20, 20,20),
    n.k               = c(3*50, 50,50),
    stop.par          = c(f=.3, gf=1),
    do.futi.stop      = FALSE)


  trialData = data.frame(Y = sapply(1:N, function(i) Trial$TrialData[i,-c(1:2)][Trial$TrialData[i,2] ]),
               A = Trial$TrialData[,2],
               T = Trial$TrialData[,1],
               T0 = enrollmentTime)

              
  ## save data
  write.csv(trialData,sprintf("data/balanced_pos_trend_dataset%i.csv",i), row.names = FALSE) 
  ## progress
  if(i%%100 == 0)cat(sprintf("Done iteration %i \n",i))
} 



#------ SEASONAL TREND  ---------------------#
## PLATFORM BAR SEASONAL TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 300
rates         = rep(2:6, each =50)
rates = c(rates, rep(6, N - length(rates)))

## design parameters
rand_par      = c(5,1.5,1)
arm_added_at  = c(0,0,0, 50,100)

## time trend
time_function        =  function(t) sin(t*pi/6)
time_function_01     =  function(t,M,m) (M-m)*time_function((M-m)*t + m)

theta_0 = -0.01275168
arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  ## enrollment times
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))
  f_t              = sapply(enrollmentTime, time_function)

  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t)))
  )

## Platform BAR  
Trial = simulate_BAR_add_trial( 
    Potential.Outcome = potential,
    added.at          = arm_added_at,
    Outcome.delay     = 0,
    n.a.max           = rep(65, 4), 
    rand              = rand_par,
    q.vec             = c(0.5,7, 15),
    m.k               = c(20, 20,20),
    n.k               = c(3*50, 50,50),
    stop.par          = c(f=.3, gf=1),
    do.futi.stop      = FALSE)


  trialData = data.frame(Y = sapply(1:N, function(i) Trial$TrialData[i,-c(1:2)][Trial$TrialData[i,2] ]),
               A = Trial$TrialData[,2],
               T = Trial$TrialData[,1],
               T0 = enrollmentTime)
  

              
  ## save data
  write.csv(trialData,sprintf("data/bar_ses_trend_dataset%i.csv",i), row.names = FALSE) 
  ## progress
  if(i%%100 == 0)cat(sprintf("Done iteration %i \n",i))
} 


## PLATFORM BALANCED SEASONAL TREND
gamma_a        = c(0.8472979,0,1.098612,0)
A             = length(gamma_a)+1
N             = 300
rates         = rep(2:6, each =50)
rates = c(rates, rep(6, N - length(rates)))

## design parameters
rand_par      = c(0,1.5,1)
arm_added_at  = c(0,0,0, 50,100)

## time trend
time_function        =  function(t) sin(t*pi/6)
time_function_01     =  function(t,M,m) (M-m)*time_function((M-m)*t + m)

theta_0 = -0.01275168
arm_param =  c(theta_0,theta_0 + gamma_a)

for(i in 1:5000)
{
  set.seed(123 + i )
  ## enrollment times
  enrollmentTime = cumsum(rexp(N, rate = rates))
  enrollmentTime_01 = (enrollmentTime -min( enrollmentTime))/(max(enrollmentTime) -min( enrollmentTime))

  f_t               = sapply(enrollmentTime, time_function)

  potential = cbind(T = enrollmentTime_01,A =  0, Y = 
                    sapply(1:A, function(a) rbinom(N,1, prob = plogis(arm_param[a] + f_t)))
  )

## Platform BAR  
Trial = simulate_BAR_add_trial( 
    Potential.Outcome = potential,
    added.at          = arm_added_at,
    Outcome.delay     = 0,
    n.a.max           = rep(65, 4), 
    rand              = rand_par,
    q.vec             = c(0.5,7, 15),
    m.k               = c(20, 20,20),
    n.k               = c(3*50, 50,50),
    stop.par          = c(f=.3, gf=1),
    do.futi.stop      = FALSE)


  trialData = data.frame(Y = sapply(1:N, function(i) Trial$TrialData[i,-c(1:2)][Trial$TrialData[i,2] ]),
               A = Trial$TrialData[,2],
               T = Trial$TrialData[,1],
               T0 = enrollmentTime)
  

              
  ## save data
  write.csv(trialData,sprintf("data/balanced_ses_trend_dataset%i.csv",i), row.names = FALSE) 
  ## progress
  if(i%%100 == 0)cat(sprintf("Done iteration %i \n",i))
} 
############################################################################################
