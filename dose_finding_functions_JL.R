library(survey)
library(tmvtnorm)
library(plyr)
library(tidyverse)
library(rlist)







### JL model posterior distributions post logit transformation ###

posterior_function_E = function(theta_E, ye, level) {
  dose_post = c()
  one_minus_form = c( (1 - exp(theta_E[1])/(1+exp(theta_E[1]))), (1 - exp(theta_E[2])/(1+exp(theta_E[2]))), (1 - exp(theta_E[3])/(1+exp(theta_E[3]))), (1 - exp(theta_E[4])/(1+exp(theta_E[4]))), (1 - exp(theta_E[5])/(1+exp(theta_E[5])))   )
  
  for (j in level:5){
    if (length(ye[[j]]) < 1){
      next
    }
    else{
      for (i in 1:length(ye[[j]])){
        post_value = (1 - prod(one_minus_form[1:j]))^ye[[j]][i] * prod(one_minus_form[1:j])^(1 - ye[[j]][i])
        dose_post = list.append(dose_post, post_value)
      }
      
    }
  }  
  prod(dose_post) * exp(theta_E[level])/(1+exp(theta_E[level]))^2 *(exp(theta_E[level])/(1+exp(theta_E[level])))^(ae[level]-1) *(1 - exp(theta_E[level])/(1+exp(theta_E[level])))^(be[level]-1)   
}


posterior_function_T = function(theta_T, yt, level) {
  dose_post = c()
  one_minus_form = c( (1 - exp(theta_T[1])/(1+exp(theta_T[1]))), (1 - exp(theta_T[2])/(1+exp(theta_T[2]))), (1 - exp(theta_T[3])/(1+exp(theta_T[3]))), (1 - exp(theta_T[4])/(1+exp(theta_T[4]))), (1 - exp(theta_T[5])/(1+exp(theta_T[5])))   )
  
  for (j in level:5){
    if (length(yt[[j]]) < 1){
      next
    }
    else{
      for (i in 1:length(yt[[j]])){
        post_value = (1 - prod(one_minus_form[1:j]))^yt[[j]][i] * prod(one_minus_form[1:j])^(1-yt[[j]][i])
        dose_post = list.append(dose_post, post_value)
      }
    }
  } 
  prod(dose_post) * exp(theta_T[level])/(1 + exp(theta_T[level]))^2 * (exp(theta_T[level])/(1+exp(theta_T[level])))^(at[level]-1) *(1 - exp(theta_T[level])/(1+exp(theta_T[level])))^(bt[level]-1)
}



### Adaptive Metropolis Hastings ###
randomvalue_generator <- function(x, sigma) rnorm(1, x, sigma)

#SINGLE BETA UPDATES
step1 <- function(theta_E, f_E, randomvalue, ye, j, sigma) {
  ## Pick new point
  xp <- theta_E
  xp[j] = randomvalue(theta_E[j], sigma[j])
  old = f_E(theta_E, ye, level = j)
  new = f_E(xp, ye, level = j)
  ## Acceptance probability:
  que <- min(1, new/ old)
  if (is.na(que)){
    que = 0
  }
  ## Accept new point with probability alpha:
  if (runif(1) < que) {
    theta_E <- xp
  }
  ## Returning the point:
  theta_E
}

#DO BETAS ONE AT A TIME
#ADD ADAPTIVE PERHAPS OR OTHER EXPERMENTAION FOR PROPOSAL VARIANCE

step2 <- function(theta_T, f_T, randomvalue, yt, j, sigma) {
  ## Pick new point
  xp <- theta_T
  xp[j] = randomvalue(theta_T[j], sigma[j])
  old = f_T(theta_T, yt, level = j)
  new = f_T(xp, yt, level = j)
  
  #print(xp)
  ## Acceptance probability:
  que <- min(1, new/ old)
  if (is.na(que)){
    que = 0
  }
  #print(que)
  ## Accept new point with probability alpha:
  if (runif(1) < que) {
    theta_T <- xp
  }
  ## Returning the point:
  theta_T
}



run_MCMC <- function(theta_E, theta_T, f_E, f_T, randomvalue, nsteps, ye, yt) {
  rese <- matrix(NA, nrow = nsteps, ncol = 5)
  rest <- matrix(NA, nrow = nsteps, ncol = 5)
  count_change = 0
  countE_change = c(0,0,0,0,0)
  countT_change = c(0,0,0,0,0)
  total_change= 0
  countE = c(0,0,0,0,0)
  countT = c(0,0,0,0,0)
  total = 0
  BE_OLD = c()
  BE_NEW = c()
  BT_OLD = c()
  BT_NEW = c()
  
  sigmaE = c(3,3,3,3,3)
  sigmaT = c(3,3,3,3,3)
  for (i in 1:nsteps){
    total = total + 1
    total_change = total_change+1
    count_change = count_change +1
    BE_OLD[1] = theta_E[1]
    BE_OLD[2] = theta_E[2]
    BE_OLD[3] = theta_E[3]
    BE_OLD[4] = theta_E[4]
    BE_OLD[5] = theta_E[5]
    BT_OLD[1] = theta_T[1]
    BT_OLD[2] = theta_T[2]
    BT_OLD[3] = theta_T[3]
    BT_OLD[4] = theta_T[4]
    BT_OLD[5] = theta_T[5]
    #adjust sigma^2 for proposal distribution
    for(j in  1:5){
      if (count_change %% 50 == 0){
        changer = min(0.1, 1/total^(1/3))
        statusE = ((countE_change[j])/total_change)/0.5
        if (statusE >= 1){
          sigmaE[j] = sigmaE[j] + changer
        }else{
          sigmaE[j] = sigmaE[j] - changer
        }
        if (sigmaE[j] >8){
          sigmaE[j] = 8
        }
        if(sigmaE[j] < 1){
          sigmaE[j] = 1
        }
        if (countE_change[j]/total_change < 0.2){
          sigmaE[j] = 1
        }
        
        statusT = ((countT_change[j])/total_change)/0.5
        if (statusT >= 1){
          sigmaT[j] = sigmaT[j] + changer
        }else{
          sigmaT[j] = sigmaT[j] - changer
        }
        if (sigmaT[j] >8){
          sigmaT[j] = 8
        }
        if(sigmaT[j] < 1){
          sigmaT[j] = 1
        }
        if (countT_change[j]/total_change < 0.2){
          sigmaT[j] = 1
        }
        
        
      }
      
      
      theta_E <- step1(theta_E, f_E, randomvalue, ye, j, sigma = sigmaE)
      rese[i,j] <- theta_E[j]
      theta_T <- step2(theta_T, f_T, randomvalue, yt,j, sigma= sigmaT)
      rest[i,j] <- theta_T[j]
      
    }
    if (count_change %% 50 == 0){
      countE_change = c(0,0,0,0,0)
      countT_change = c(0,0,0,0,0)
      total_change= 0
    }
    
    BE_NEW[1] <- theta_E[1]
    BE_NEW[2] <- theta_E[2]
    BE_NEW[3] <- theta_E[3]
    BE_NEW[4] <- theta_E[4]
    BE_NEW[5] <- theta_E[5]
    BT_NEW[1] <- theta_T[1]
    BT_NEW[2] <- theta_T[2]
    BT_NEW[3] <- theta_T[3]
    BT_NEW[4] <- theta_T[4]
    BT_NEW[5] <- theta_T[5]
    #Check for acceptance rate
    for (t in 1:5){
      if (BE_OLD[t] != BE_NEW[t]){countE[t] = countE[t]+1}
      if (BT_OLD[t] != BT_NEW[t]){countT[t] = countT[t]+1}
      if (BE_OLD[t] != BE_NEW[t]){countE_change[t] = countE_change[t]+1}
      if (BT_OLD[t] != BT_NEW[t]){countT_change[t] = countT_change[t]+1}
    }
  }
 
  k = cbind(rese,rest)
  return(k)
}


### Function for simulation studies ###



JL_simulation = function(starter_dose, pe_true, pt_true, pe_prior, pt_prior, posterior_function_E, posterior_function_T, cohorts, cohort_size, burnin ,nsteps, nsims, Ce, Ct){
  #First Round with starter dose
  count = 0
  dose_choices = c()
  assigned = c()
  
  while (count<nsims) {
    #make empty data
    cohort = cohorts
    ye = list(ye1 = c(), ye2 = c(), ye3 = c(), ye4 = c(),ye5 = c())
    yt = list(yt1 = c(), yt2 = c(), yt3 = c(), yt4 = c(), yt5 = c())
    
    
    tested_levels = c()
    tested_levels = list.append(tested_levels, starter_dose)
    
    #Establish starter values for the MCMC
    Be_starter = c()
    Be_starter[1] = pe_prior[1]
    for (i in 2:length(pe_prior)){
      Be_starter[i] = (pe_prior[i] - pe_prior[i-1])/(1- pe_prior[i-1])
    }
    
    Bt_starter = c()
    Bt_starter[1] = pt_prior[1]
    for (i in 2:length(pt_prior)){
      Bt_starter[i] = (pt_prior[i] - pt_prior[i-1])/(1- pt_prior[i-1])
    }
    theta_E_start = log(Be_starter/(1-Be_starter))
    theta_T_start = log(Bt_starter/(1-Bt_starter))   
    next_dose = starter_dose
    
    
    
    #Start loop for continuing trial for the inputted number of cohorts
    
    while (cohort > 0) {
      #SAMPLE NEXT COHORT
      current_dose = next_dose
      tested_levels = list.append(tested_levels, current_dose)
      new_sample = binomial_sampler(pe_true, pt_true, dose_level = current_dose, cohort_size)
      ye[[current_dose]] = list.append(ye[[current_dose]], new_sample[,1])
      yt[[current_dose]] = list.append(yt[[current_dose]], new_sample[,2])
      
      #MCMC to simulate Be and Bt
      theta_samples = run_MCMC(theta_E_start, theta_T_start, f_E = posterior_function_E, f_T = posterior_function_T , randomvalue = randomvalue_generator, nsteps = nsteps, ye=ye, yt=yt)

      #Transform Thetas into betas
      Beta_samples = matrix(c(exp(theta_samples[,1])/(1+exp(theta_samples[,1])), exp(theta_samples[,2])/(1+exp(theta_samples[,2])), exp(theta_samples[,3])/(1+exp(theta_samples[,3])), exp(theta_samples[,4])/(1+exp(theta_samples[,4])), exp(theta_samples[,5])/(1+exp(theta_samples[,5])), exp(theta_samples[,6])/(1+exp(theta_samples[,6])), exp(theta_samples[,7])/(1+exp(theta_samples[,7])), exp(theta_samples[,8])/(1+exp(theta_samples[,8])), exp(theta_samples[,9])/(1+exp(theta_samples[,9])), exp(theta_samples[,10])/(1+exp(theta_samples[,10]))), ncol = 10)
      Beta_samples = Beta_samples[burnin:nsteps,]
      
      #Get pe and pt from Be and Bt
      pe_sims = Be_to_pe(Beta_samples[,1:5])
      pt_sims = Bt_to_pt(Beta_samples[,6:10])
      
      #STOPPING RULES
      pe_counts = c(0,0,0,0,0)
      pt_counts = c(0,0,0,0,0)
      
      for (j in 1:ncol(pe_sims)){
        for (i in 1:nrow(pe_sims)){
          if (pe_sims[i,j] < phi_e){
            pe_counts[j] = pe_counts[j] +1
          }
          if (pt_sims[i,j] > phi_t){
            pt_counts[j] = pt_counts[j]+1
          }
        }
      }
      pe_stop = pe_counts/nrow(pe_sims)
      pt_stop = pt_counts/nrow(pt_sims)

      
      #Calculate Utilities from Pe and Pt simulated
      final_utilities = matrix(NA, nrow = nrow(pe_sims), ncol = ncol(pt_sims))
      mean_util = c()
      for (i in 1:ncol(pe_sims)){
        for (j in 1:nrow(pe_sims)){
          final_utilities[j,i] = Utility(pe_sims[j,i], pt_sims[j,i])
        }
        mean_util[i] = mean(final_utilities[,i])
      }
      
      
      #Next_dose is the the dose level that had the most instances of highest utility, make sure that it only chooses a dose level if the previous level has already been tested
      #check for dose levels that do not fit equation 2.8/2.9
      
      #Stop
      if (min(pe_stop)>Ce|min(pt_stop)>Ct){
        final_dose = 0
        break
      }
      
      #After dropping levels that do not fit in equation 2.8/2.9 update max_util_dose
      #sample_mean_utils is for the dose choice process, picking the next dose with probability equal to the utility proportion
      sample_mean_utils = 1.42+mean_util
      max_tested = max(tested_levels)
      for (i in 1:length(mean_util)){
        if (i > max_tested +1){
          mean_util[i] = NA
        }
      }

      max_util_dose = which.max(mean_util)

      
      #creating dose sampling group for picking next dose
      if (max_util_dose  == 1){
        dose_sample = c(1, 2)
      }else {if (max_util_dose == 5){
        dose_sample = c(4,5)
      }else {
        dose_sample = c(max_util_dose -1, max_util_dose, max_util_dose +1)
        }
      }
      
      dose_sample = dose_sample[dose_sample - max_tested<=1]
      
      if (length(dose_sample) == 0){
        dose_sample = c(max(tested_levels)+1)
        
      }
      
      dose_sample_counts = c()
      for (i in 1:length(dose_sample)){
        dose_sample_counts[i] = sample_mean_utils[dose_sample[i]]
      }

      if (max(dose_sample_counts == 0)){
        for (i in 1:length(dose_sample_counts)){
          dose_sample_counts[i] = 1
        }
      }
      if (length(dose_sample) == 1){
        next_dose = dose_sample
      }else{
        next_dose = sample(dose_sample, size = 1, prob = dose_sample_counts/sum(dose_sample_counts))
      }
      
      assigned = list.append(assigned, next_dose)
      cohort = cohort - 1
      
    }
    
    #After all data has been collected, the dose of highest utility is picked amongst tested dose levels
    best_dose_selection = unique(tested_levels)
    for (item in best_dose_selection){
      if (pe_stop[item]>Ce|pt_stop[item]>Ct){
        best_dose_selection = best_dose_selection[best_dose_selection != item]
      }
    }
    
    if (min(pe_stop)>Ce|min(pt_stop)>Ct){
      best_dose_selection = c()
    }
    
    last_max_col = c()
    
    if (length(best_dose_selection) == 0){
      final_dose = 0
    }else{
      for (i in 1:5){
        if (i %in% best_dose_selection){
          last_max_col[i] = mean_util[i]
        }else{
          last_max_col[i] = NA
        }
      }
      final_dose = which.max(last_max_col)
    }
    
    
    if (is.na(final_dose)){
      final_dose = 0
    }
    
    #dose_choices keeps track of the dose with the highest utility at the end of each simulation
    cat("\n\nbest dose chosen", final_dose , "\n")
    pe_empirical1= c()
    pt_empirical1=c()
    util_emp = c()
    for (i in 1:5){
      pe_empirical1[i] = sum(ye[[i]])/(length(ye[[i]]))
      pt_empirical1[i]  = sum(yt[[i]])/(length(yt[[i]]))
      pe_empirical1[is.na(pe_empirical1)] <- 0
      pt_empirical1[is.na(pt_empirical1)] <- 0
      util_emp[i] = Utility(pe_empirical1[i], pt_empirical1[i])
    }
    print(pe_stop)
    print(pt_stop)
    
    cat("\n pe values: ", colMeans(pe_sims))
    cat("\n pt values: ", colMeans(pt_sims))
    cat("\nmean util: ", mean_util)
    cat("\nemp_util:", util_emp, "\n\n\n")
    dose_choices = list.append(dose_choices, final_dose)
    count = count+1
    if (count%%10 == 0){print(count)}
  }
  
  #output the frequency each dose-level was selcted as the optimal dose and number of doses administered to patients
  end_util = data.frame(table(factor(dose_choices, levels = 0:5)))
  colnames(end_util) = c("dose_level", "Final_choice")
  
  assigned = data.frame(table(factor(assigned, levels = 0:5)))
  assigned[,2] = assigned[,2]/sum(assigned[,2])
  colnames(assigned) = c("dose_level", "assigned")
  
  end_util$Final_choice = end_util$Final_choice/nsims
  end_util$assigned = assigned$assigned
  
  end_util
  
}















