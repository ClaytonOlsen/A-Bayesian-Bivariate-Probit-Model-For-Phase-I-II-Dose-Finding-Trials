

library(tmvtnorm)
library(rlist)
library(truncnorm)

#Conditional posteriors for alpha_e, alpha_t, y_tildE, y_tildT
#alpha E

alpha_E_Gibbs = function(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt, level){
  data_counts_e = c(length(ye[[1]]), length(ye[[2]]), length(ye[[3]]), length(ye[[4]]), length(ye[[5]]))
  if (level ==1){
    t_E = t_E1
  }else{
    t_E = t_E2
  }
  a_e_prior =  alpha_E_bar[level]/t_E^2
  sig_alpha_e = ( sum(data_counts_e[level:5])/(1-(2*p-1)^2) + 1/t_E^2 )^(-1)
  
  part_mean = c()
  for (j in level:5){
    if(length(ye[[j]]) == 0){
      next
    }
    for (i in 1:length(ye[[j]])){
      z_ei = y_tild_E[[j]][i] - (2*p-1)*(y_tild_T[[j]][i] - sum(alpha_T[1:j]))
      mu_alpha_e =  (z_ei - (sum(alpha_E[1:j])- alpha_E[level])/(1 - (2*p-1)^2))
      part_mean = list.append(part_mean, mu_alpha_e)
      
    }
  }
  alpha_mean_e = sig_alpha_e *(sum(part_mean) + a_e_prior)
  if (level > 1){
    return(rtruncnorm(1, a=0, b= Inf, alpha_mean_e, sig_alpha_e))
  }else{
    return(rnorm(1, alpha_mean_e, sig_alpha_e))
  }
}


#alpha T

alpha_T_Gibbs = function(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt, level){
  data_counts_t = c(length(yt[[1]]), length(yt[[2]]), length(yt[[3]]), length(yt[[4]]), length(yt[[5]]))
  if (level ==1){
    t_T = t_T1
  }else{
    t_T = t_T2
  }
  
  a_t_prior =  alpha_T_bar[level]/t_T^2
  sig_alpha_t = ( sum(data_counts_t[level:5])/(1-(2*p-1)^2) + 1/t_T^2 )^(-1)
  
  part_mean = c()
  for (j in level:5){
    if(length(yt[[j]]) ==0){
      next
    }
    for (i in 1:length(yt[[j]])){
      z_ti = y_tild_T[[j]][i] - (2*p-1)*(y_tild_E[[j]][i] - sum(alpha_E[1:j]))
      mu_alpha_t =  (z_ti - (sum(alpha_T[1:j])-alpha_T[level])/(1 - (2*p-1)^2))
      part_mean = list.append(part_mean, mu_alpha_t)
      
    }
  }
  alpha_mean_t = sig_alpha_t *(sum(part_mean) + a_t_prior)
  if (level > 1){
    return(rtruncnorm(1, a=0, b= Inf, alpha_mean_t, sig_alpha_t))
  }else{
    return(rnorm(1, alpha_mean_t, sig_alpha_t))
  }
}

#y_tilde_E
y_tild_E_Gibbs = function(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt){
  
  y_tild_E = list(ye1 = c(), ye2 = c(), ye3 = c(), ye4 = c(),ye5 = c())
  
  for (j in 1:5){
    if(length(ye[[j]]) ==0){
      next
    }
    for (i in 1:length(ye[[j]])){
      mu_ytild_E = sum(alpha_E[1:j]) + (2*p-1)*(y_tild_T[[j]][i] - sum(alpha_T[1:j]))
      sig_tild_E = 1 -(2*p -1)^2
      if (ye[[j]][i] == 1){
        y_tild_E[[j]][i] = rtruncnorm(1, a=0, b= Inf, mu_ytild_E, sig_tild_E)
      }else{
        y_tild_E[[j]][i] = rtruncnorm(1, a=-Inf, b= 0, mu_ytild_E, sig_tild_E)
      }
    }
  }
  return(y_tild_E)
}


#y_tilde_T
y_tild_T_Gibbs = function(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt){
  
  y_tild_T = list(yt1 = c(), yt2 = c(), yt3 = c(), yt4 = c(), yt5 = c())
  
  for (j in 1:5){
    if(length(yt[[j]]) ==0){
      next
    }
    for (i in 1:length(yt[[j]])){
      mu_ytild_T = sum(alpha_T[1:j]) + (2*p-1)*(y_tild_E[[j]][i] - sum(alpha_E[1:j]))
      sig_tild_T = 1 -(2*p -1)^2
      if (yt[[j]][i] == 1){
        y_tild_T[[j]][i] = rtruncnorm(1, a=0, b= Inf, mu_ytild_T, sig_tild_T)
      }else{
        y_tild_T[[j]][i] = rtruncnorm(1, a=-Inf, b= 0, mu_ytild_T, sig_tild_T)
      }
    }
  }
  return(y_tild_T)
}


#Posterior for rho and Metropolis Hastig step

p_posterior = function(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt){
  post_prod = c()
  for (j in 1:5){
    if (length(ye[[j]]) <1){
      next
    }
    for (i in 1:length(ye[[j]])){
      posterior_val = dmvnorm(c(y_tild_E[[j]][i], y_tild_T[[j]][i] ), c(sum(alpha_E[1:j]), sum(alpha_T[1:j])), matrix(c(1, 2*p-1, 2*p-1, 1), nrow = 2), log = TRUE)
      #print(posterior_val)
      post_prod = list.append(post_prod, posterior_val)
    }
  }
  sum(post_prod) * dbeta(p, ap, bp, log = TRUE)
}

randomvalue_generator <- function(x) {rtruncnorm(1, a=0, b=1, x, 0.08)}

#SINGLE BETA UPDATES
stepMH <- function(f, randomvalue, alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt) {
  ## Pick new point
  xp = randomvalue(p)
  old = f(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt)
  new = f(alpha_E, alpha_T, y_tild_E, y_tild_T, xp, ye, yt)
  ## Acceptance probability:
  que <- min(1, new- old)
  if (is.na(que)){
    que = 0
  }
  ## Accept new point with probability alpha:
  if (log(runif(1)) < que) {
    p <- xp
  }
  ## Returning the point:
  p
}



#GIBBS SAMPLER
gibbs_process = function(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt, nsteps){
  final_alphas_E = matrix(NA, nrow = nsteps, ncol = 5)
  final_alphas_T = matrix(NA, nrow = nsteps, ncol = 5)
  final_ye_tilde = matrix(NA, nrow = nsteps, ncol = 5)
  final_yt_tilde = matrix(NA, nrow = nsteps, ncol = 5)
  p_samples = c()
  
  #to get starter values
  
  #MCMC
  
  for (i in 1:nsteps){
    y_tild_E = y_tild_E_Gibbs(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt)
    y_tild_T = y_tild_T_Gibbs(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt)
    
    for (j in 1:5){
      alpha_E[j] = alpha_E_Gibbs(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt, j)
      alpha_T[j] = alpha_T_Gibbs(alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt, j)
    }
    p <- stepMH(p_posterior, randomvalue_generator, alpha_E, alpha_T, y_tild_E, y_tild_T, p, ye, yt)
    p_samples[i] = p
    final_alphas_E[i, ] = alpha_E
    final_alphas_T[i, ] = alpha_T
  }
  
  return(cbind(final_alphas_E, final_alphas_T, p_samples))
}




#FULL ALGORITHM


full_algorithm_probit = function(starter_dose, pe_true, pt_true, alpha_E_bar, alpha_T_bar, p_start, cohorts, cohort_size, burnin ,nsteps, nsims, Ce, Ct){
  #First Round with starter dose
  count = 0
  assigned = c()
  dose_choices = c()
  while (count<nsims) {
    #make empty data
    cohort = cohorts
    ye = list(ye1 = c(), ye2 = c(), ye3 = c(), ye4 = c(),ye5 = c())
    yt = list(yt1 = c(), yt2 = c(), yt3 = c(), yt4 = c(), yt5 = c())
    tested_levels = c(starter_dose)
    next_dose = starter_dose
    
    
    #Start loop for each cohort
    while (cohort > 0) {
      #SAMPLE NEXT COHORT
      current_dose = next_dose
      tested_levels = list.append(tested_levels, current_dose)
      
      new_sample = binomial_sampler(pe_true, pt_true, dose_level = current_dose, cohort_size)
      ye[[current_dose]] = list.append(ye[[current_dose]], new_sample[,1])
      yt[[current_dose]] = list.append(yt[[current_dose]], new_sample[,2])
      
      #Starter y_tilde for the gibbs sampler
      y_tild_E_start = ye
      y_tild_T_start = yt
      y_tild_E_start = y_tild_E_Gibbs(alpha_E_bar, alpha_T_bar, y_tild_E_start, y_tild_T_start, p=p_start, ye, yt )
      y_tild_T_start = y_tild_T_Gibbs(alpha_E_bar, alpha_T_bar, y_tild_E_start, y_tild_T_start, p=p_start, ye, yt )
      
      
      #Gibbs sampling
      gibbs_samps = gibbs_process(alpha_E_bar, alpha_T_bar, y_tild_E_start, y_tild_T_start, p_start, ye, yt, nsteps)
      gibbs_samps = gibbs_samps[burnin:nsteps, ]

      pe_sims = matrix(NA, nrow = nrow(gibbs_samps), ncol = 5)
      pt_sims = matrix(NA, nrow = nrow(gibbs_samps), ncol = 5)
      for (i in 1:nrow(gibbs_samps)){
        pe_sims[i,1] = 1 - pnorm(0, gibbs_samps[i,1], 1)
        pe_sims[i,2] = 1 - pnorm(0, gibbs_samps[i,2]+ gibbs_samps[i,1], 1)
        pe_sims[i,3] = 1 - pnorm(0, gibbs_samps[i,3]+gibbs_samps[i,2]+ gibbs_samps[i,1], 1) 
        pe_sims[i,4] = 1 - pnorm(0, gibbs_samps[i,4]+gibbs_samps[i,3]+gibbs_samps[i,2]+ gibbs_samps[i,1], 1)
        pe_sims[i,5] = 1 - pnorm(0, gibbs_samps[i,5]+gibbs_samps[i,4]+gibbs_samps[i,3]+gibbs_samps[i,2]+ gibbs_samps[i,1], 1)
        pt_sims[i,1] = 1 - pnorm(0, gibbs_samps[i,6], 1)
        pt_sims[i,2] = 1 - pnorm(0, gibbs_samps[i,7]+ gibbs_samps[i,6], 1)
        pt_sims[i,3] = 1 - pnorm(0, gibbs_samps[i,8]+gibbs_samps[i,7]+ gibbs_samps[i,6], 1) 
        pt_sims[i,4] = 1 - pnorm(0, gibbs_samps[i,9]+gibbs_samps[i,8]+gibbs_samps[i,7]+ gibbs_samps[i,6], 1)
        pt_sims[i,5] = 1 - pnorm(0, gibbs_samps[i,10]+gibbs_samps[i,9]+gibbs_samps[i,8]+gibbs_samps[i,7]+ gibbs_samps[i,6], 1)
      }

      
      #STOPPING RULES
      phi_e = 0.2
      phi_t = 0.3
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
      
      if (min(pe_stop)>Ce|min(pt_stop)>Ct){
        final_dose = 0
        break
      }
      
      #Utility Function
      final_utilities = matrix(NA, nrow = nrow(pe_sims), ncol = ncol(pt_sims))
      mean_util = c()
      for (i in 1:ncol(pe_sims)){
        for (j in 1:nrow(pe_sims)){
          final_utilities[j,i] = Utility(pe_sims[j,i], pt_sims[j,i])
        }
        mean_util[i] = mean(final_utilities[,i])
      }
      
      #After dropping levels that do not fit in equation 2.8/2.9 update max_util_dose
      sample_mean_utils = mean_util + 1.42
      
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
      
      #Pick next dose level
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
    
    #After all cohorts are tested, select the optimal dose level
    best_dose_selection = unique(tested_levels)
    
    
    for (item in best_dose_selection){
      if (pe_stop[item]>Ce|pt_stop[item]>Ct){
        best_dose_selection = best_dose_selection[best_dose_selection != item]
      }
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

  end_util = data.frame(table(factor(dose_choices, levels = 0:5)))
  colnames(end_util) = c("dose_level", "Final_choice")
  
  assigned = data.frame(table(factor(assigned, levels = 0:5)))
  assigned[,2] = assigned[,2]/sum(assigned[,2])
  colnames(assigned) = c("dose_level", "assigned")
  
  end_util$Final_choice = end_util$Final_choice/nsims
  end_util$assigned = assigned$assigned
  
  end_util
}



















