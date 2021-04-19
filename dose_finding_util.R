library(survey)
library(plyr)
library(tidyverse)
library(rlist)

### Turn Beta_kj's into p_kj form ###
Be_to_pe = function(Be){
  pe_sim = matrix(NA, nrow = length(Be[,1]), ncol = 5)
  pe_sim[,1] = Be[,1]
  for (i in 2:length(Be[1,])){
    pe_sim[,i] = pe_sim[,i-1] + (1 - pe_sim[, i-1])*Be[,i]
  }
  pe_sim
}

Bt_to_pt = function(Bt){
  pt_sim = matrix(NA, nrow = length(Bt[,1]), ncol = 5)
  pt_sim[,1] = Bt[,1]
  for (i in 2:length(Bt[1,])){
    pt_sim[,i] = pt_sim[,i-1] + (1 - pt_sim[, i-1])*Bt[,i]
  }
  pt_sim
}




### Utility Function ###

Utility = function(pe, pt){
  util = c()
  w1 = 0.33 #33% peanlty for toxicities less than threshold
  w2 = 1.09 #142% peanlty for toxicities greater than threshold
  phi_t = 0.3
  for (i in 1:length(pe)){
    if (pt[i] >= phi_t){
      util[i] = pe[i] - w1*pt[i]  - w2*pt[i]
    }
    else {
      util[i] = pe[i] - w1*pt[i]
    }
  }
  util
}


#Mode funciton ###
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



### Gumbel Distribution Data sampling functions ###

prob_sampler = function(pe, pt, dose_level, gam =3){
  pi_ab = c()
  a = c(1,1,0,0)
  b = c(1,0,1,0)
  j = dose_level
  #Row 1 = 1,1    Row 2 = 1,0(efficacy but no toxicity)   Row 3 = 0,1(toxicity but no efficacy)  Row 4 = 0,0
  for (i in 1:length(a)){
    pi_ab[i] = pe[j]^a[i] * (1-pe[j])^(1-a[i]) * pt[j]^b[i] * (1-pt[j])^(1-b[i]) + (-1)^(a[i]+b[i]) * pe[j] * (1-pe[j]) * pt[j] * (1-pt[j]) * ((exp(gam) - 1)/(exp(gam)+1))
  }
  pi_ab
}

binomial_sampler = function(pe_true, pt_true, dose_level, cohort_size){
  pi_ab = prob_sampler(pe_true, pt_true, dose_level, gam = 3)
  binomial_samples = replicate(cohort_size, sample(c(1,2,3,4), size = 1, replace = FALSE, prob = pi_ab) )
  samples_final = matrix(NA, nrow = cohort_size, ncol = 2)
  for (i in 1:length(binomial_samples)){
    if (binomial_samples[i] == 1) {samples_final[i,] = c(1,1)}
    if (binomial_samples[i] == 2) {samples_final[i,] = c(1,0)}
    if (binomial_samples[i] == 3) {samples_final[i,] = c(0,1)}
    if (binomial_samples[i] == 4) {samples_final[i,] = c(0,0)}
  }
  samples_final
}