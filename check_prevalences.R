library(modelr)
library(tidyverse)
library(tidyr)
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("simulate_MSM_treatment_switch.R")
set.seed(NULL)
library(MASS)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)


treat_prev <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)
sample_size <- c(500,1000)
alpha_y <- c(-4.7,-3.8,-3)
summary_table <- array(,dim = c(270,9))

scenarios <- tidyr::crossing(sample_size,conf,treat_prev)
for (j in 1:10){
for (i in 1:27){
  data <- DATA_GEN_treatment_switch(as.numeric(scenarios[i,1]), 5, conf = as.numeric(scenarios[i,2]), treat_prev = as.numeric(scenarios[i,3]) , 
                                      outcome_prev = -3.8, censor = F)
  
  initiators <- data %>% 
    filter(t == 0) %>% 
    summarise(CA = sum(A))
  events <- data %>% 
    summarise(CY = sum(Y))

  PP_prep <- TrialEmulation::data_preparation(data, id='ID', period='t', treatment='A', outcome='Y', 
                                              eligible ='eligible',
                                              switch_d_cov = ~X1 + X2 + X3,
                                              outcome_cov = ~X1 + X2 , model_var = c('assigned_treatment'),
                                              use_weight=T, use_censor=T, quiet = T,
                                              save_weight_models = F,
                                              data_dir = getwd())
  data <- data %>% 
    mutate(Z1 = X1^3/9,Z2 = X1*X2, Z3 = log(abs(X3))+4)

  
  PP_prep_miss <- TrialEmulation::data_preparation(data, id='ID', period='t', treatment='A', outcome='Y', 
                                              eligible ='eligible',
                                              switch_d_cov = ~Z1 + Z2 + Z3 ,
                                              outcome_cov = ~Z1 + Z2 + Z3 , model_var = c('assigned_treatment'),
                                              use_weight=T, use_censor=T, quiet = T,
                                              save_weight_models = F,
                                              data_dir = getwd())
  
  summary_table[i,] <- do.call(rbind, 
                               lapply(c(scenarios[i,1], scenarios[i,2],scenarios[i,3], 
                                        initiators$CA/scenarios[i,1],events$CY/scenarios[i,1],
                                        sd(PP_prep$data$weight),sd(PP_prep_miss$data$weight), 
                                        max(PP_prep$data$weight), max(PP_prep_miss$data$weight)), as.numeric))
}
}
colnames(summary_table) <- c('Sample size', 'Conf', 'treat_prev', 'Initiators', 'Events', 'SD(weights)', 'SD(weights_miss)', 'max(weights)', 'max(weightss_miss)')


for (i in 1:27){
  data <- DATA_GEN_treatment_switch(as.numeric(scenarios[i,1]), 5, conf = as.numeric(scenarios[i,2]), treat_prev = 1 , 
                                    outcome_prev = as.numeric(scenarios[i,3]), all_treat = FALSE, 
                                    all_control = FALSE, censor = F)
  
  initiators <- data %>% 
    filter(t == 0) %>% 
    summarise(CA = sum(A))
  events <- data %>% 
    summarise(CY = sum(Y))
  
  summary_table[i,] <- do.call(rbind, lapply(c(scenarios[i,1], scenarios[i,2],scenarios[i,3], initiators$CA/scenarios[i,1],events$CY/scenarios[i,1]), as.numeric))
}
