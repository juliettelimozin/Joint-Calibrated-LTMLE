library(modelr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("simulate_MSM_treatment_switch.R")
set.seed(NULL)
library(MASS)
library(TrialEmulation)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)

treat <- c(-1,0,1)
outcome_prev <- c(-4.7,-3.8,-3)

scenarios <- tidyr::crossing(treat)

true_HR <- array(,dim = c(20,3,3))
true_MRD <- array(,dim = c(5,2,3,3))

for (l in 1){
  for (j in 3){
    simdata_censored<-DATA_GEN_treatment_switch(1000000, 5, 
                                                treat_prev = as.numeric(scenarios[l,1]),
                                                outcome_prev = as.numeric(outcome_prev[j]),
                                                censor = F)
    simdata_censored$C <- 0.0
    PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                                eligible ='eligible', cense = 'C',
                                                switch_d_cov = ~X1 + X2 + X3,
                                                outcome_cov = ~X1 + X2 + X3, model_var = c('assigned_treatment'),
                                                estimand_type = 'PP', quiet = T,
                                                save_weight_models = F,
                                                data_dir = getwd())
    switch_data <- PP_prep$data %>% 
      dplyr::filter(trial_period == 0) %>% 
      dplyr::mutate(t_1 = ifelse(followup_time == 1,1,0),
                    t_2 = ifelse(followup_time == 2,1,0),
                    t_3 = ifelse(followup_time == 3,1,0),
                    t_4 = ifelse(followup_time == 4,1,0),
                    t_1A = t_1*assigned_treatment,
                    t_2A = t_2*assigned_treatment,
                    t_3A = t_3*assigned_treatment,
                    t_4A = t_4*assigned_treatment,
                    t_1X1 = t_1*X1,
                    t_2X1 = t_2*X1,
                    t_3X1 = t_3*X1,
                    t_4X1 = t_4*X1,
                    t_1X2 = t_1*X2,
                    t_2X2 = t_2*X2,
                    t_3X2 = t_3*X2,
                    t_4X2 = t_4*X2,
                    t_1X3 = t_1*X3,
                    t_2X3 = t_2*X3,
                    t_3X3 = t_3*X3,
                    t_4X3 = t_4*X3)
    PP <- TrialEmulation::trial_msm(data = switch_data,
                                    outcome_cov = ~ X1 + X2+ X3+ assigned_treatment+
                                      t_1 + t_2 + t_3 + t_4 +
                                      t_1A + t_2A + t_3A + t_4A +
                                      t_1X1 + t_2X1 + t_3X1 + t_4X1 +
                                      t_1X2 + t_2X2 + t_3X2 + t_4X2 +
                                      t_1X3 + t_2X3 + t_3X3 + t_4X3,
                                    model_var = c('assigned_treatment'),
                                    glm_function = 'glm',
                                    include_trial_period = ~1, include_followup_time = ~1,
                                    estimand_type = 'PP',, quiet = T, use_sample_weights =  F)
    true_HR[,l,j] <- PP$model$coefficients
    
    simdata_censored_treat<-DATA_GEN_treatment_switch(1000000, 5, 
                                                      treat_prev = as.numeric(scenarios[l,1]),
                                                      outcome_prev = outcome_prev[j],
                                                      all_treat = T,
                                                      censor = F)
    simdata_censored_control<-DATA_GEN_treatment_switch(1000000,5, 
                                                        treat_prev = as.numeric(scenarios[l,1]),
                                                        outcome_prev = outcome_prev[j],
                                                        all_control = T,
                                                        censor = F)
    
    surv_data_treat <- simdata_censored_treat[ !duplicated(simdata_censored_treat[, c("ID")], fromLast=T),] %>% 
      dplyr::mutate(status = Y) %>% 
      dplyr::select(ID, t, status)
    
    f1 <- survfit(Surv(t, status) ~ 1, data = surv_data_treat)
    
    surv_data_control <- simdata_censored_control[ !duplicated(simdata_censored_control[, c("ID")], fromLast=T),] %>% 
      dplyr::mutate(status = Y) %>% 
      dplyr::select(ID, t, status)
    
    f2 <- survfit(Surv(t, status) ~ 1, data = surv_data_control)
    
    true_MRD[,1,l,j] <- f2$surv
    
    true_MRD[,2,l,j] <- f1$surv
    
  }
  save(true_HR, file = "Simulation results/true_HR_singletrial.rda")
  save(true_MRD, file = "Simulation results/true_MRD_singletrial.rda")
}


