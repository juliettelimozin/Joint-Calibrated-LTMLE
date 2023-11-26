#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("simulate_MSM_treatment_switch.R")
library(TrialEmulation)
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)
library(nleqslv)
source('calibration_func_trials.R')
set.seed(NULL)

iters <- 200
l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
size <- c(500,1000,5000)
treat <- c(-1,0,1)

scenarios <- tidyr::crossing(size, treat)
objectives <- array(,dim = c(4,8,iters))
hr_estimates <- array(,dim = c(4,iters))
mr_estimates <- array(,dim = c(4,2,5,iters))
meandiffs <- array(, dim = c(5,4,6,iters))
max_weight <- array(,dim = c(4,iters))
# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 67)

for (i in 1:iters){
  tryCatch({
    simdata_censored<-DATA_GEN_treatment_switch(as.numeric(scenarios[l,1]), 5, 
                                                treat_prev = as.numeric(scenarios[l,2]),
                                                outcome_prev = -3.8,
                                                censor = F)
    
    ######### Correctly specified ############### 
    PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                                eligible ='eligible',
                                                switch_d_cov = ~X1 + X2 + X3,
                                                outcome_cov = ~X1 + X2 + X3, model_var = c('assigned_treatment'),
                                                use_weight=T, use_censor=T, quiet = T,
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
    
    
    data_restric <- simdata_censored %>% 
      dplyr::group_by(ID) %>% 
      dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1), nextX1 = lead(X1),nextX3 = lead(X3)) %>% 
      dplyr::mutate(CRA = cumsum(RA),
                    nextX1 = ifelse(is.na(nextX1), 0, nextX1),
                    nextX3 = ifelse(is.na(nextX3), 0, nextX3),
                    RA = ifelse(CRA == t+1,1,0),
                    RC = ifelse(lag(C) == 0,1,0),
                    A1X2 = A_0*X2,
                    A0X2 = (1-A_0)*X2,
                    A1nextX1 = A_0*nextX1,
                    A0nextX1 = (1-A_0)*nextX1,
                    A1nextX3= A_0*nextX3,
                    A0nextX3 = (1-A_0)*nextX3,
                    A1 = A_0,
                    A0 = 1-A_0,
                    sub = ID,
                    tall = t) %>% 
      merge(dplyr::select(switch_data,id, followup_time, weight), 
            by.x = c('ID', 't'), by.y = c('id', 'followup_time'), all.x = T) %>% 
      dplyr::mutate(weights = ifelse(!is.na(weight), weight,0)) %>% 
      dplyr::arrange(ID, t) 
    
    simdatafinal <- calibration(simdatafinal = data_restric, 
                                var = c('A1','A1nextX1', 'A1X2','A1nextX3',
                                        'A0','A0nextX1','A0X2','A0nextX3'))
    
    objectives[1,,i] <- simdatafinal$objective.IPW
    objectives[2,,i] <- simdatafinal$objective.Cali
    
    meandiffs_summary <- simdatafinal$data %>% 
      dplyr::mutate(RAX2 = RA*X2,
                    RAX3 = RA*X3,
                    RAX1 = RA*X1) %>% 
      dplyr::select(t,A1,A0,X2,X3,X1, RAX2, RAX3, RAX1, weights, Cweights)
    
    for (h in 1:4){
      try({
        meandiffs[1,h,1:3,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 1 ,c('RAX1', 'RAX2', 'RAX3')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 1 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
      try({
        meandiffs[2,h,1:3,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 1 ,c('RAX1', 'RAX2', 'RAX3')]*meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 1,c('weights')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 1 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
      try({
        meandiffs[3,h,1:3,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 1 ,c('RAX1', 'RAX2', 'RAX3')]*meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 1,c('Cweights')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 1 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
      
      try({
        meandiffs[1,h,4:6,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 0 ,c('RAX1', 'RAX2', 'RAX3')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 0 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
      try({
        meandiffs[2,h,4:6,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 0 ,c('RAX1', 'RAX2', 'RAX3')]*meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 0,c('weights')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 0 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
      try({
        meandiffs[3,h,4:6,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 0 ,c('RAX1', 'RAX2', 'RAX3')]*meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 0,c('Cweights')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 0 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
    }
    
    
    
    PP <- TrialEmulation::trial_msm(data = switch_data,
                                    outcome_cov = ~ X1 + X2+ X3 + assigned_treatment+
                                      t_1 + t_2 + t_3 + t_4 +
                                      t_1A + t_2A + t_3A + t_4A +
                                      t_1X1 + t_2X1 + t_3X1 + t_4X1 +
                                      t_1X2 + t_2X2 + t_3X2 + t_4X2 +
                                      t_1X3 + t_2X3 + t_3X3 + t_4X3,
                                    model_var = c('assigned_treatment'),
                                    glm_function = 'glm',
                                    include_trial_period = ~1, include_followup_time = ~1,
                                    use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
    hr_estimates[1,i] <- PP$model$coefficients[2]
    
    max_weight[1,i] <- max(switch_data$weight)
    
    switch_data$weight <- simdatafinal$data %>% dplyr::filter(RA == 1) %>% dplyr::select(Cweights)
    
    max_weight[2,i] <- max(switch_data$weight)
    
    PP_calibrated <- TrialEmulation::trial_msm(data = switch_data,
                                               outcome_cov = ~ X1 + X2+ + X3 + assigned_treatment+
                                                 t_1 + t_2 + t_3 + t_4 +
                                                 t_1A + t_2A + t_3A + t_4A +
                                                 t_1X1 + t_2X1 + t_3X1 + t_4X1 +
                                                 t_1X2 + t_2X2 + t_3X2 + t_4X2 +
                                                 t_1X3 + t_2X3 + t_3X3 + t_4X3,
                                               model_var = c('assigned_treatment'),
                                               glm_function = 'glm',
                                               include_trial_period = ~1, include_followup_time = ~1,
                                               use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
    hr_estimates[2,i] <- PP_calibrated$model$coefficients[2]
    
    design_mat <- expand.grid(id = 1:as.numeric(dim(switch_data)[1]),
                              trial_period = 0:4,
                              followup_time = 0:4)
    design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]
    
    fitting_data_treatment <-  switch_data %>%
      dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>%
      dplyr::select(id,trial_period, followup_time, X1,  X2, X3, assigned_treatment) %>%
      merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>%
      dplyr::group_by(id) %>%
      tidyr::fill( X1,X2,X3, assigned_treatment,.direction = "down") %>%
      dplyr::ungroup() %>%
      dplyr::select(id, trial_period, followup_time, X1, X2, X3, assigned_treatment) %>%
      merge(data.frame(id = switch_data$id, trial_period = switch_data$trial_period), by = c("id", "trial_period"), all.y = TRUE) %>%
      dplyr::arrange(id, trial_period, followup_time) %>%
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
                    t_4X3 = t_4*X3) %>%
      dplyr::filter(trial_period == 0)
    
    fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]
    
    fitting_data_control <- fitting_data_treatment %>%
      dplyr::mutate(assigned_treatment = assigned_treatment*0,
                    t_1A = t_1*0,
                    t_2A = t_2*0,
                    t_3A = t_3*0,
                    t_4A = t_4*0)
    
    Y_pred_PP_treatment <- predict.glm(PP$model,
                                       fitting_data_treatment,
                                       type = "response")
    Y_pred_PP_control <- predict.glm(PP$model,
                                     fitting_data_control,
                                     type = "response")
    Y_pred_PP_treatment_cali <- predict.glm(PP_calibrated$model,
                                            fitting_data_treatment,
                                            type = "response")
    Y_pred_PP_control_cali <- predict.glm(PP_calibrated$model,
                                          fitting_data_control,
                                          type = "response")
    predicted_probas_PP <- fitting_data_treatment %>%
      dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment,
                    predicted_proba_control = Y_pred_PP_control,
                    predicted_proba_treatment_cali = Y_pred_PP_treatment_cali,
                    predicted_proba_control_cali = Y_pred_PP_control_cali) %>%
      dplyr::group_by(id, trial_period) %>%
      dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                    cum_hazard_control = cumprod(1-predicted_proba_control),
                    cum_hazard_treatment_cali = cumprod(1-predicted_proba_treatment_cali),
                    cum_hazard_control_cali = cumprod(1-predicted_proba_control_cali)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(followup_time) %>%
      dplyr::summarise(risk_treatment = mean(cum_hazard_treatment),
                       risk_control = mean(cum_hazard_control),
                       mrd = risk_control - risk_treatment,
                       risk_treatment_cali = mean(cum_hazard_treatment_cali),
                       risk_control_cali = mean(cum_hazard_control_cali),
                       mrd_cali = risk_control_cali - risk_treatment_cali)
    
    mr_estimates[1,1,,i] <- pull(predicted_probas_PP,risk_control)
    mr_estimates[1,2,,i] <- pull(predicted_probas_PP,risk_treatment)
    mr_estimates[2,1,,i] <- pull(predicted_probas_PP,risk_control_cali)
    mr_estimates[2,2,,i] <- pull(predicted_probas_PP,risk_treatment_cali)
    
    
    #########  Misspecified ############### 
    simdata_censored <- simdata_censored %>% 
      mutate(Z1 = X1^3/9,Z2 = X1*X2, Z3 = log(abs(X3))+4)
    PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                                eligible ='eligible',
                                                switch_d_cov = ~ Z1 + Z2 + Z3,
                                                outcome_cov = ~Z1 + Z2 + Z3, model_var = c('assigned_treatment'),
                                                use_weight=T, use_censor=T, quiet = T,
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
                    t_1Z1 = t_1*Z1,
                    t_2Z1 = t_2*Z1,
                    t_3Z1 = t_3*Z1,
                    t_4Z1 = t_4*Z1,
                    t_1Z2 = t_1*Z2,
                    t_2Z2 = t_2*Z2,
                    t_3Z2 = t_3*Z2,
                    t_4Z2 = t_4*Z2,
                    t_1Z3 = t_1*Z3,
                    t_2Z3 = t_2*Z3,
                    t_3Z3 = t_3*Z3,
                    t_4Z3 = t_4*Z3)
    
    
    data_restric <- simdata_censored %>% 
      dplyr::group_by(ID) %>% 
      dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1), 
                    nextZ1 = lead(Z1),nextZ2 = lead(Z2),nextZ3 = lead(Z3)) %>% 
      dplyr::mutate(CRA = cumsum(RA),
                    nextZ1 = ifelse(is.na(nextZ1), 0, nextZ1),
                    nextZ2 = ifelse(is.na(nextZ2), 0, nextZ2),
                    nextZ3 = ifelse(is.na(nextZ3), 0, nextZ3),
                    RA = ifelse(CRA == t+1,1,0),
                    RC = ifelse(lag(C) == 0,1,0),
                    A1nextZ1 = A_0*exp(nextZ1),
                    A0nextZ1 = (1-A_0)*exp(nextZ1),
                    A1nextZ2 = A_0*exp(nextZ2),
                    A0nextZ2 = (1-A_0)*exp(nextZ2),
                    A1nextZ3 = A_0*exp(nextZ3),
                    A0nextZ3 = (1-A_0)*exp(nextZ3),
                    A1 = A_0,
                    A0 = 1-A_0,
                    sub = ID,
                    tall = t) %>% 
      merge(dplyr::select(switch_data,id, followup_time, weight), 
            by.x = c('ID', 't'), by.y = c('id', 'followup_time')) %>% 
      dplyr::mutate(weights = weight) %>% 
      dplyr::arrange(ID, t) 
    
    simdatafinal <- calibration(simdatafinal = data_restric, 
                                var = c('A1', 'A1nextZ1', 'A1nextZ2','A1nextZ3',
                                        'A0','A0nextZ1', 'A0nextZ2','A0nextZ3'))
    objectives[3,,i] <- simdatafinal$objective.IPW
    objectives[4,,i] <- simdatafinal$objective.Cali
    
    meandiffs_summary <- simdatafinal$data %>%  
      dplyr::mutate(RAX2 = RA*X2,
                    RAX3 = RA*X3,
                    RAX1 = RA*X1) %>% 
      dplyr::select(t,A1,A0,X2,X3,X1, RAX2, RAX3, RAX1, weights, Cweights)
    
    for (h in 1:4){
      try({
        meandiffs[4,h,1:3,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 1 ,c('RAX1', 'RAX2', 'RAX3')]*meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 1,c('weights')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 1 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
      try({
        meandiffs[5,h,1:3,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 1 ,c('RAX1', 'RAX2', 'RAX3')]*meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 1,c('Cweights')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 1 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
      try({
        meandiffs[4,h,4:6,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 0 ,c('RAX1', 'RAX2', 'RAX3')]*meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 0,c('weights')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 0 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
      try({
        meandiffs[5,h,4:6,i] <-colMeans(meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 0 ,c('RAX1', 'RAX2', 'RAX3')]*meandiffs_summary[meandiffs_summary$t == h & meandiffs_summary$A1 == 0,c('Cweights')]) - colMeans(meandiffs_summary[meandiffs_summary$t == 1 & meandiffs_summary$A1 == 0 ,c('X1', 'X2', 'X3')])
      },
      silent = T)
    }
    
    
    PP <- TrialEmulation::trial_msm(data = switch_data,
                                    outcome_cov = ~ Z1 + Z2 + Z3+ assigned_treatment+
                                      t_1 + t_2 + t_3 + t_4 +
                                      t_1A + t_2A + t_3A + t_4A +
                                      t_1Z1 + t_2Z1 + t_3Z1 + t_4Z1 +
                                      t_1Z2 + t_2Z2 + t_3Z2 + t_4Z2+
                                      t_1Z3 + t_2Z3 + t_3Z3 + t_4Z3,
                                    model_var = c('assigned_treatment'),
                                    glm_function = 'glm',
                                    include_trial_period = ~1, include_followup_time = ~1,
                                    use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
    hr_estimates[3,i] <- PP$model$coefficients[2]
    max_weight[3,i] <- max(switch_data$weight)
    switch_data$weight <- simdatafinal$data %>% dplyr::filter(RA == 1) %>% dplyr::select(Cweights)
    max_weight[4,i] <- max(switch_data$weight)
    PP_calibrated <- TrialEmulation::trial_msm(data = switch_data,
                                               outcome_cov = ~ Z1 + Z2 + Z3+ assigned_treatment+
                                                 t_1 + t_2 + t_3 + t_4 +
                                                 t_1A + t_2A + t_3A + t_4A +
                                                 t_1Z1 + t_2Z1 + t_3Z1 + t_4Z1 +
                                                 t_1Z2 + t_2Z2 + t_3Z2 + t_4Z2+
                                                 t_1Z3 + t_2Z3 + t_3Z3 + t_4Z3,
                                               model_var = c('assigned_treatment'),
                                               glm_function = 'glm',
                                               include_trial_period = ~1, include_followup_time = ~1,
                                               use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
    hr_estimates[4,i] <- PP_calibrated$model$coefficients[2]
    
    design_mat <- expand.grid(id = 1:as.numeric(dim(switch_data)[1]),
                              trial_period = 0:4,
                              followup_time = 0:4)
    design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]
    
    fitting_data_treatment <-  switch_data %>%
      dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>%
      dplyr::select(id,trial_period, followup_time, Z1,  Z2, Z3, assigned_treatment) %>%
      merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>%
      dplyr::group_by(id) %>%
      tidyr::fill( Z1,Z2,Z3,assigned_treatment,.direction = "down") %>%
      dplyr::ungroup() %>%
      dplyr::select(id, trial_period, followup_time, Z1, Z2,Z3, assigned_treatment) %>%
      merge(data.frame(id = switch_data$id, trial_period = switch_data$trial_period), by = c("id", "trial_period"), all.y = TRUE) %>%
      dplyr::arrange(id, trial_period, followup_time) %>%
      dplyr::mutate(t_1 = ifelse(followup_time == 1,1,0),
                    t_2 = ifelse(followup_time == 2,1,0),
                    t_3 = ifelse(followup_time == 3,1,0),
                    t_4 = ifelse(followup_time == 4,1,0),
                    t_1A = t_1*assigned_treatment,
                    t_2A = t_2*assigned_treatment,
                    t_3A = t_3*assigned_treatment,
                    t_4A = t_4*assigned_treatment,
                    t_1Z1 = t_1*Z1,
                    t_2Z1 = t_2*Z1,
                    t_3Z1 = t_3*Z1,
                    t_4Z1 = t_4*Z1,
                    t_1Z2 = t_1*Z2,
                    t_2Z2 = t_2*Z2,
                    t_3Z2 = t_3*Z2,
                    t_4Z2 = t_4*Z2,
                    t_1Z3 = t_1*Z3,
                    t_2Z3 = t_2*Z3,
                    t_3Z3 = t_3*Z3,
                    t_4Z3 = t_4*Z3) %>%
      dplyr::filter(trial_period == 0)
    
    fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]
    
    fitting_data_control <- fitting_data_treatment %>%
      dplyr::mutate(assigned_treatment = assigned_treatment*0,
                    t_1A = t_1*0,
                    t_2A = t_2*0,
                    t_3A = t_3*0,
                    t_4A = t_4*0)
    
    Y_pred_PP_treatment <- predict.glm(PP$model,
                                       fitting_data_treatment,
                                       type = "response")
    Y_pred_PP_control <- predict.glm(PP$model,
                                     fitting_data_control,
                                     type = "response")
    Y_pred_PP_treatment_cali <- predict.glm(PP_calibrated$model,
                                            fitting_data_treatment,
                                            type = "response")
    Y_pred_PP_control_cali <- predict.glm(PP_calibrated$model,
                                          fitting_data_control,
                                          type = "response")
    predicted_probas_PP <- fitting_data_treatment %>%
      dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment,
                    predicted_proba_control = Y_pred_PP_control,
                    predicted_proba_treatment_cali = Y_pred_PP_treatment_cali,
                    predicted_proba_control_cali = Y_pred_PP_control_cali) %>%
      dplyr::group_by(id, trial_period) %>%
      dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                    cum_hazard_control = cumprod(1-predicted_proba_control),
                    cum_hazard_treatment_cali = cumprod(1-predicted_proba_treatment_cali),
                    cum_hazard_control_cali = cumprod(1-predicted_proba_control_cali)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(followup_time) %>%
      dplyr::summarise(risk_treatment = mean(cum_hazard_treatment),
                       risk_control = mean(cum_hazard_control),
                       mrd = risk_control - risk_treatment,
                       risk_treatment_cali = mean(cum_hazard_treatment_cali),
                       risk_control_cali = mean(cum_hazard_control_cali),
                       mrd_cali = risk_control_cali - risk_treatment_cali)
    mr_estimates[3,1,,i] <- pull(predicted_probas_PP,risk_control)
    mr_estimates[3,2,,i] <- pull(predicted_probas_PP,risk_treatment)
    mr_estimates[4,1,,i] <- pull(predicted_probas_PP,risk_control_cali)
    mr_estimates[4,2,,i] <- pull(predicted_probas_PP,risk_treatment_cali)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
save(meandiffs, file = paste("Simulation results/meandiffs_singletrial_med_scneario2_",as.character(l),".rda", sep = ""))
save(objectives, file = paste("Simulation results/objectives_singletrial_med_scneario2_",as.character(l),".rda", sep = ""))
save(hr_estimates, file = paste("Simulation results/hr_estimates_singletrial_med_scneario2_",as.character(l),".rda", sep = ""))
save(mr_estimates, file = paste("Simulation results/mr_estimates_singletrial_med_scneario2_",as.character(l),".rda", sep = ""))
