#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("simulate_MSM_treatment_switch.R")
library(TrialEmulation)
library(MASS)
library(sandwich)
library(foreach)
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
# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 2)
multiResultClass <- function(objectiveIPW = NULL,
                             objectiveCali = NULL,
                             objectiveIPWseq = NULL,
                             objectiveCaliseq = NULL,
                             objectiveIPWmis = NULL,
                             objectiveCalimis = NULL,
                             objectiveIPWseqmis = NULL,
                             objectiveCaliseqmis = NULL,
                             hr_estimates=NULL,mr_1_estimates = NULL, mr_0_estimates = NULL, balance_summary = NULL)
{
  me <- list(
    objectiveIPW = objectiveIPW,
    objectiveCali = objectiveCali,
    objectiveIPWseq = objectiveIPWseq,
    objectiveCaliseq = objectiveCaliseq,
    objectiveIPWmis = objectiveIPWmis,
    objectiveCalimis = objectiveCalimis,
    objectiveIPWseqmis = objectiveIPWseqmis,
    objectiveCaliseqmis = objectiveCaliseqmis,
    hr_estimates = hr_estimates,
    mr_1_estimates = mr_1_estimates,
    mr_0_estimates = mr_0_estimates,
    balance_summary = balance_summary
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}
oper <- foreach(i = 1:iters,.combine=cbind) %dopar% {
  tryCatch({
    result <- multiResultClass()
    simdata_censored<-DATA_GEN_treatment_switch(as.numeric(scenarios[l,1]), 5, 
                                                treat_prev = as.numeric(scenarios[l,2]),
                                                outcome_prev = -3,
                                                censor = F)
    simdata_censored$C <- 0.0
    ######### Correctly specified ############### 
    PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y',cense = 'C', 
                                                eligible ='eligible',
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
    
    
    data_restric <- simdata_censored %>% 
      dplyr::group_by(ID) %>% 
      dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1)) %>% 
      dplyr::mutate(CRA = cumsum(RA),
                    RA = ifelse(CRA == t+1,1,0),
                    RC = ifelse(lag(C) == 0,1,0),
                    A1X2 = A_0*X2,
                    A0X2 = (1-A_0)*X2,
                    A1X1 = A_0*X1,
                    A0X1 = (1-A_0)*X1,
                    A1X3 = A_0*X3,
                    A0X3 = (1-A_0)*X3,
                    A1 = A_0,
                    A0 = 1-A_0,
                    sub = ID,
                    tall = t) %>% 
      merge(dplyr::select(switch_data,id, followup_time, weight), 
            by.x = c('ID', 't'), by.y = c('id', 'followup_time'), all.x = T) %>% 
      dplyr::mutate(weights = ifelse(!is.na(weight), weight,0)) %>% 
      dplyr::arrange(ID, t) 
    
    ################### Calibration aggregated#######################
    simdatafinal1 <- calibration(simdatafinal = data_restric, 
                                        var = c('A1','A1X1', 'A1X2','A1X3',
                                                'A0','A0X1','A0X2','A0X3'))
    
    
    result$objectiveIPW <- simdatafinal1$objective.IPW
    result$objectiveCali <- simdatafinal1$objective.Cali
    
    ################### Calibration by time #######################
    simdatafinal2 <- calibration_by_time(simdatafinal = data_restric, 
                                        var = c('A1','A1X1', 'A1X2','A1X3',
                                                'A0','A0X1','A0X2','A0X3'))
    
  
    result$objectiveIPWseq <- simdatafinal2$objective.IPW
    result$objectiveCaliseq <- simdatafinal2$objective.Cali
    
    meandiffs_summary <- simdatafinal2$data %>% 
      dplyr::mutate(RAX2 = RA*X2,
                    RAX3 = RA*X3,
                    RAX1 = RA*X1,
                    RAA1 = RA*A1,
                    RAA0 = RA*A0,
                    Cweights_sequential = Cweights,
                    Cweights = simdatafinal1$data$Cweights) %>% 
      dplyr::select(t,RA,A1,A0,X2,X3,X1, RAX2, RAX3, RAX1,RAA1, RAA0, weights, Cweights,Cweights_sequential)
    
    treatment_numbers <- meandiffs_summary %>% 
      dplyr::group_by(t) %>% 
      dplyr::summarise(Treated_unadjusted = sum(RAA1),
                       Control_unadjusted = sum(RAA0),
                       Treated_IPW = sum(RAA1*weights),
                       Control_IPW = sum(RAA0*weights),
                       Treated_Cali = sum(RAA1 * Cweights),
                       Control_Cali = sum(RAA0*Cweights),
                       Treated_Cali_sequential = sum(RAA1 * Cweights_sequential),
                       Control_Cali_sequential = sum(RAA0*Cweights_sequential),
                       X1_treated = sum(RAX1*A1)/sum(RAA1),
                       X1_control = sum(RAX1*A0)/sum(RAA0),
                       X1_treated_IPW = sum(RAX1*weights*A1)/sum(RAA1*weights),
                       X1_control_IPW = sum(RAX1*weights*A0)/sum(RAA0*weights),
                       X1_treated_Cali_sequential = sum(RAX1*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X1_control_Cali_sequential = sum(RAX1*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X1_treated_Cali = sum(RAX1*Cweights*A1)/sum(RAA1*Cweights),
                       X1_control_Cali = sum(RAX1*Cweights*A0)/sum(RAA0*Cweights),
                       X2_treated = sum(RAX2*A1)/sum(RAA1),
                       X2_control = sum(RAX2*A0)/sum(RAA0),
                       X2_treated_IPW = sum(RAX2*weights*A1)/sum(RAA1*weights),
                       X2_control_IPW = sum(RAX2*weights*A0)/sum(RAA0*weights),
                       X2_treated_Cali_sequential = sum(RAX2*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X2_control_Cali_sequential = sum(RAX2*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X2_treated_Cali = sum(RAX2*Cweights*A1)/sum(RAA1*Cweights),
                       X2_control_Cali = sum(RAX2*Cweights*A0)/sum(RAA0*Cweights),
                       X3_treated = sum(RAX3*A1)/sum(RAA1),
                       X3_control = sum(RAX3*A0)/sum(RAA0),
                       X3_treated_IPW = sum(RAX3*weights*A1)/sum(RAA1*weights),
                       X3_control_IPW = sum(RAX3*weights*A0)/sum(RAA0*weights),
                       X3_treated_Cali_sequential = sum(RAX3*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X3_control_Cali_sequential = sum(RAX3*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X3_treated_Cali = sum(RAX3*Cweights*A1)/sum(RAA1*Cweights),
                       X3_control_Cali = sum(RAX3*Cweights*A0)/sum(RAA0*Cweights),
                       min_IPW = min(weights),
                       min_Cali = min(Cweights),
                       min_Cali_sequential = min(Cweights_sequential),
                       max_IPW = max(weights),
                       max_Cali = max(Cweights),
                       max_Cali_sequential = max(Cweights_sequential))
    
    PP_IPW <- TrialEmulation::trial_msm(data = switch_data,
                                    outcome_cov = ~ X1 + X2+ X3 + assigned_treatment+
                                      t_1 + t_2 + t_3 + t_4 +
                                      t_1A + t_2A + t_3A + t_4A +
                                      t_1X1 + t_2X1 + t_3X1 + t_4X1 +
                                      t_1X2 + t_2X2 + t_3X2 + t_4X2 +
                                      t_1X3 + t_2X3 + t_3X3 + t_4X3,
                                    model_var = c('assigned_treatment'),
                                    glm_function = 'glm',
                                    include_trial_period = ~1, include_followup_time = ~1,
                                    estimand_type = 'PP', quiet = T, use_sample_weights =  F)
    result$hr_estimates[1] <- PP_IPW$model$coefficients['assigned_treatment']

    switch_data$weight <- simdatafinal1$data %>% dplyr::filter(RA == 1) %>% dplyr::select(Cweights)
    
    
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
                                               estimand_type = 'PP', quiet = T, use_sample_weights =  F)
    
    result$hr_estimates[2] <- PP_calibrated$model$coefficients['assigned_treatment']
    
    switch_data$weight <- simdatafinal2$data %>% dplyr::filter(RA == 1) %>% dplyr::select(Cweights)
    
    PP_calibrated_seq <- TrialEmulation::trial_msm(data = switch_data,
                                               outcome_cov = ~ X1 + X2+ + X3 + assigned_treatment+
                                                 t_1 + t_2 + t_3 + t_4 +
                                                 t_1A + t_2A + t_3A + t_4A +
                                                 t_1X1 + t_2X1 + t_3X1 + t_4X1 +
                                                 t_1X2 + t_2X2 + t_3X2 + t_4X2 +
                                                 t_1X3 + t_2X3 + t_3X3 + t_4X3,
                                               model_var = c('assigned_treatment'),
                                               glm_function = 'glm',
                                               include_trial_period = ~1, include_followup_time = ~1,
                                               estimand_type = 'PP', quiet = T, use_sample_weights =  F)
    
    result$hr_estimates[3] <- PP_calibrated_seq$model$coefficients['assigned_treatment']
    
    switch_data$weight <- 1.0
    
    PP_naive <- TrialEmulation::trial_msm(data = switch_data,
                                                   outcome_cov = ~ X1 + X2+ + X3 + assigned_treatment+
                                                     t_1 + t_2 + t_3 + t_4 +
                                                     t_1A + t_2A + t_3A + t_4A +
                                                     t_1X1 + t_2X1 + t_3X1 + t_4X1 +
                                                     t_1X2 + t_2X2 + t_3X2 + t_4X2 +
                                                     t_1X3 + t_2X3 + t_3X3 + t_4X3,
                                                   model_var = c('assigned_treatment'),
                                                   glm_function = 'glm',
                                                   include_trial_period = ~1, include_followup_time = ~1,
                                                   estimand_type = 'PP', quiet = T, use_sample_weights =  F)
    
    result$hr_estimates[4] <- PP_naive$model$coefficients['assigned_treatment']
    
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
    
    Y_pred_PP_treatment_IPW <- predict.glm(PP_IPW$model,
                                       fitting_data_treatment,
                                       type = "response")
    Y_pred_PP_control_IPW <- predict.glm(PP_IPW$model,
                                     fitting_data_control,
                                     type = "response")
    Y_pred_PP_treatment_cali <- predict.glm(PP_calibrated$model,
                                            fitting_data_treatment,
                                            type = "response")
    Y_pred_PP_control_cali <- predict.glm(PP_calibrated$model,
                                          fitting_data_control,
                                          type = "response")
    
    Y_pred_PP_treatment_naive <- predict.glm(PP_naive$model,
                                       fitting_data_treatment,
                                       type = "response")
    Y_pred_PP_control_naive <- predict.glm(PP_naive$model,
                                     fitting_data_control,
                                     type = "response")
    Y_pred_PP_treatment_cali_seq <- predict.glm(PP_calibrated_seq$model,
                                            fitting_data_treatment,
                                            type = "response")
    Y_pred_PP_control_cali_seq <- predict.glm(PP_calibrated_seq$model,
                                          fitting_data_control,
                                          type = "response")
    predicted_probas_PP <- fitting_data_treatment %>%
      dplyr::mutate(predicted_proba_treatment_IPW = Y_pred_PP_treatment_IPW,
                    predicted_proba_control_IPW = Y_pred_PP_control_IPW,
                    predicted_proba_treatment_cali = Y_pred_PP_treatment_cali,
                    predicted_proba_control_cali = Y_pred_PP_control_cali,
                    predicted_proba_treatment_naive = Y_pred_PP_treatment_naive,
                    predicted_proba_control_naive = Y_pred_PP_control_naive,
                    predicted_proba_treatment_cali_seq = Y_pred_PP_treatment_cali_seq,
                    predicted_proba_control_cali_seq = Y_pred_PP_control_cali_seq) %>%
      dplyr::group_by(id, trial_period) %>%
      dplyr::mutate(cum_hazard_treatment_IPW = cumprod(1-predicted_proba_treatment_IPW),
                    cum_hazard_control_IPW = cumprod(1-predicted_proba_control_IPW),
                    cum_hazard_treatment_cali = cumprod(1-predicted_proba_treatment_cali),
                    cum_hazard_control_cali = cumprod(1-predicted_proba_control_cali),
                    cum_hazard_treatment_naive = cumprod(1-predicted_proba_treatment_naive),
                    cum_hazard_control_naive = cumprod(1-predicted_proba_control_naive),
                    cum_hazard_treatment_cali_seq = cumprod(1-predicted_proba_treatment_cali_seq),
                    cum_hazard_control_cali_seq = cumprod(1-predicted_proba_control_cali_seq)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(followup_time) %>%
      dplyr::summarise(risk_treatment_IPW = mean(cum_hazard_treatment_IPW),
                       risk_treatment_cali = mean(cum_hazard_treatment_cali),
                       risk_treatment_cali_seq = mean(cum_hazard_treatment_cali_seq),
                       risk_treatment_naive = mean(cum_hazard_treatment_naive),
                       
                       risk_control_IPW = mean(cum_hazard_control_IPW),
                       risk_control_cali = mean(cum_hazard_control_cali),
                       risk_control_cali_seq = mean(cum_hazard_control_cali_seq),
                       risk_control_naive = mean(cum_hazard_control_naive))
    result$mr_1_estimates<- predicted_probas_PP[,2:5]
    result$mr_0_estimates<- predicted_probas_PP[,6:9]
    
    ######### Misspecified ############### 
    simdata_censored <- simdata_censored %>% 
      mutate(Z1 = X1^3/9,Z2 = X1*X2, Z3 = log(abs(X3))+4)
    PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y',cense = 'C', 
                                                eligible ='eligible',
                                                switch_d_cov = ~Z1 + Z2 + Z3,
                                                outcome_cov = ~Z1 + Z2 + Z3, model_var = c('assigned_treatment'),
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
      dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1)) %>% 
      dplyr::mutate(CRA = cumsum(RA),
                    RA = ifelse(CRA == t+1,1,0),
                    RC = ifelse(lag(C) == 0,1,0),
                    A1Z2 = A_0*Z2,
                    A0Z2 = (1-A_0)*Z2,
                    A1Z1 = A_0*Z1,
                    A0Z1 = (1-A_0)*Z1,
                    A1Z3 = A_0*Z3,
                    A0Z3 = (1-A_0)*Z3,
                    A1 = A_0,
                    A0 = 1-A_0,
                    sub = ID,
                    tall = t) %>% 
      merge(dplyr::select(switch_data,id, followup_time, weight), 
            by.x = c('ID', 't'), by.y = c('id', 'followup_time'), all.x = T) %>% 
      dplyr::mutate(weights = ifelse(!is.na(weight), weight,0)) %>% 
      dplyr::arrange(ID, t) 
    
    ################### Calibration aggregated#######################
    simdatafinal1 <- calibration(simdatafinal = data_restric, 
                                 var = c('A1','A1Z1', 'A1Z2','A1Z3',
                                         'A0','A0Z1','A0Z2','A0Z3'))
    
    
    result$objectiveIPWmis <- simdatafinal1$objective.IPW
    result$objectiveCalimis <- simdatafinal1$objective.Cali
    
    ################### Calibration by time #######################
    simdatafinal2 <- calibration_by_time(simdatafinal = data_restric, 
                                         var = c('A1','A1Z1', 'A1Z2','A1Z3',
                                                 'A0','A0Z1','A0Z2','A0Z3'))
    
    
    result$objectiveIPWseqmis <- simdatafinal2$objective.IPW
    result$objectiveCaliseqmis <- simdatafinal2$objective.Cali
    
    meandiffs_summary <- simdatafinal2$data %>% 
      dplyr::mutate(RAX2 = RA*X2,
                    RAX3 = RA*X3,
                    RAX1 = RA*X1,
                    RAA1 = RA*A1,
                    RAA0 = RA*A0,
                    Cweights_sequential = Cweights,
                    Cweights = simdatafinal1$data$Cweights) %>% 
      dplyr::select(t,RA,A1,A0,X2,X3,X1, RAX2, RAX3, RAX1,RAA1, RAA0, weights, Cweights,Cweights_sequential)
    
    treatment_numbers2 <- meandiffs_summary %>% 
      dplyr::group_by(t) %>% 
      dplyr::summarise(Treated_unadjusted_mis = sum(RAA1),
                       Control_unadjusted_mis = sum(RAA0),
                       Treated_IPW_mis = sum(RAA1*weights),
                       Control_IPW_mis = sum(RAA0*weights),
                       Treated_Cali_mis = sum(RAA1 * Cweights),
                       Control_Cali_mis = sum(RAA0*Cweights),
                       Treated_Cali_sequential_mis = sum(RAA1 * Cweights_sequential),
                       Control_Cali_sequential_mis = sum(RAA0*Cweights_sequential),
                       X1_treated_mis = sum(RAX1*A1)/sum(RAA1),
                       X1_control_mis = sum(RAX1*A0)/sum(RAA0),
                       X1_treated_IPW_mis = sum(RAX1*weights*A1)/sum(RAA1*weights),
                       X1_control_IPW_mis = sum(RAX1*weights*A0)/sum(RAA0*weights),
                       X1_treated_Cali_sequential_mis = sum(RAX1*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X1_control_Cali_sequential_mis = sum(RAX1*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X1_treated_Cali_mis = sum(RAX1*Cweights*A1)/sum(RAA1*Cweights),
                       X1_control_Cali_mis = sum(RAX1*Cweights*A0)/sum(RAA0*Cweights),
                       X2_treated_mis = sum(RAX2*A1)/sum(RAA1),
                       X2_control_mis = sum(RAX2*A0)/sum(RAA0),
                       X2_treated_IPW_mis = sum(RAX2*weights*A1)/sum(RAA1*weights),
                       X2_control_IPW_mis = sum(RAX2*weights*A0)/sum(RAA0*weights),
                       X2_treated_Cali_sequential_mis = sum(RAX2*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X2_control_Cali_sequential_mis = sum(RAX2*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X2_treated_Cali_mis = sum(RAX2*Cweights*A1)/sum(RAA1*Cweights),
                       X2_control_Cali_mis = sum(RAX2*Cweights*A0)/sum(RAA0*Cweights),
                       X3_treated = sum(RAX3*A1)/sum(RAA1),
                       X3_control = sum(RAX3*A0)/sum(RAA0),
                       X3_treated_IPW = sum(RAX3*weights*A1)/sum(RAA1*weights),
                       X3_control_IPW = sum(RAX3*weights*A0)/sum(RAA0*weights),
                       X3_treated_Cali_sequential = sum(RAX3*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X3_control_Cali_sequential = sum(RAX3*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X3_treated_Cali = sum(RAX3*Cweights*A1)/sum(RAA1*Cweights),
                       X3_control_Cali = sum(RAX3*Cweights*A0)/sum(RAA0*Cweights),
                       min_IPW_mis = min(weights),
                       min_Cali_mis = min(Cweights),
                       min_Cali_sequential_mis = min(Cweights_sequential),
                       max_IPW_mis = max(weights),
                       max_Cali_mis = max(Cweights),
                       max_Cali_sequential_mis = max(Cweights_sequential))
    
    result$balance_summary <- cbind(treatment_numbers, treatment_numbers2)
    PP_IPW <- TrialEmulation::trial_msm(data = switch_data,
                                        outcome_cov = ~ Z1 + Z2+ Z3 + assigned_treatment+
                                          t_1 + t_2 + t_3 + t_4 +
                                          t_1A + t_2A + t_3A + t_4A +
                                          t_1Z1 + t_2Z1 + t_3Z1 + t_4Z1 +
                                          t_1Z2 + t_2Z2 + t_3Z2 + t_4Z2 +
                                          t_1Z3 + t_2Z3 + t_3Z3 + t_4Z3,
                                        model_var = c('assigned_treatment'),
                                        glm_function = 'glm',
                                        include_trial_period = ~1, include_followup_time = ~1,
                                        estimand_type = 'PP', quiet = T, use_sample_weights =  F)
    result$hr_estimates[5] <- PP_IPW$model$coefficients['assigned_treatment']
    
    switch_data$weight <- simdatafinal1$data %>% dplyr::filter(RA == 1) %>% dplyr::select(Cweights)
    
    
    PP_calibrated <- TrialEmulation::trial_msm(data = switch_data,
                                               outcome_cov = ~ Z1 + Z2+ + Z3 + assigned_treatment+
                                                 t_1 + t_2 + t_3 + t_4 +
                                                 t_1A + t_2A + t_3A + t_4A +
                                                 t_1Z1 + t_2Z1 + t_3Z1 + t_4Z1 +
                                                 t_1Z2 + t_2Z2 + t_3Z2 + t_4Z2 +
                                                 t_1Z3 + t_2Z3 + t_3Z3 + t_4Z3,
                                               model_var = c('assigned_treatment'),
                                               glm_function = 'glm',
                                               include_trial_period = ~1, include_followup_time = ~1,
                                               estimand_type = 'PP', quiet = T, use_sample_weights =  F)
    
    result$hr_estimates[6] <- PP_calibrated$model$coefficients['assigned_treatment']
    
    switch_data$weight <- simdatafinal2$data %>% dplyr::filter(RA == 1) %>% dplyr::select(Cweights)
    
    PP_calibrated_seq <- TrialEmulation::trial_msm(data = switch_data,
                                                   outcome_cov = ~ Z1 + Z2+ + Z3 + assigned_treatment+
                                                     t_1 + t_2 + t_3 + t_4 +
                                                     t_1A + t_2A + t_3A + t_4A +
                                                     t_1Z1 + t_2Z1 + t_3Z1 + t_4Z1 +
                                                     t_1Z2 + t_2Z2 + t_3Z2 + t_4Z2 +
                                                     t_1Z3 + t_2Z3 + t_3Z3 + t_4Z3,
                                                   model_var = c('assigned_treatment'),
                                                   glm_function = 'glm',
                                                   include_trial_period = ~1, include_followup_time = ~1,
                                                   estimand_type = 'PP', quiet = T, use_sample_weights =  F)
    
    result$hr_estimates[7] <- PP_calibrated_seq$model$coefficients['assigned_treatment']
    
    switch_data$weight <- 1.0
    
    PP_naive <- TrialEmulation::trial_msm(data = switch_data,
                                          outcome_cov = ~ Z1 + Z2+ + Z3 + assigned_treatment+
                                            t_1 + t_2 + t_3 + t_4 +
                                            t_1A + t_2A + t_3A + t_4A +
                                            t_1Z1 + t_2Z1 + t_3Z1 + t_4Z1 +
                                            t_1Z2 + t_2Z2 + t_3Z2 + t_4Z2 +
                                            t_1Z3 + t_2Z3 + t_3Z3 + t_4Z3,
                                          model_var = c('assigned_treatment'),
                                          glm_function = 'glm',
                                          include_trial_period = ~1, include_followup_time = ~1,
                                          estimand_type = 'PP', quiet = T, use_sample_weights =  F)
    
    result$hr_estimates[8] <- PP_naive$model$coefficients['assigned_treatment']
    
    design_mat <- expand.grid(id = 1:as.numeric(dim(switch_data)[1]),
                              trial_period = 0:4,
                              followup_time = 0:4)
    design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]
    
    fitting_data_treatment <-  switch_data %>%
      dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>%
      dplyr::select(id,trial_period, followup_time, Z1,  Z2, Z3, assigned_treatment) %>%
      merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>%
      dplyr::group_by(id) %>%
      tidyr::fill( Z1,Z2,Z3, assigned_treatment,.direction = "down") %>%
      dplyr::ungroup() %>%
      dplyr::select(id, trial_period, followup_time, Z1, Z2, Z3, assigned_treatment) %>%
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
    
    Y_pred_PP_treatment_IPW <- predict.glm(PP_IPW$model,
                                           fitting_data_treatment,
                                           type = "response")
    Y_pred_PP_control_IPW <- predict.glm(PP_IPW$model,
                                         fitting_data_control,
                                         type = "response")
    Y_pred_PP_treatment_cali <- predict.glm(PP_calibrated$model,
                                            fitting_data_treatment,
                                            type = "response")
    Y_pred_PP_control_cali <- predict.glm(PP_calibrated$model,
                                          fitting_data_control,
                                          type = "response")
    
    Y_pred_PP_treatment_naive <- predict.glm(PP_naive$model,
                                             fitting_data_treatment,
                                             type = "response")
    Y_pred_PP_control_naive <- predict.glm(PP_naive$model,
                                           fitting_data_control,
                                           type = "response")
    Y_pred_PP_treatment_cali_seq <- predict.glm(PP_calibrated_seq$model,
                                                fitting_data_treatment,
                                                type = "response")
    Y_pred_PP_control_cali_seq <- predict.glm(PP_calibrated_seq$model,
                                              fitting_data_control,
                                              type = "response")
    predicted_probas_PP <- fitting_data_treatment %>%
      dplyr::mutate(predicted_proba_treatment_IPW = Y_pred_PP_treatment_IPW,
                    predicted_proba_control_IPW = Y_pred_PP_control_IPW,
                    predicted_proba_treatment_cali = Y_pred_PP_treatment_cali,
                    predicted_proba_control_cali = Y_pred_PP_control_cali,
                    predicted_proba_treatment_naive = Y_pred_PP_treatment_naive,
                    predicted_proba_control_naive = Y_pred_PP_control_naive,
                    predicted_proba_treatment_cali_seq = Y_pred_PP_treatment_cali_seq,
                    predicted_proba_control_cali_seq = Y_pred_PP_control_cali_seq) %>%
      dplyr::group_by(id, trial_period) %>%
      dplyr::mutate(cum_hazard_treatment_IPW = cumprod(1-predicted_proba_treatment_IPW),
                    cum_hazard_control_IPW = cumprod(1-predicted_proba_control_IPW),
                    cum_hazard_treatment_cali = cumprod(1-predicted_proba_treatment_cali),
                    cum_hazard_control_cali = cumprod(1-predicted_proba_control_cali),
                    cum_hazard_treatment_naive = cumprod(1-predicted_proba_treatment_naive),
                    cum_hazard_control_naive = cumprod(1-predicted_proba_control_naive),
                    cum_hazard_treatment_cali_seq = cumprod(1-predicted_proba_treatment_cali_seq),
                    cum_hazard_control_cali_seq = cumprod(1-predicted_proba_control_cali_seq)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(followup_time) %>%
      dplyr::summarise(risk_treatment_IPW_mis = mean(cum_hazard_treatment_IPW),
                       risk_treatment_cali_mis = mean(cum_hazard_treatment_cali),
                       risk_treatment_cali_seq_mis = mean(cum_hazard_treatment_cali_seq),
                       risk_treatment_naive_mis = mean(cum_hazard_treatment_naive),
                       
                       risk_control_IPW_mis = mean(cum_hazard_control_IPW),
                       risk_control_cali_mis = mean(cum_hazard_control_cali),
                       risk_control_cali_seq_mis = mean(cum_hazard_control_cali_seq),
                       risk_control_naive_mis = mean(cum_hazard_control_naive))
    result$mr_1_estimates<- cbind( result$mr_1_estimates, predicted_probas_PP[,2:5])
    result$mr_0_estimates<-  cbind(result$mr_0_estimates,predicted_probas_PP[,6:9])
    return(result)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
save(oper, file = paste("Simulation results/simulation_ipw_cali_cali_seq_singletrial",as.character(l),".rda", sep = ""))
