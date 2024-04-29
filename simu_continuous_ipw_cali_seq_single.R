#!/usr/bin R
library(dplyr)
library(tidyr)
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("continuous_outcome_datagen.R")
library(TrialEmulation)
library(MASS)
library(sandwich)
library(foreach)
library(doParallel)
library(doRNG)
library(nleqslv)
source('calibration_func_trials.R')
set.seed(NULL)

iters <- 500
#l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
size <- c(200,500,1000)
treat <- c(-1,0,1)
conf <- c(1,1.5,2)

scenarios <- tidyr::crossing(size, treat,conf)
# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 10)
multiResultClass <- function(objectiveIPW = NULL,
                             objectiveCali = NULL,
                             objectiveIPWseq = NULL,
                             objectiveCaliseq = NULL,
                             objectiveIPWCaliT = NULL,
                             objectiveCaliT = NULL,
                             objectiveIPWmis = NULL,
                             objectiveCalimis = NULL,
                             objectiveIPWseqmis = NULL,
                             objectiveCaliseqmis = NULL,
                             objectiveIPWCaliTmiss = NULL,
                             objectiveCaliTmiss = NULL,
                             hr_estimates=NULL,predict_estimates = NULL, balance_summary = NULL)
{
  me <- list(
    objectiveIPW = objectiveIPW,
    objectiveCali = objectiveCali,
    objectiveIPWseq = objectiveIPWseq,
    objectiveCaliseq = objectiveCaliseq,
    objectiveIPWCaliT = objectiveIPWCaliT,
    objectiveCaliT = objectiveCaliT,
    objectiveIPWmis = objectiveIPWmis,
    objectiveCalimis = objectiveCalimis,
    objectiveIPWseqmis = objectiveIPWseqmis,
    objectiveCaliseqmis = objectiveCaliseqmis,
    objectiveIPWCaliTmiss = objectiveIPWCaliTmiss,
    objectiveCaliTmiss = objectiveCaliTmiss,
    hr_estimates = hr_estimates,
    predict_estimates = predict_estimates,
    balance_summary = balance_summary
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}
for (l in 1:27){
oper <- foreach(i = 1:iters,.combine=cbind) %dopar% {
  tryCatch({
    result <- multiResultClass()

    simdata<-DATA_GEN_continous_outcome_treatment_switch(ns = as.numeric(scenarios[l,1]),nv = 5,treat_prev =  as.numeric(scenarios[l,2]),
                                                         conf =  as.numeric(scenarios[l,3]),
                                                         censor = F)
    simdata <- simdata %>% 
      mutate(switch = ifelse(t == 0, 0,ifelse(A!=Ap,1,0)))
    #con4<-xtabs(~t + switch, data=simdata)
    #ftable(con4)
    ######### Correctly specified ############### 
    PP_prep <- TrialEmulation::data_preparation(simdata, id='ID', period='t', treatment='A', outcome='Y', cense = 'C',
                                                eligible ='eligible',
                                                estimand_type = 'PP',
                                                switch_d_cov = ~ X1 + X2 + X3,
                                                switch_n_cov = ~ -1,
                                                outcome_cov = ~ X1 + X2 + X3, model_var = c('assigned_treatment'),
                                                quiet = T,
                                                save_weight_models = F,
                                                data_dir = getwd())
   
    switch_data <- PP_prep$data %>% 
      dplyr::mutate(t = trial_period + followup_time) %>%
      merge(simdata[,c('ID', 't', 'Y')], by.x = c('id', 't'), by.y = c('ID', 't')) %>% 
      dplyr::filter(trial_period == 0) %>% 
      dplyr::arrange(id, followup_time) %>% 
      dplyr::group_by(id) %>% 
      dplyr::mutate(Ap = ifelse(followup_time == 0, 0,lag(assigned_treatment)))
    
    #con4<-xtabs(~followup_time + assigned_treatment, data=switch_data)
    #ftable(con4)
    
    data_restric <- simdata %>% 
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
                    tA1X2 =t* A_0*X2,
                    tA0X2 = t*(1-A_0)*X2,
                    tA1X1 = t*A_0*X1,
                    tA0X1 = t*(1-A_0)*X1,
                    tA1X3 = t*A_0*X3,
                    tA0X3 = t*(1-A_0)*X3,
                    tA1 = t*A_0,
                    tA0 = t*(1-A_0),
                    sub = ID,
                    tall = t,
                    One = 1.0) %>% 
      merge(dplyr::select(switch_data,id, followup_time, weight), 
            by.x = c('ID', 't'), by.y = c('id', 'followup_time'), all.x = T) %>% 
      dplyr::mutate(weights = ifelse(!is.na(weight), weight, 0)) %>% 
      dplyr::arrange(ID, t) 
    
    weight_training_data <- data_restric %>% 
      group_by(ID) %>% 
      filter(cumsum(1-RA) <= 1, t!= 0)
    
    weight_model <- glm(data = weight_training_data, 
                        formula = A ~ Ap + X1 + X3, family = 'binomial')
    #summary(weight_model)
    data_restric$p_1 <- 1.0
    data_restric[data_restric$t != 0,]$p_1 <- predict.glm(weight_model, data_restric[data_restric$t != 0,], type = 'response')
    data_restric <- data_restric %>% 
      arrange(ID, t) %>% 
      group_by(ID) %>% 
      dplyr::mutate(
        wt = ifelse( t == 0, 1.0,ifelse(A == 1, 1/p_1, 1/(1-p_1))),
                     wtprod = cumprod(wt),
        weights = ifelse(weights !=0.0,wtprod,0.0))
   
    ################### Calibration aggregated#######################
    simdatafinal1 <- calibration(simdatafinal = data_restric, 
                                 var = c('A1','A1X1','A1X3',
                                         'A0','A0X1','A0X3'))
    
    
    result$objectiveIPW <- simdatafinal1$objective.IPW
    result$objectiveCali <- simdatafinal1$objective.Cali
    
    ################### Calibration by time #######################
    simdatafinal2 <- calibration_by_time(simdatafinal = data_restric, 
                                         var = c('A1','A1X1', 'A1X3',
                                         'A0','A0X1','A0X3'))
    
    
    result$objectiveIPWseq <- simdatafinal2$objective.IPW
    result$objectiveCaliseq <- simdatafinal2$objective.Cali
    
    ################## Calibration aggregated with time interaction ###########################
    simdatafinal3 <- calibration(simdatafinal = data_restric, 
                                         var = c('tA1','tA1X1','tA1X3',
                                                 'tA0','tA0X1','tA0X3'))
    
    
    result$objectiveIPWCaliT <- simdatafinal3$objective.IPW
    result$objectiveCaliT <- simdatafinal3$objective.Cali
    
    
    meandiffs_summary <- simdatafinal2$data %>% 
      dplyr::mutate(RAX2 = RA*X2,
                    RAX3 = RA*X3,
                    RAX1 = RA*X1,
                    RAA1 = RA*A1,
                    RAA0 = RA*A0,
                    Cweights_sequential = Cweights) %>% 
      dplyr::select(t,RA,A1,A0,X2,X3,X1, RAX2, RAX3, RAX1,RAA1, RAA0, weights, Cweights,Cweights_sequential)
    meandiffs_summary$Cweights <- simdatafinal1$data$Cweights
    meandiffs_summary$CweightsT <- simdatafinal3$data$Cweights

    treatment_numbers <- meandiffs_summary %>% 
      dplyr::group_by(t) %>% 
      dplyr::summarise(Treated_unadjusted = sum(RAA1),
                       Control_unadjusted = sum(RAA0),
                       Treated_notcensored = sum(A1),
                       Control_notcensored = sum(A0),
                       Treated_IPW = sum(RAA1*weights),
                       Control_IPW = sum(RAA0*weights),
                       Treated_Cali = sum(RAA1 * Cweights),
                       Control_Cali = sum(RAA0*Cweights),
                       Treated_Cali_sequential = sum(RAA1 * Cweights_sequential),
                       Control_Cali_sequential = sum(RAA0*Cweights_sequential),
                       Treated_CaliT = sum(RAA1 * CweightsT),
                       Control_CaliT = sum(RAA0*CweightsT),
                       X1_treated_notcensored = sum(X1*A1)/sum(A1),
                       X1_control_notcensored = sum(X1*A0)/sum(A0),
                       X1_treated = sum(RAX1*A1)/sum(RAA1),
                       X1_control = sum(RAX1*A0)/sum(RAA0),
                       X1_treated_IPW = sum(RAX1*weights*A1)/sum(RAA1*weights),
                       X1_control_IPW = sum(RAX1*weights*A0)/sum(RAA0*weights),
                       X1_treated_Cali = sum(RAX1*Cweights*A1)/sum(RAA1*Cweights),
                       X1_control_Cali = sum(RAX1*Cweights*A0)/sum(RAA0*Cweights),
                       X1_treated_Cali_sequential = sum(RAX1*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X1_control_Cali_sequential = sum(RAX1*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X1_treated_CaliT = sum(RAX1*CweightsT*A1)/sum(RAA1*CweightsT),
                       X1_control_CaliT = sum(RAX1*CweightsT*A0)/sum(RAA0*CweightsT),
                       X2_treated_notcensored = sum(X2*A1)/sum(A1),
                       X2_control_notcensored = sum(X2*A0)/sum(A0),
                       X2_treated = sum(RAX2*A1)/sum(RAA1),
                       X2_control = sum(RAX2*A0)/sum(RAA0),
                       X2_treated_IPW = sum(RAX2*weights*A1)/sum(RAA1*weights),
                       X2_control_IPW = sum(RAX2*weights*A0)/sum(RAA0*weights),
                       X2_treated_Cali = sum(RAX2*Cweights*A1)/sum(RAA1*Cweights),
                       X2_control_Cali = sum(RAX2*Cweights*A0)/sum(RAA0*Cweights),
                       X2_treated_Cali_sequential = sum(RAX2*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X2_control_Cali_sequential = sum(RAX2*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X2_treated_CaliT = sum(RAX2*CweightsT*A1)/sum(RAA1*CweightsT),
                       X2_control_CaliT = sum(RAX2*CweightsT*A0)/sum(RAA0*CweightsT),
                       X3_treated_notcensored = sum(X3*A1)/sum(A1),
                       X3_control_notcensored = sum(X3*A0)/sum(A0),
                       X3_treated = sum(RAX3*A1)/sum(RAA1),
                       X3_control = sum(RAX3*A0)/sum(RAA0),
                       X3_treated_IPW = sum(RAX3*weights*A1)/sum(RAA1*weights),
                       X3_control_IPW = sum(RAX3*weights*A0)/sum(RAA0*weights),
                       X3_treated_Cali = sum(RAX3*Cweights*A1)/sum(RAA1*Cweights),
                       X3_control_Cali = sum(RAX3*Cweights*A0)/sum(RAA0*Cweights),
                       X3_treated_Cali_sequential = sum(RAX3*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X3_control_Cali_sequential = sum(RAX3*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X3_treated_CaliT = sum(RAX3*CweightsT*A1)/sum(RAA1*CweightsT),
                       X3_control_CaliT = sum(RAX3*CweightsT*A0)/sum(RAA0*CweightsT),
                       max_IPW = max(weights),
                       max_Cali = max(Cweights),
                       max_Cali_sequential = max(Cweights_sequential),
                       max_CaliT = max(CweightsT))
    
    #con4<-xtabs(~followup_time + assigned_treatment, data=switch_data)
    #ftable(con4)
    switch_data$weight <- simdatafinal1$data[simdatafinal1$data$RA == 1,]$weights
    switch_data$Cweights <- simdatafinal1$data[simdatafinal1$data$RA == 1,]$Cweights
    switch_data$Cweights_sequential <- simdatafinal2$data[simdatafinal2$data$RA == 1,]$Cweights
    switch_data$CweightsT <- simdatafinal3$data[simdatafinal3$data$RA == 1,]$Cweights


    switch_data <- switch_data %>% 
      dplyr::mutate(kAAp = as.numeric(followup_time != 0)*(assigned_treatment + Ap)/2,
             a0 = as.numeric(followup_time == 0)*assigned_treatment,
             X10 = as.numeric(followup_time == 0)*X1,
             X30 = as.numeric(followup_time == 0)*X3,
             X20 = as.numeric(followup_time == 0)*X2)
    
    PP_naive <- glm(data = switch_data,
                    formula = Y ~ a0 + kAAp + X10 + X20 + X30,
                    weights = NULL, family = 'gaussian')
    #summary(PP_naive)
    
    PP_IPW <- glm(data = switch_data,
                  formula = Y ~ a0 + kAAp + X10 + X20 + X30,
                  weights = weight, family = 'gaussian')
    summary(PP_IPW)
    PP_calibrated <- glm(data = switch_data,
                         formula = Y ~ a0 + kAAp + X10 + X20 + X30,
                         weights = Cweights, family = 'gaussian')
    summary(PP_calibrated)
    PP_calibrated_sequential <- glm(data = switch_data,
                                    formula = Y ~ a0 + kAAp + X10 + X20 + X30,
                                    weights = Cweights_sequential, family = 'gaussian')
    summary(PP_calibrated_sequential)
    PP_calibrated_T <- glm(data = switch_data,
                                    formula = Y ~ a0 + kAAp + X10 + X20 + X30,
                                    weights = CweightsT, family = 'gaussian')
    summary(PP_calibrated_T)
  
    
    result$hr_estimates<- array(,dim = c(10,6))
    result$hr_estimates[1,] <- PP_naive$coefficients
    result$hr_estimates[2,] <- PP_IPW$coefficients
    result$hr_estimates[3,] <- PP_calibrated$coefficients
    result$hr_estimates[4,] <- PP_calibrated_sequential$coefficients
    result$hr_estimates[5,] <- PP_calibrated_T$coefficients

    design_mat <- expand.grid(id = 1:as.numeric(dim(switch_data)[1]),
                              trial_period = 0:4,
                              followup_time = 0:4)
    design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]
    
    fitting_data_treatment <-  switch_data %>%
      dplyr::mutate(assigned_treatment = 1.0) %>%
      dplyr::select(id,trial_period, followup_time, X1,  X2, X3, assigned_treatment) %>%
      merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>%
      dplyr::group_by(id) %>%
      tidyr::fill( X1,X2,X3, assigned_treatment,.direction = "down") %>%
      dplyr::ungroup() %>%
      dplyr::select(id, trial_period, followup_time, X1, X2, X3, assigned_treatment) %>%
      merge(data.frame(id = switch_data$id, trial_period = switch_data$trial_period), by = c("id", "trial_period"), all.y = TRUE) %>%
      dplyr::arrange(id, trial_period, followup_time) %>%
      dplyr::distinct(id, trial_period, followup_time, .keep_all = T) %>% 
      dplyr::group_by(id,trial_period) %>%
      dplyr::mutate(kAAp = as.numeric(followup_time != 0)*(assigned_treatment + 1.0)/2,
                    a0 = as.numeric(followup_time == 0)*assigned_treatment,
                    X10 = as.numeric(followup_time == 0)*X1,
                    X30 = as.numeric(followup_time == 0)*X3,
                    X20 = as.numeric(followup_time == 0)*X2) %>%
      dplyr::filter(trial_period == 0) %>% 
      dplyr::ungroup()
    
    fitting_data_control <- fitting_data_treatment %>%
      dplyr::mutate(kAAp = as.numeric(followup_time != 0)*(0.0 + 0.0)/2,
                    a0 = as.numeric(followup_time == 0)*0.0,
                    X10 = as.numeric(followup_time == 0)*X1,
                    X30 = as.numeric(followup_time == 0)*X3,
                    X20 = as.numeric(followup_time == 0)*X2)
    
    Y_pred_PP_treatment_naive <- predict.glm(PP_naive,
                                           fitting_data_treatment)
    Y_pred_PP_control_naive <- predict.glm(PP_naive,
                                         fitting_data_control)
    Y_pred_PP_treatment_IPW <- predict.glm(PP_IPW,
                                           fitting_data_treatment)
    Y_pred_PP_control_IPW <- predict.glm(PP_IPW,
                                         fitting_data_control)
    Y_pred_PP_treatment_cali <- predict.glm(PP_calibrated,
                                            fitting_data_treatment)
    Y_pred_PP_control_cali <- predict.glm(PP_calibrated,
                                          fitting_data_control)
    Y_pred_PP_treatment_cali_sequential <- predict.glm(PP_calibrated_sequential,
                                                       fitting_data_treatment)
    Y_pred_PP_control_cali_sequential <- predict.glm(PP_calibrated_sequential,
                                                     fitting_data_control)
    Y_pred_PP_treatment_cali_T <- predict.glm(PP_calibrated_T,
                                            fitting_data_treatment)
    Y_pred_PP_control_cali_T <- predict.glm(PP_calibrated_T,
                                          fitting_data_control)
    
    predicted_probas_PP <- fitting_data_treatment %>%
      dplyr::mutate(predicted_proba_treatment_naive = Y_pred_PP_treatment_naive,
                    predicted_proba_control_naive = Y_pred_PP_control_naive,
                    predicted_proba_treatment_IPW = Y_pred_PP_treatment_IPW,
                    predicted_proba_control_IPW = Y_pred_PP_control_IPW,
                    predicted_proba_treatment_cali = Y_pred_PP_treatment_cali,
                    predicted_proba_control_cali = Y_pred_PP_control_cali,
                    predicted_proba_treatment_cali_sequential = Y_pred_PP_treatment_cali_sequential,
                    predicted_proba_control_cali_sequential = Y_pred_PP_control_cali_sequential,
                    predicted_proba_treatment_cali_T = Y_pred_PP_treatment_cali_T,
                    predicted_proba_control_cali_T= Y_pred_PP_control_cali_T) %>%
      dplyr::group_by(trial_period, followup_time) %>%
      dplyr::summarise(expected_treatment_naive = mean(predicted_proba_treatment_naive),
                       expected_control_naive = mean(predicted_proba_control_naive),
                       ate_naive = expected_treatment_naive - expected_control_naive,
                       expected_treatment_IPW = mean(predicted_proba_treatment_IPW),
                       expected_control_IPW = mean(predicted_proba_control_IPW),
                       ate_IPW = expected_treatment_IPW - expected_control_IPW,
                       expected_treatment_cali = mean(predicted_proba_treatment_cali),
                       expected_control_cali = mean(predicted_proba_control_cali),
                       ate_cali = expected_treatment_cali - expected_control_cali,
                       expected_treatment_cali_sequential = mean(predicted_proba_treatment_cali_sequential),
                       expected_control_cali_sequential = mean(predicted_proba_control_cali_sequential),
                       ate_cali_sequential = expected_treatment_cali_sequential - expected_control_cali_sequential,
                       expected_treatment_cali_T = mean(predicted_proba_treatment_cali_T),
                       expected_control_cali_T = mean(predicted_proba_control_cali_T),
                       ate_cali_T = expected_treatment_cali_T - expected_control_cali_T)
    
    
    result$predict_estimates <- predicted_probas_PP[,3:17]
    
    ######### Misspecified ############### 
    simdata <- simdata %>% 
      mutate(Z1 =X1/(1+exp(X3))+10,Z2 =(X2/25+0.6)^3, Z3 = exp(X3/2))
    
    PP_prep <- TrialEmulation::data_preparation(simdata, id='ID', period='t', treatment='A', outcome='Y', cense = 'C',
                                                eligible ='eligible',
                                                estimand_type = 'PP',
                                                switch_d_cov = ~Z1 +Z2 + Z3,
                                                switch_n_cov = ~ -1,
                                                outcome_cov = ~Z1 + Z2 + Z3, model_var = c('assigned_treatment'),
                                                quiet = T,
                                                save_weight_models = F,
                                                data_dir = getwd())
    
    switch_data <- PP_prep$data %>% 
      mutate(t = trial_period + followup_time) %>%
      merge(simdata[,c('ID', 't', 'Y')], by.x = c('id', 't'), by.y = c('ID', 't')) %>% 
      dplyr::filter(trial_period == 0) %>% 
      dplyr::arrange(id, followup_time) %>% 
      dplyr::group_by(id) %>% 
      dplyr::mutate(CA = cumsum(assigned_treatment),
                    Ap = ifelse(followup_time == 0, 0,lag(assigned_treatment)))
    
    
    data_restric <- simdata %>% 
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
                    A1X2 = A_0*X2,
                    A0X2 = (1-A_0)*X2,
                    A1X1 = A_0*X1,
                    A0X1 = (1-A_0)*X1,
                    A1X3 = A_0*X3,
                    A0X3 = (1-A_0)*X3,
                    A1 = A_0,
                    A0 = 1-A_0,
                    tA1Z2 =t* A_0*Z2,
                    tA0Z2 = t*(1-A_0)*Z2,
                    tA1Z1 = t*A_0*Z1,
                    tA0Z1 = t*(1-A_0)*Z1,
                    tA1Z3 = t*A_0*Z3,
                    tA0Z3 = t*(1-A_0)*Z3,
                    tA1 = t*A_0,
                    tA0 = t*(1-A_0),
                    sub = ID,
                    tall = t,
                    One = 1.0) %>% 
      merge(dplyr::select(switch_data,id, followup_time, weight), 
            by.x = c('ID', 't'), by.y = c('id', 'followup_time'), all.x = T) %>% 
      dplyr::mutate(weights = ifelse(!is.na(weight), weight, 0)) %>% 
      dplyr::arrange(ID, t) 
    
    weight_training_data <- data_restric %>% 
      group_by(ID) %>% 
      filter(cumsum(1-RA) <= 1, t!= 0)
    
    weight_model <- glm(data = weight_training_data, 
                        formula = A ~ Ap + Z1 + Z3, family = 'binomial')
    #summary(weight_model)
    data_restric$p_1 <- 1.0
    data_restric[data_restric$t != 0,]$p_1 <- predict.glm(weight_model, data_restric[data_restric$t != 0,], type = 'response')
    data_restric <- data_restric %>% 
      arrange(ID, t) %>% 
      group_by(ID) %>% 
      dplyr::mutate(
        wt = ifelse( t == 0, 1.0,ifelse(A == 1, 1/p_1, 1/(1-p_1))),
        wtprod = cumprod(wt),
        weights = ifelse(weights !=0.0,wtprod,0.0))
    
    ################### Calibration aggregated#######################
    simdatafinal1 <- calibration(simdatafinal = data_restric, 
                                 var = c('A1','A1Z1', 'A1Z3',
                                         'A0','A0Z1','A0Z3'))
    
    
    result$objectiveIPWmis <- simdatafinal1$objective.IPW
    result$objectiveCalimis <- simdatafinal1$objective.Cali
    
    ################### Calibration by time #######################
    simdatafinal2 <- calibration_by_time(simdatafinal = data_restric, 
                                         var = c('A1','A1Z1', 'A1Z3',
                                                 'A0','A0Z1','A0Z3'))
    
    
    result$objectiveIPWseqmis <- simdatafinal2$objective.IPW
    result$objectiveCaliseqmis <- simdatafinal2$objective.Cali
    
    ################## Calibration aggregated with time interaction ###########################
    simdatafinal3 <- calibration(simdatafinal = data_restric, 
                                 var = c('tA1','tA1Z1', 'tA1Z3',
                                         'tA0','tA0Z1','tA0Z3'))
    
    
    result$objectiveIPWCaliT <- simdatafinal3$objective.IPW
    result$objectiveCaliT <- simdatafinal3$objective.Cali
    
  
    
    meandiffs_summary <- simdatafinal2$data %>% 
      dplyr::mutate(RAX2 = RA*X2,
                    RAX3 = RA*X3,
                    RAX1 = RA*X1,
                    RAA1 = RA*A1,
                    RAA0 = RA*A0,
                    Cweights_sequential = Cweights) %>% 
      dplyr::select(t,RA,A1,A0,X2,X3,X1, RAX2, RAX3, RAX1,RAA1, RAA0, weights, Cweights,Cweights_sequential)
    meandiffs_summary$Cweights <- simdatafinal1$data$Cweights
    meandiffs_summary$CweightsT <- simdatafinal3$data$Cweights

    treatment_numbers_mis <- meandiffs_summary %>% 
      dplyr::group_by(t) %>% 
      dplyr::summarise(
                       Treated_IPWmis = sum(RAA1*weights),
                       Control_IPWmis = sum(RAA0*weights),
                       Treated_Calimis = sum(RAA1 * Cweights),
                       Control_Calimis = sum(RAA0*Cweights),
                       Treated_Cali_sequentialmis = sum(RAA1 * Cweights_sequential),
                       Control_Cali_sequentialmis = sum(RAA0*Cweights_sequential),
                       Treated_CaliTmis = sum(RAA1 * CweightsT),
                       Control_CaliTmis = sum(RAA0*CweightsT),
                       X1_treated_IPWmis = sum(RAX1*weights*A1)/sum(RAA1*weights),
                       X1_control_IPWmis = sum(RAX1*weights*A0)/sum(RAA0*weights),
                       X1_treated_Calimis = sum(RAX1*Cweights*A1)/sum(RAA1*Cweights),
                       X1_control_Calimis = sum(RAX1*Cweights*A0)/sum(RAA0*Cweights),
                       X1_treated_Cali_sequentialmis = sum(RAX1*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X1_control_Cali_sequentialmis = sum(RAX1*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X1_treated_CaliTmis = sum(RAX1*CweightsT*A1)/sum(RAA1*CweightsT),
                       X1_control_CaliTmis = sum(RAX1*CweightsT*A0)/sum(RAA0*CweightsT),
                       X2_treated_IPWmis = sum(RAX2*weights*A1)/sum(RAA1*weights),
                       X2_control_IPWmis = sum(RAX2*weights*A0)/sum(RAA0*weights),
                       X2_treated_Cali_sequentialmis = sum(RAX2*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X2_control_Cali_sequentialmis = sum(RAX2*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X2_treated_Calimis = sum(RAX2*Cweights*A1)/sum(RAA1*Cweights),
                       X2_control_Calimis = sum(RAX2*Cweights*A0)/sum(RAA0*Cweights),
                       X2_treated_CaliTmis = sum(RAX2*CweightsT*A1)/sum(RAA1*CweightsT),
                       X2_control_CaliTmis = sum(RAX2*CweightsT*A0)/sum(RAA0*CweightsT),
                       X3_treated_IPWmis = sum(RAX3*weights*A1)/sum(RAA1*weights),
                       X3_control_IPWmis = sum(RAX3*weights*A0)/sum(RAA0*weights),
                       X3_treated_Cali_sequentialmis = sum(RAX3*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                       X3_control_Cali_sequentialmis = sum(RAX3*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                       X3_treated_Calimis = sum(RAX3*Cweights*A1)/sum(RAA1*Cweights),
                       X3_control_Calimis = sum(RAX3*Cweights*A0)/sum(RAA0*Cweights),
                       X3_treated_CaliTmis = sum(RAX3*CweightsT*A1)/sum(RAA1*CweightsT),
                       X3_control_CaliTmis = sum(RAX3*CweightsT*A0)/sum(RAA0*CweightsT),
                       max_IPWmis = max(weights),
                       max_Calimis = max(Cweights),
                       max_Cali_sequentialmis = max(Cweights_sequential),
                       max_CaliTmis = max(CweightsT),)
    
    result$balance_summary <- cbind(treatment_numbers,treatment_numbers_mis)
    
    
    switch_data$weight <- simdatafinal1$data[simdatafinal1$data$RA == 1,]$weights
    switch_data$Cweights <- simdatafinal1$data[simdatafinal1$data$RA == 1,]$Cweights
    switch_data$Cweights_sequential <- simdatafinal2$data[simdatafinal2$data$RA == 1,]$Cweights
    switch_data$CweightsT <- simdatafinal3$data[simdatafinal3$data$RA == 1,]$Cweights
    
    switch_data <- switch_data %>% 
      dplyr::mutate(kAAp = as.numeric(followup_time != 0)*(assigned_treatment + Ap)/2,
                    a0 = as.numeric(followup_time == 0)*assigned_treatment,
                    Z10 = as.numeric(followup_time == 0)*Z1,
                    Z30 = as.numeric(followup_time == 0)*Z3,
                    Z20 = as.numeric(followup_time == 0)*Z2)
    
    PP_naive <- glm(data = switch_data,
                    formula = Y ~ a0 + kAAp + Z10 + Z20 + Z30,
                    weights = NULL, family = 'gaussian')
    summary(PP_naive)
    
    PP_IPW <- glm(data = switch_data,
                  formula = Y ~ a0 + kAAp + Z10 + Z20 + Z30,
                  weights = weight, family = 'gaussian')
    summary(PP_IPW)
    PP_calibrated <- glm(data = switch_data,
                         formula = Y ~ a0 + kAAp + Z10 + Z20 + Z30,
                         weights = Cweights, family = 'gaussian')
    summary(PP_calibrated)
    PP_calibrated_sequential <- glm(data = switch_data,
                                    formula = Y ~ a0 + kAAp + Z10 + Z20 + Z30,
                                    weights = Cweights_sequential, family = 'gaussian')
    summary(PP_calibrated_sequential)
    PP_calibrated_T <- glm(data = switch_data,
                           formula = Y ~ a0 + kAAp + Z10 + Z20 + Z30,
                           weights = CweightsT, family = 'gaussian')
    summary(PP_calibrated_T)
    
    
    result$hr_estimates[6,] <- PP_naive$coefficients
    result$hr_estimates[7,] <- PP_IPW$coefficients
    result$hr_estimates[8,] <- PP_calibrated$coefficients
    result$hr_estimates[9,] <- PP_calibrated_sequential$coefficients
    result$hr_estimates[10,] <- PP_calibrated_T$coefficients

    design_mat <- expand.grid(id = 1:as.numeric(dim(switch_data)[1]),
                              trial_period = 0:4,
                              followup_time = 0:4)
    design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]
    
    fitting_data_treatment <-  switch_data %>%
      dplyr::mutate(assigned_treatment = 1.0) %>%
      dplyr::select(id,trial_period, followup_time, Z1,  Z2, Z3, assigned_treatment) %>%
      merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>%
      dplyr::group_by(id) %>%
      tidyr::fill( Z1,Z2,Z3, assigned_treatment,.direction = "down") %>%
      dplyr::ungroup() %>%
      dplyr::select(id, trial_period, followup_time, Z1, Z2, Z3, assigned_treatment) %>%
      merge(data.frame(id = switch_data$id, trial_period = switch_data$trial_period), by = c("id", "trial_period"), all.y = TRUE) %>%
      dplyr::arrange(id, trial_period, followup_time) %>%
      dplyr::distinct(id, trial_period, followup_time, .keep_all = T) %>% 
      dplyr::group_by(id,trial_period) %>%
      dplyr::mutate(kAAp = as.numeric(followup_time != 0)*(assigned_treatment + 1.0)/2,
                    a0 = as.numeric(followup_time == 0)*assigned_treatment,
                    Z10 = as.numeric(followup_time == 0)*Z1,
                    Z30 = as.numeric(followup_time == 0)*Z3,
                    Z20 = as.numeric(followup_time == 0)*Z2) %>%
      dplyr::filter(trial_period == 0) %>% 
      dplyr::ungroup()
    
    fitting_data_control <- fitting_data_treatment %>%
      dplyr::mutate(kAAp = as.numeric(followup_time != 0)*(0.0 + 0.0)/2,
                    a0 = as.numeric(followup_time == 0)*0.0,
                    Z10 = as.numeric(followup_time == 0)*Z1,
                    Z30 = as.numeric(followup_time == 0)*Z3,
                    Z20 = as.numeric(followup_time == 0)*Z2)
    
    Y_pred_PP_treatment_naive <- predict.glm(PP_naive,
                                             fitting_data_treatment)
    Y_pred_PP_control_naive <- predict.glm(PP_naive,
                                           fitting_data_control)
    Y_pred_PP_treatment_IPW <- predict.glm(PP_IPW,
                                           fitting_data_treatment)
    Y_pred_PP_control_IPW <- predict.glm(PP_IPW,
                                         fitting_data_control)
    Y_pred_PP_treatment_cali <- predict.glm(PP_calibrated,
                                            fitting_data_treatment)
    Y_pred_PP_control_cali <- predict.glm(PP_calibrated,
                                          fitting_data_control)
    Y_pred_PP_treatment_cali_sequential <- predict.glm(PP_calibrated_sequential,
                                                       fitting_data_treatment)
    Y_pred_PP_control_cali_sequential <- predict.glm(PP_calibrated_sequential,
                                                     fitting_data_control)
    Y_pred_PP_treatment_cali_T <- predict.glm(PP_calibrated_T,
                                              fitting_data_treatment)
    Y_pred_PP_control_cali_T <- predict.glm(PP_calibrated_T,
                                            fitting_data_control)
    predicted_probas_PP <- fitting_data_treatment %>%
      dplyr::mutate(predicted_proba_treatment_naive = Y_pred_PP_treatment_naive,
                    predicted_proba_control_naive = Y_pred_PP_control_naive,
                    predicted_proba_treatment_IPW = Y_pred_PP_treatment_IPW,
                    predicted_proba_control_IPW = Y_pred_PP_control_IPW,
                    predicted_proba_treatment_cali = Y_pred_PP_treatment_cali,
                    predicted_proba_control_cali = Y_pred_PP_control_cali,
                    predicted_proba_treatment_cali_sequential = Y_pred_PP_treatment_cali_sequential,
                    predicted_proba_control_cali_sequential = Y_pred_PP_control_cali_sequential,
                    predicted_proba_treatment_cali_T = Y_pred_PP_treatment_cali_T,
                    predicted_proba_control_cali_T= Y_pred_PP_control_cali_T) %>%
      dplyr::group_by(trial_period, followup_time) %>%
      dplyr::summarise(expected_treatment_naive_miss = mean(predicted_proba_treatment_naive),
                       expected_control_naive_miss = mean(predicted_proba_control_naive),
                       ate_naive_miss = expected_treatment_naive_miss - expected_control_naive_miss,
                       expected_treatment_IPW_miss = mean(predicted_proba_treatment_IPW),
                       expected_control_IPW_miss = mean(predicted_proba_control_IPW),
                       ate_IPW_miss = expected_treatment_IPW_miss - expected_control_IPW_miss,
                       expected_treatment_cali_miss = mean(predicted_proba_treatment_cali),
                       expected_control_cali_miss = mean(predicted_proba_control_cali),
                       ate_cali_miss = expected_treatment_cali_miss - expected_control_cali_miss,
                       expected_treatment_cali_sequential_miss = mean(predicted_proba_treatment_cali_sequential),
                       expected_control_cali_sequential_miss = mean(predicted_proba_control_cali_sequential),
                       ate_cali_sequential_miss = expected_treatment_cali_sequential_miss - expected_control_cali_sequential_miss,
                       expected_treatment_cali_T_miss = mean(predicted_proba_treatment_cali_T),
                       expected_control_cali_T_miss = mean(predicted_proba_control_cali_T),
                       ate_cali_T_miss = expected_treatment_cali_T_miss - expected_control_cali_T_miss)
    
    
    
    
    result$predict_estimates <- cbind(result$predict_estimates,predicted_probas_PP[,3:17])
    return(result)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
save(oper, file = paste("Simulation results/result_simu_continuous_ipw_cali_seq_single_",as.character(l),".rda", sep = ""))
}
