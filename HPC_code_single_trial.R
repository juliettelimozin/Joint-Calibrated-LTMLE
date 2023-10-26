#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
#setwd("~/rds/hpc-work")
source("simulate_MSM_simplified.R")
library(TrialEmulation)
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)
library(cobalt)
library(nleqslv)
source('calibration_func_trials.R')
set.seed(NULL)

iters <- 500
#l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(1,3,5)

scenarios <- tidyr::crossing(size,conf, treat)
for (l in 20:27){
  absmeandiffs <- array(,dim = c(3,2,4,iters))
  objectives <- array(,dim = c(2,6,iters))
  hr_estimates <- array(,dim = c(2,iters))
  mrd_estimates <- array(,dim = c(2,5,iters))
  
  #data_direction <- paste("~/rds/hpc-work/models_scenario_med_",l,sep = "")
  # Set number of cores. 67 is sufficient for 200 cores.
  #registerDoParallel(cores = 67)
  
  for (i in 1:iters){
    tryCatch({
      simdata_censored<-DATA_GEN_censored_reduced(as.numeric(scenarios[l,1]), 5, 
                                                  conf = as.numeric(scenarios[l,2]), 
                                                  treat_prev = as.numeric(scenarios[l,3]),
                                                  outcome_prev = -4.7,
                                                  censor = F)
          
      PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                                  eligible ='eligible',
                                                  switch_d_cov = ~X2 + X4,
                                                  outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
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
                      t_1X2 = t_1*X2,
                      t_2X2 = t_2*X2,
                      t_3X2 = t_3*X2,
                      t_4X2 = t_4*X2,
                      t_1X4 = t_1*X4,
                      t_2X4 = t_2*X4,
                      t_3X4 = t_3*X4,
                      t_4X4 = t_4*X4)
      
      
      data_restric <- simdata_censored %>% 
        dplyr::group_by(ID) %>% 
        dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1), nextX2 = lead(X2),prevX2 = lag(X2)) %>% 
        dplyr::mutate(CRA = cumsum(RA),
                      nextX2 = ifelse(is.na(nextX2), 0, nextX2),
                      RA = ifelse(CRA == t+1,1,0),
                      RC = ifelse(lag(C) == 0,1,0),
                      t0 = ifelse(t == 0,1,0),
                      t1 = ifelse(t == 1,1,0),
                      t2 = ifelse(t == 2, 1, 0),
                      A1X4 = A_0*X4,
                      A0X4 = (1-A_0)*X4,
                      A1X2 = A_0*X2,
                      A0X2 = (1-A_0)*X2,
                      A1nextX2 = A_0*nextX2,
                      A0nextX2 = (1-A_0)*nextX2,
                      A1 = A_0,
                      A0 = 1-A_0,
                      sub = ID,
                      tall = t) %>% 
        dplyr::filter(RA == 1) %>% 
        merge(dplyr::select(switch_data,id, followup_time, weight), 
              by.x = c('ID', 't'), by.y = c('id', 'followup_time')) %>% 
        dplyr::mutate(weights = weight) %>% 
        dplyr::arrange(ID, t) 
      
      simdatafinal <- calibration(simdatafinal = data_restric, 
                                  var = c('A1', 'A0', 'A1X4', 'A0X4', 
                                          'A1nextX2', 'A0nextX2'))
      
      objectives[1,,i] <- simdatafinal$objective.IPW
      objectives[2,,i] <- simdatafinal$objective.Cali
      
      balance.tab1 <- bal.tab(simdatafinal$data[simdatafinal$data$t == 1,c('X2', 'X4')],treat = simdatafinal$data[simdatafinal$data$t == 1,]$A,
                              stats = c('m'),var.name = c('X2', 'X4'), un = TRUE,
                              weights = list(MLE = simdatafinal$data[simdatafinal$data$t == 1,]$weights, 
                                             Cali= simdatafinal$data[simdatafinal$data$t == 1,]$Cweights),
                              abs = TRUE)
      absmeandiffs[,,1,i] <- do.call(rbind, lapply(balance.tab1$Balance[,-1], as.numeric))
      
      try({balance.tab2 <- bal.tab(simdatafinal$data[simdatafinal$data$t == 2,c('X2', 'X4')],treat = simdatafinal$data[simdatafinal$data$t == 2,]$A,
                                  stats = c('m'),var.name = c('X2', 'X4'), un = TRUE,
                                  weights = list(MLE = simdatafinal$data[simdatafinal$data$t == 2,]$weights, 
                                                 Cali= simdatafinal$data[simdatafinal$data$t == 2,]$Cweights),
                                  abs = TRUE)
      absmeandiffs[,,2,i] <-do.call(rbind, lapply(balance.tab2$Balance[,-1], as.numeric))},
                                silent = T)
      
       try({balance.tab3 <- bal.tab(simdatafinal$data[simdatafinal$data$t == 3,c('X2', 'X4')],treat = simdatafinal$data[simdatafinal$data$t == 3,]$A,
                                  stats = c('m'),var.name = c('X2', 'X4'), un = TRUE,
                                  weights = list(MLE = simdatafinal$data[simdatafinal$data$t == 3,]$weights, 
                                                 Cali= simdatafinal$data[simdatafinal$data$t == 3,]$Cweights),
                                  abs = TRUE)
       absmeandiffs[,,3,i] <-do.call(rbind, lapply(balance.tab3$Balance[,-1], as.numeric))},
                                 silent = T)
      
      try({balance.tab4 <- bal.tab(simdatafinal$data[simdatafinal$data$t == 4,c('X2', 'X4')],treat = simdatafinal$data[simdatafinal$data$t == 4,]$A,
                                  stats = c('m'),var.name = c('X2', 'X4'), un = TRUE,
                                  weights = list(MLE = simdatafinal$data[simdatafinal$data$t == 4,]$weights, 
                                                 Cali= simdatafinal$data[simdatafinal$data$t == 4,]$Cweights),
                                  abs = TRUE)
        
      absmeandiffs[,,4,i] <- do.call(rbind, lapply(balance.tab4$Balance[,-1], as.numeric))},
                                silent = T)
      
      
      PP <- TrialEmulation::trial_msm(data = switch_data,
                                      outcome_cov = ~ X2 + X4+ assigned_treatment+
                                        t_1 + t_2 + t_3 + t_4 +
                                        t_1A + t_2A + t_3A + t_4A +
                                        t_1X2 + t_2X2 + t_3X2 + t_4X2 +
                                        t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                      model_var = c('assigned_treatment'),
                                      glm_function = 'glm',
                                      include_trial_period = ~1, include_followup_time = ~1,
                                      use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
      hr_estimates[1,i] <- PP$model$coefficients[2]
      
      switch_data$weight <- simdatafinal$data$Cweights
      PP_calibrated <- TrialEmulation::trial_msm(data = switch_data,
                                                 outcome_cov = ~ X2 + X4+ assigned_treatment+
                                                   t_1 + t_2 + t_3 + t_4 +
                                                   t_1A + t_2A + t_3A + t_4A +
                                                   t_1X2 + t_2X2 + t_3X2 + t_4X2 +
                                                   t_1X4 + t_2X4 + t_3X4 + t_4X4,
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
        dplyr::select(id,trial_period, followup_time, X2,  X4, assigned_treatment) %>%
        merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>%
        dplyr::group_by(id) %>%
        tidyr::fill( X2,X4,assigned_treatment,.direction = "down") %>%
        dplyr::ungroup() %>%
        dplyr::select(id, trial_period, followup_time, X2, X4, assigned_treatment) %>%
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
                      t_1X2 = t_1*X2,
                      t_2X2 = t_2*X2,
                      t_3X2 = t_3*X2,
                      t_4X2 = t_4*X2,
                      t_1X4 = t_1*X4,
                      t_2X4 = t_2*X4,
                      t_3X4 = t_3*X4,
                      t_4X4 = t_4*X4) %>%
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
        dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                         survival_control = mean(cum_hazard_control),
                         mrd = survival_control - survival_treatment,
                         survival_treatment_cali = mean(cum_hazard_treatment_cali),
                         survival_control_cali = mean(cum_hazard_control_cali),
                         mrd_cali = survival_control_cali - survival_treatment_cali)
      mrd_estimates[1,,i] <- pull(predicted_probas_PP,mrd)
      mrd_estimates[2,,i] <- pull(predicted_probas_PP,mrd_cali)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  save(absmeandiffs, file = paste("absmeandiffs_singletrial_",as.character(l),".rda", sep = ""))
  save(objectives, file = paste("objectives_singletrial_",as.character(l),".rda", sep = ""))
  save(hr_estimates, file = paste("hr_estimates_singletrial_",as.character(l),".rda", sep = ""))
  save(mrd_estimates, file = paste("mrd_estimates_singletrial_",as.character(l),".rda", sep = ""))
}