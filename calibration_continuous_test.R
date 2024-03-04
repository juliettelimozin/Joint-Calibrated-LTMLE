library(modelr)
library(tidyverse)
library(tidyr)
source("continuous_outcome_datagen.R")
library(TrialEmulation)
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)
library(nleqslv)
source('calibration_func_trials.R')
set.seed(NULL)
library(ggplot2)
library(ggpubr)
library(calculus)

simdata<-DATA_GEN_continous_outcome_treatment_switch(5000, 5, 
                                            treat_prev = -1,
                                            outcome_prev = -3,
                                            censor = F)

######### Correctly specified ############### 
PP_prep <- TrialEmulation::data_preparation(simdata, id='ID', period='t', treatment='A', outcome='Y', cense = 'C',
                                            eligible ='eligible',
                                            estimand_type = 'PP',
                                            switch_d_cov = ~X1 + X2 + X3,
                                            outcome_cov = ~X1 + X2 + X3, model_var = c('assigned_treatment'),
                                            quiet = T,
                                            save_weight_models = F,
                                            data_dir = getwd())

switch_data <- PP_prep$data %>% 
  mutate(t = trial_period + followup_time) %>%
  merge(simdata[,c('ID', 't', 'Y')], by.x = c('id', 't'), by.y = c('ID', 't')) %>% 
  dplyr::filter(trial_period == 0) %>% 
  dplyr::arrange(id, followup_time) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(CA = cumsum(assigned_treatment))


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
                sub = ID,
                tall = t) %>% 
  merge(dplyr::select(switch_data,id, followup_time, weight), 
        by.x = c('ID', 't'), by.y = c('id', 'followup_time'), all.x = T) %>% 
  dplyr::mutate(weights = ifelse(!is.na(weight), weight,0)) %>% 
  dplyr::arrange(ID, t) 

simdatafinal <- calibration_by_time(simdatafinal = data_restric, 
                                    var = c('A1','A1X1', 'A1X2','A1X3',
                                            'A0','A0X1','A0X2','A0X3'))
simdatafinal2 <- calibration(simdatafinal = data_restric, 
                             var = c('A1','A1X1', 'A1X2','A1X3',
                                     'A0','A0X1','A0X2','A0X3'))

meandiffs_summary <- simdatafinal$data %>% 
  dplyr::mutate(RAX2 = RA*X2,
                RAX3 = RA*X3,
                RAX1 = RA*X1,
                RAA1 = RA*A1,
                RAA0 = RA*A0,
                Cweights_sequential = Cweights,
                Cweights = simdatafinal2$data$Cweights) %>% 
  dplyr::select(t,RA,A1,A0,X2,X3,X1,RAA1,RAA0, RAX2, RAX3, RAX1, weights, Cweights,Cweights_sequential)


treatment_numbers <- meandiffs_summary %>% 
  dplyr::group_by(t) %>% 
  dplyr::summarise(Treated_notcensored = sum(A1),
                   Control_notcensored = sum(A0),
                   Treated_unadjusted = sum(RAA1),
                   Control_unadjusted = sum(RAA0),
                   Treated_IPW = sum(RAA1*weights),
                   Control_IPW = sum(RAA0*weights),
                   Treated_Cali = sum(RAA1 * Cweights),
                   Control_Cali = sum(RAA0*Cweights),
                   Treated_Cali_sequential = sum(RAA1*Cweights_sequential),
                   Control_Cali_sequential = sum(RAA0*Cweights_sequential),
                   X1_treated_notcensored = sum(X1*A1)/sum(A1),
                   X1_control_notcensored = sum(X1*A0)/sum(A0),
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
                   min_IPW = min(weights),
                   min_Cali = min(Cweights),
                   min_Cali_sequential = min(Cweights_sequential),
                   max_IPW = max(weights),
                   max_Cali = max(Cweights),
                   max_Cali_sequential = max(Cweights_sequential))

switch_data$Cweights <- simdatafinal2$data[simdatafinal2$data$RA == 1,'Cweights']
switch_data$Cweights_sequential <- simdatafinal$data[simdatafinal$data$RA == 1,'Cweights']

PP_IPW <- glm(data = switch_data,
                                formula = Y ~ assigned_treatment + X1 + X2+ X3,
                                weights = weight, family = 'gaussian')

PP_calibrated <- glm(data = switch_data,
                    formula = Y ~ assigned_treatment + X1 + X2+ X3,
                    weights = Cweights, family = 'gaussian')

PP_calibrated_sequential <- glm(data = switch_data,
                    formula = Y ~ assigned_treatment + X1 + X2+ X3,
                    weights = Cweights_sequential, family = 'gaussian')

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
  dplyr::distinct(id, trial_period, followup_time, .keep_all = T) %>% 
  dplyr::group_by(id,trial_period) %>%
  dplyr::mutate(CA = cumsum(assigned_treatment)) %>%
  dplyr::filter(trial_period == 0) %>% 
  dplyr::ungroup()

fitting_data_control <- fitting_data_treatment %>%
  dplyr::mutate(assigned_treatment = 0.0)

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

predicted_probas_PP <- fitting_data_treatment %>%
  dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment_IPW,
                predicted_proba_control = Y_pred_PP_control_IPW,
                predicted_proba_treatment_cali = Y_pred_PP_treatment_cali,
                predicted_proba_control_cali = Y_pred_PP_control_cali,
                predicted_proba_treatment_cali_sequential = Y_pred_PP_treatment_cali_sequential,
                predicted_proba_control_cali_sequential = Y_pred_PP_control_cali_sequential) %>%
  dplyr::group_by(trial_period, followup_time) %>%
  dplyr::summarise(expected_treatment_IPW = mean(Y_pred_PP_treatment_IPW),
                   expected_control_IPW = mean(predicted_proba_control),
                   ate_IPW = expected_treatment_IPW - expected_control_IPW,
                   expected_treatment_cali = mean(predicted_proba_treatment_cali),
                   expected_control_cali = mean(predicted_proba_control_cali),
                   ate_cali = expected_treatment_cali - expected_control_cali,
                   expected_treatment_cali_sequential = mean(predicted_proba_treatment_cali_sequential),
                   expected_control_cali_sequential = mean(predicted_proba_control_cali_sequential),
                   ate_cali_sequential = expected_treatment_cali_sequential - expected_control_cali_sequential)


#########  Misspecified ############### 
simdata <- simdata %>% 
  mutate(Z1 = X1^3/9,Z2 = X1*X2, Z3 = log(abs(X3))+4)
PP_prep <- TrialEmulation::data_preparation(simdata, id='ID', period='t', treatment='A', outcome='Y', cense = 'C',
                                            eligible ='eligible',
                                            estimand_type = 'PP',
                                            switch_d_cov = ~Z1 + Z2 + Z3,
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
  dplyr::mutate(CA = cumsum(assigned_treatment))


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
                A1 = A_0,
                A0 = 1-A_0,
                sub = ID,
                tall = t) %>% 
  merge(dplyr::select(switch_data,id, followup_time, weight), 
        by.x = c('ID', 't'), by.y = c('id', 'followup_time'), all.x = T) %>% 
  dplyr::mutate(weights = ifelse(!is.na(weight), weight,0)) %>% 
  dplyr::arrange(ID, t) 

simdatafinal <- calibration_by_time(simdatafinal = data_restric, 
                                    var = c('A1','A1Z1', 'A1Z2','A1Z3',
                                            'A0','A0Z1','A0Z2','A0Z3'))
simdatafinal2 <- calibration(simdatafinal = data_restric, 
                             var = c('A1','A1Z1', 'A1Z2','A1Z3',
                                     'A0','A0Z1','A0Z2','A0Z3'))

meandiffs_summary <- simdatafinal$data %>% 
  dplyr::mutate(RAZ2 = RA*Z2,
                RAZ3 = RA*Z3,
                RAZ1 = RA*Z1,
                RAA1 = RA*A1,
                RAA0 = RA*A0,
                Cweights_sequential = Cweights,
                Cweights = simdatafinal2$data$Cweights) %>% 
  dplyr::select(t,RA,A1,A0,Z2,Z3,Z1,RAA1,RAA0, RAZ2, RAZ3, RAZ1, weights, Cweights,Cweights_sequential)


treatment_numbers <- meandiffs_summary %>% 
  dplyr::group_by(t) %>% 
  dplyr::summarise(Treated_notcensored = sum(A1),
                   Control_notcensored = sum(A0),
                   Treated_unadjusted = sum(RAA1),
                   Control_unadjusted = sum(RAA0),
                   Treated_IPW = sum(RAA1*weights),
                   Control_IPW = sum(RAA0*weights),
                   Treated_Cali = sum(RAA1 * Cweights),
                   Control_Cali = sum(RAA0*Cweights),
                   Treated_Cali_sequential = sum(RAA1*Cweights_sequential),
                   Control_Cali_sequential = sum(RAA0*Cweights_sequential),
                   Z1_treated = sum(RAZ1*A1)/sum(RAA1),
                   Z1_control = sum(RAZ1*A0)/sum(RAA0),
                   Z1_treated_IPW = sum(RAZ1*weights*A1)/sum(RAA1*weights),
                   Z1_control_IPW = sum(RAZ1*weights*A0)/sum(RAA0*weights),
                   Z1_treated_Cali_sequential = sum(RAZ1*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                   Z1_control_Cali_sequential = sum(RAZ1*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                   Z1_treated_Cali = sum(RAZ1*Cweights*A1)/sum(RAA1*Cweights),
                   Z1_control_Cali = sum(RAZ1*Cweights*A0)/sum(RAA0*Cweights),
                   Z2_treated = sum(RAZ2*A1)/sum(RAA1),
                   Z2_control = sum(RAZ2*A0)/sum(RAA0),
                   Z2_treated_IPW = sum(RAZ2*weights*A1)/sum(RAA1*weights),
                   Z2_control_IPW = sum(RAZ2*weights*A0)/sum(RAA0*weights),
                   Z2_treated_Cali_sequential = sum(RAZ2*Cweights_sequential*A1)/sum(RAA1*Cweights_sequential),
                   Z2_control_Cali_sequential = sum(RAZ2*Cweights_sequential*A0)/sum(RAA0*Cweights_sequential),
                   Z2_treated_Cali = sum(RAZ2*Cweights*A1)/sum(RAA1*Cweights),
                   Z2_control_Cali = sum(RAZ2*Cweights*A0)/sum(RAA0*Cweights),
                   min_IPW = min(weights),
                   min_Cali = min(Cweights),
                   min_Cali_sequential = min(Cweights_sequential),
                   max_IPW = max(weights),
                   max_Cali = max(Cweights),
                   max_Cali_sequential = max(Cweights_sequential))

switch_data$Cweights <- simdatafinal2$data[simdatafinal2$data$RA == 1,'Cweights']
switch_data$Cweights_sequential <- simdatafinal$data[simdatafinal$data$RA == 1,'Cweights']

PP_IPW <- glm(data = switch_data,
              formula = Y ~ assigned_treatment + Z1 + Z2+ Z3,
              weights = weight, family = 'gaussian')

PP_calibrated <- glm(data = switch_data,
                     formula = Y ~ assigned_treatment + Z1 + Z2+ Z3,
                     weights = Cweights, family = 'gaussian')

PP_calibrated_sequential <- glm(data = switch_data,
                                formula = Y ~ assigned_treatment + Z1 + Z2+ Z3,
                                weights = Cweights_sequential, family = 'gaussian')

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
  dplyr::distinct(id, trial_period, followup_time, .keep_all = T) %>% 
  dplyr::group_by(id,trial_period) %>%
  dplyr::mutate(CA = cumsum(assigned_treatment)) %>%
  dplyr::filter(trial_period == 0) %>% 
  dplyr::ungroup()

fitting_data_control <- fitting_data_treatment %>%
  dplyr::mutate(assigned_treatment = 0.0)

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

predicted_probas_PP <- fitting_data_treatment %>%
  dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment_IPW,
                predicted_proba_control = Y_pred_PP_control_IPW,
                predicted_proba_treatment_cali = Y_pred_PP_treatment_cali,
                predicted_proba_control_cali = Y_pred_PP_control_cali,
                predicted_proba_treatment_cali_sequential = Y_pred_PP_treatment_cali_sequential,
                predicted_proba_control_cali_sequential = Y_pred_PP_control_cali_sequential) %>%
  dplyr::group_by(trial_period, followup_time) %>%
  dplyr::summarise(expected_treatment_IPW = mean(Y_pred_PP_treatment_IPW),
                   expected_control_IPW = mean(predicted_proba_control),
                   ate_IPW = expected_treatment_IPW - expected_control_IPW,
                   expected_treatment_cali = mean(predicted_proba_treatment_cali),
                   expected_control_cali = mean(predicted_proba_control_cali),
                   ate_cali = expected_treatment_cali - expected_control_cali,
                   expected_treatment_cali_sequential = mean(predicted_proba_treatment_cali_sequential),
                   expected_control_cali_sequential = mean(predicted_proba_control_cali_sequential),
                   ate_cali_sequential = expected_treatment_cali_sequential - expected_control_cali_sequential)

