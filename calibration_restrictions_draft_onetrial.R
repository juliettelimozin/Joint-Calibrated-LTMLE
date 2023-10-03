#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
source("simulate_MSM_simplified.R")
library(TrialEmulation)
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)
library(cobalt)
source('calibration_func_trials.R')
simdata_censored<-DATA_GEN_censored_reduced(2500,5, conf = 1.5, censor = T) 

PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                            eligible ='eligible',
                                            switch_d_cov = ~X2 + X4,
                                            cense_d_cov = ~X2 + X4,
                                            outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                            use_weight=T, use_censor=T, quiet = T,
                                            save_weight_models = T,
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
                                    'A1X2', 'A0X2', 'A1nextX2', 'A0nextX2'))
simdatafinal$treat <- simdatafinal$A

bal.tab(simdatafinal[simdatafinal$t == 1,c('X2', 'X4')],treat = simdatafinal[simdatafinal$t == 1,]$A,
        stats = c('m', 'v'),var.name = c('X2', 'X4'))

bal.plot(simdatafinal[simdatafinal$t == 1,],stats = c('m', 'v'),var.name = c('X2', 'X4'),
         treat = simdatafinal[simdatafinal$t == 1,]$A )

bal.tab(simdatafinal[simdatafinal$t == 1,c('X2', 'X4')],treat = simdatafinal[simdatafinal$t == 1,]$A,
        stats = c('m', 'v'),var.name = c('X2', 'X4'), weights = simdatafinal[simdatafinal$t == 1,]$weights )

bal.plot(simdatafinal[simdatafinal$t == 1,],stats = c('m', 'v'),var.name = c('X2', 'X4'),
         treat = simdatafinal[simdatafinal$t == 1,]$A, weights = simdatafinal[simdatafinal$t == 1,]$weights )

bal.tab(simdatafinal[simdatafinal$t == 1,c('X2', 'X4')],treat = simdatafinal[simdatafinal$t == 1,]$A,
        stats = c('m', 'v'),var.name = c('X2', 'X4'), weights = simdatafinal[simdatafinal$t == 1,]$Cweights )
bal.plot(simdatafinal[simdatafinal$t == 1,],stats = c('m', 'v'),var.name = c('X2', 'X4'),
         treat = simdatafinal[simdatafinal$t == 1,]$A, weights = simdatafinal[simdatafinal$t == 1,]$Cweights )







########PREVIOUS VERSION (NOT SEAN'S CODE)#################################
lVEC <- 4*cbind(sum(switch_data[switch_data$followup_time == 0,]$assigned_treatment),
                sum(1-switch_data[switch_data$followup_time == 0,]$assigned_treatment),
                sum(switch_data[switch_data$followup_time == 1,]$assigned_treatment*switch_data[switch_data$followup_time == 1,]$X4),
                sum(1-switch_data[switch_data$followup_time == 1,]$assigned_treatment*switch_data[switch_data$followup_time == 1,]$X4),
                sum(switch_data[switch_data$followup_time == 1,]$assigned_treatment*switch_data[switch_data$followup_time == 1,]$X2),
                sum(1-switch_data[switch_data$followup_time == 1,]$assigned_treatment*switch_data[switch_data$followup_time == 1,]$X2),
                sum(switch_data[switch_data$followup_time == 1,]$assigned_treatment*switch_data[switch_data$followup_time == 1,]$X2),
                sum((1-switch_data[switch_data$followup_time == 1,]$assigned_treatment)*switch_data[switch_data$followup_time == 1,]$X2),
                0,
                0,
                sum(switch_data[switch_data$followup_time == 0,]$assigned_treatment),
                sum(1-switch_data[switch_data$followup_time == 0,]$assigned_treatment),
                sum(switch_data[switch_data$followup_time == 0,]$assigned_treatment*switch_data[switch_data$followup_time == 0,]$X4),
                sum((1-switch_data[switch_data$followup_time == 0,]$assigned_treatment)*switch_data[switch_data$followup_time == 0,]$X4),
                sum(switch_data[switch_data$followup_time == 0,]$assigned_treatment*switch_data[switch_data$followup_time == 0,]$X2),
                sum((1-switch_data[switch_data$followup_time == 0,]$assigned_treatment)*switch_data[switch_data$followup_time == 0,]$X2),
                0,
                0,
                0,
                0)

KMAT <-  cbind(data_restric$A*data_restric$RA,
               (1-data_restric$A)*data_restric$RA,
               data_restric$A*data_restric$RA*data_restric$X4,
               (1-data_restric$A)*data_restric$RA*data_restric$X4,
               data_restric$A*data_restric$RA*((4 - data_restric$t - 1)*data_restric$X2 - (4 - data_restric$t)*data_restric$nextX2),
               (1-data_restric$A)*data_restric$RA*((4 - data_restric$t - 1)*data_restric$X2 - (4 - data_restric$t)*data_restric$nextX2),
               data_restric$A*data_restric$RA*((4 - 1 - 1)*data_restric$X2*data_restric$t1 - (4 - 1)*data_restric$nextX2*data_restric$t1),
               (1-data_restric$A)*data_restric$RA*((4 - 1 - 1)*data_restric$X2*data_restric$t1 - (4 - 1)*data_restric$nextX2*data_restric$t1),
               data_restric$A*data_restric$RA*((4 - 2 - 1)*data_restric$X2*data_restric$t2 - (4 - 1)*data_restric$nextX2*data_restric$t2),
               (1-data_restric$A)*data_restric$RA*((4 - 2 - 1)*data_restric$X2*data_restric$t2 - (4 - 2)*data_restric$nextX2*data_restric$t2),
               data_restric$A,
               (1-data_restric$A),
               data_restric$A*data_restric$X4,
               (1-data_restric$A)*data_restric$X4,
               data_restric$A*((4 - data_restric$t - 1)*data_restric$prevX2 - (4 - data_restric$t)*data_restric$X2),
               (1-data_restric$A)*((4 - data_restric$t - 1)*data_restric$prevX2 - (4 - data_restric$t)*data_restric$X2),
               data_restric$A*data_restric$RA*((4 - 1 - 1)*data_restric$prevX2*data_restric$t1 - (4 -1)*data_restric$X2*data_restric$t1),
               (1-data_restric$A)*data_restric$RA*((4 - 1 - 1)*data_restric$prevX2*data_restric$t1 - (4 - 1)*data_restric$X2*data_restric$t1),
               data_restric$A*data_restric$RA*((4 - 1 - 1)*data_restric$prevX2*data_restric$t2 - (4 - 1)*data_restric$X2*data_restric$t2),
               (1-data_restric$A)*data_restric$RA*((4 - 2 - 1)*data_restric$prevX2*data_restric$t2 - (4 - 2)*data_restric$X2*data_restric$t2))


KMAT_old <-  cbind(data_restric$A*data_restric$RA,
               (1-data_restric$A)*data_restric$RA,
               data_restric$A*data_restric$RA*data_restric$X4,
               (1-data_restric$A)*data_restric$RA*data_restric$X4,
               data_restric$A*data_restric$RA*((4 - data_restric$t - 1)*data_restric$X2 - (4 - data_restric$t)*data_restric$nextX2),
               (1-data_restric$A)*data_restric$RA*((4 - data_restric$t - 1)*data_restric$X2 - (4 - data_restric$t)*data_restric$nextX2),
               data_restric$A*(1-data_restric$C),
               (1-data_restric$A)*(1-data_restric$C),
               data_restric$A*data_restric$X4*(1-data_restric$C),
               (1-data_restric$A)*data_restric$X4*(1-data_restric$C),
               data_restric$A*(1-data_restric$C)*((4 - data_restric$t - 1)*data_restric$prevX2 - (4 - data_restric$t)*data_restric$X2),
               (1-data_restric$A)*(1-data_restric$C)*((4 - data_restric$t - 1)*data_restric$prevX2 - (4 - data_restric$t)*data_restric$X2))

constraints <- function(w){
  t(KMAT)%*%(switch_data[switch_data$followup_time != 0,]$weight*exp(KMAT%*%w)) - t(lVEC)
}


cali_restriction_check <- function(w){
  t(KMAT)%*%w - t(lVEC)
}

library(nleqslv)
weioptAR<-nleqslv(rep(0,20),constraints,method="Broyden",
                  control=list(maxit=50000,ftol=10^(-16),xtol=10^(-16), cndtol = 10^(-15), allowSingular = TRUE))
calibrated_weights <- cbind(data_restric$ID, data_restric$t, switch_data[switch_data$followup_time != 0,]$weight*exp(KMAT%*%weioptAR$x))
colnames(calibrated_weights) <- c('id', 'followup_time', 'calibrated_weight')

switch_data <- switch_data %>% 
  merge(calibrated_weights, by = c('id', 'followup_time'), all.x = T) %>% 
  dplyr::mutate(assigned_treatment = as.factor(assigned_treatment),
                calibrated_weight = ifelse(!is.na(calibrated_weight), calibrated_weight, 1))

eligible_data <- simdata_censored %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1), nextX2 = lead(X2)) %>% 
  dplyr::mutate(CRA = cumsum(RA),
                nextX2 = ifelse(is.na(nextX2), 0, nextX2)) %>% 
  dplyr::filter(CRA == t+1) %>% 
  dplyr::ungroup() %>% 
  merge(dplyr::select(switch_data,id, followup_time, weight, calibrated_weight), 
        by.x = c('ID', 't'), by.y = c('id', 'followup_time')) %>% 
  dplyr::arrange(ID, t)

print(cali_restriction_check(eligible_data[eligible_data$t != 0,]$weight))
print(cali_restriction_check(eligible_data[eligible_data$t != 0,]$calibrated_weight))

bal.plot(eligible_data[eligible_data$t == 1,],stats = c('m', 'v'),var.name = c('X2', 'X4'),
         treat = eligible_data[eligible_data$t == 1,]$A_0 )


bal.tab(eligible_data[eligible_data$t == 1,],stats = c('m', 'v'),
        treat = eligible_data[eligible_data$t == 1,]$A_0, weights = eligible_data[eligible_data$t == 1,]$weight )

bal.plot(eligible_data[eligible_data$t == 1,],stats = c('m', 'v'),var.name = c('X2', 'X4'),
         treat = eligible_data[eligible_data$t == 1,]$A_0, weights = eligible_data[eligible_data$t == 1,]$weight )


bal.tab(eligible_data[eligible_data$t == 1,],stats = c('m', 'v'),
        treat = eligible_data[eligible_data$t == 1,]$A_0, weights = eligible_data[eligible_data$t == 1,]$calibrated_weight )
bal.plot(eligible_data[eligible_data$t == 1,],stats = c('m', 'v'),var.name = c('X2', 'X4'),
        treat = eligible_data[eligible_data$t == 1,]$A_0, weights = eligible_data[eligible_data$t == 1,]$calibrated_weight )




# PP <- TrialEmulation::data_modelling(data = switch_data,
#                                      outcome_cov = ~ X2 + X4+ assigned_treatment+
#                                        t_1 + t_2 + t_3 + t_4 +
#                                        t_1A + t_2A + t_3A + t_4A + 
#                                        t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
#                                        t_1X4 + t_2X4 + t_3X4 + t_4X4,
#                                      model_var = c('assigned_treatment'),
#                                      glm_function = 'glm',
#                                      include_expansion_time = ~1, include_followup_time = ~1,
#                                      use_weight=1, use_censor=1, quiet = T, use_sample_weights =  F)
# switch_data$p_i <- predict.glm(PP$model, switch_data,type = 'response')
# 
# switch_d0 <- readRDS(paste(data_direction,'/weight_model_switch_d0.rds', sep = ""))
# switch_n0 <- readRDS(paste(data_direction,'/weight_model_switch_n0.rds', sep = ""))
# switch_d1 <- readRDS(paste(data_direction,'/weight_model_switch_d1.rds', sep = ""))
# switch_n1 <- readRDS(paste(data_direction,'/weight_model_switch_n1.rds', sep = ""))
# 
# design_mat <- expand.grid(id = 1:as.numeric(scenarios[l,1]),
#                           trial_period = 0:4,
#                           followup_time = 0:4) 
# design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]
# 
# fitting_data_treatment <-  switch_data %>% 
#   dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
#   dplyr::select(id,trial_period, followup_time, X2,  X4, assigned_treatment) %>% 
#   merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>% 
#   dplyr::group_by(id) %>% 
#   tidyr::fill( X2,X4,assigned_treatment,.direction = "down") %>% 
#   dplyr::ungroup() %>% 
#   dplyr::select(id, trial_period, followup_time, X2, X4, assigned_treatment) %>% 
#   merge(data.frame(id = switch_data$id, trial_period = switch_data$trial_period), by = c("id", "trial_period"), all.y = TRUE) %>% 
#   dplyr::arrange(id, trial_period, followup_time) %>% 
#   dplyr::mutate(t_1 = ifelse(followup_time == 1,1,0),
#                 t_2 = ifelse(followup_time == 2,1,0),
#                 t_3 = ifelse(followup_time == 3,1,0),
#                 t_4 = ifelse(followup_time == 4,1,0),
#                 t_1A = t_1*assigned_treatment,
#                 t_2A = t_2*assigned_treatment,
#                 t_3A = t_3*assigned_treatment,
#                 t_4A = t_4*assigned_treatment,
#                 t_1X2 = t_1*X2,
#                 t_2X2 = t_2*X2,
#                 t_3X2 = t_3*X2,
#                 t_4X2 = t_4*X2,
#                 t_1X4 = t_1*X4,
#                 t_2X4 = t_2*X4,
#                 t_3X4 = t_3*X4,
#                 t_4X4 = t_4*X4) %>% 
#   dplyr::filter(trial_period == 0)
# 
# fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]
# 
# fitting_data_control <- fitting_data_treatment %>% 
#   dplyr::mutate(assigned_treatment = assigned_treatment*0,
#                 t_1A = t_1*0,
#                 t_2A = t_2*0,
#                 t_3A = t_3*0,
#                 t_4A = t_4*0)
# 
# Y_pred_PP_treatment <- predict.glm(PP$model, 
#                                    fitting_data_treatment, 
#                                    type = "response")
# Y_pred_PP_control <- predict.glm(PP$model, 
#                                  fitting_data_control,
#                                  type = "response")
# predicted_probas_PP <- fitting_data_treatment %>% 
#   dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment,
#                 predicted_proba_control = Y_pred_PP_control) %>% 
#   dplyr::group_by(id, trial_period) %>% 
#   dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
#                 cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::group_by(followup_time) %>% 
#   dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
#                    survival_control = mean(cum_hazard_control),
#                    survival_difference = survival_treatment - survival_control)

estimates[,i] <- pull(predicted_probas_PP,survival_difference)