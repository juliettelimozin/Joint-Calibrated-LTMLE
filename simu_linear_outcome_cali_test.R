#!/usr/bin R
library(dplyr)
library(tidyr)
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("continuous_timevarying_outcome_datagen.R")
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
size <- c(200,500,1000,5000)

# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 10)
multiResultClass <- function(weights = NULL,
                             hr_estimates=NULL,predict_estimates = NULL)
{
  me <- list(
    weights = weights,
    hr_estimates = hr_estimates,
    predict_estimates = predict_estimates
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}
for (l in 1:4){
  oper <- foreach(i = 1:iters,.combine=cbind) %dopar% {
    tryCatch({
      result <- multiResultClass()
      
      simdata<-DATA_GEN_continous_timevarying_outcome_treatment_switch(ns = as.numeric(size[l]),nv = 5,treat_prev =  0.5,
                                                                       conf =  1.5,
                                                                       censor = F)
      simdata <- simdata %>% 
        mutate(switch = ifelse(t == 0, 0,ifelse(A!=Ap,1,0)),
               missCX1 = log(abs(CX1))/4,
               missCX2 = sqrt(abs(CX2))/3,
               missX1 = log(abs(X1))/4,
               missX2 = sqrt(abs(X2))/3,
               CX1sq = CX1^2,
               CX2sq = CX2^2,
               X1sq = X1^2,
               X2sq = X2^2)
      #con4<-xtabs(~t + switch, data=simdata)
      #ftable(con4)
      
      
      
      ######### Correctly specified ############### 
      PP_prep <- TrialEmulation::data_preparation(simdata, id='ID', period='t', treatment='A', outcome='Y', cense = 'C',
                                                  eligible ='eligible',
                                                  estimand_type = 'PP',
                                                  switch_d_cov = ~ X1 + X2,
                                                  switch_n_cov = ~ -1,
                                                  outcome_cov = ~ X1 + X2, model_var = c('assigned_treatment'),
                                                  quiet = T,
                                                  save_weight_models = F,
                                                  data_dir = getwd())
      
      switch_data <- PP_prep$data %>% 
        dplyr::mutate(t = trial_period + followup_time) %>%
        merge(simdata[,c('ID', 't', 'Y')], by.x = c('id', 't'), by.y = c('ID', 't')) %>% 
        dplyr::filter(trial_period == 0) %>% 
        dplyr::arrange(id, followup_time) %>% 
        dplyr::group_by(id) %>% 
        dplyr::mutate(Ap = ifelse(followup_time == 0, 0,lag(assigned_treatment)),
                      CAp = cumsum(Ap),
                      CA = cumsum(assigned_treatment))
      
      #con4<-xtabs(~followup_time + assigned_treatment, data=switch_data)
      #ftable(con4)
      
      ########## Outcome imputation correct and miss ############
      library(data.table)
      wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, value.var = c("A", "X1", "X2", "missX1", "missX2", "Y"))
      
      wideSimdata$y4pred <- wideSimdata$Y_4
      
      q4 <- glm(y4pred ~ X1_4 + X2_4 + A_4,
                 data = wideSimdata, family = 'gaussian')
      q4miss <- glm(y4pred ~ missX1_4 + missX2_4  + A_4,
                 data = wideSimdata, family = 'gaussian')
      
      wideSimdata$y4pred <- predict.glm(q4,type = "response", newdata = wideSimdata %>% mutate(A_4 ))
      wideSimdata$y4predmiss <- predict.glm(q4miss,type = "response", newdata = wideSimdata %>% mutate(A_4 = A_0))
      
      q3 <- glm(y4pred ~ X1_3 + X2_3 + A_3,
                 data = wideSimdata, family = 'gaussian')
      q3miss <- glm(y4predmiss ~ missX1_3 + missX2_3  + A_3,
                data = wideSimdata, family = 'gaussian')
      
      wideSimdata$y4pred <- predict.glm(q3,type = "response", newdata = wideSimdata %>% mutate(A_3 = A_0))
      wideSimdata$y4predmiss <- predict.glm(q3miss,type = "response", newdata = wideSimdata %>% mutate(A_3 = A_0))
      
      q2 <- glm(y4pred ~  X1_2 + X2_2 + A_2,
                 data = wideSimdata, family = 'gaussian')
      q2miss <- glm(y4predmiss ~  missX1_2 + missX2_2 + A_2,
                data = wideSimdata, family = 'gaussian')
      
      wideSimdata$y4pred <- predict.glm(q2,type = "response", newdata = wideSimdata %>% mutate(A_2 = A_0))
      wideSimdata$y4predmiss <- predict.glm(q2miss,type = "response", newdata = wideSimdata %>% mutate(A_2 = A_0))
      
      q1 <- glm(y4pred ~  X1_1 + X2_1 + A_1,
                data = wideSimdata, family = 'gaussian')
      q1miss <- glm(y4predmiss ~  missX1_1 + missX2_1  + A_1,
                data = wideSimdata, family = 'gaussian')
      
      wideSimdata$y4pred <- predict.glm(q1,type = "response", newdata = wideSimdata %>% mutate(A_1 = A_0))
      wideSimdata$y4predmiss <- predict.glm(q1miss,type = "response", newdata = wideSimdata %>% mutate(A_1 = A_0))
      
      q0 <- glm(y4pred ~  X1_0 + X2_0 + A_0,
                data = wideSimdata, family = 'gaussian')
      q0miss <- glm(y4predmiss ~  missX1_0 + missX2_0 + A_0,
                data = wideSimdata, family = 'gaussian')
      
      wideSimdata$y4pred <- predict.glm(q0,type = "response", newdata = wideSimdata)
      wideSimdata$y4predmiss <- predict.glm(q0miss,type = "response", newdata = wideSimdata)
      
      wideSimdata$Y0 <- mean(wideSimdata$A_0*predict.glm(q0, type = "response", newdata = wideSimdata %>% mutate(A_0 = 1))) + mean((1-wideSimdata$A_0)*predict.glm(q0, type = "response", newdata = wideSimdata %>% mutate(A_0 = 0)))
      wideSimdata$Y1 <- mean(wideSimdata$A_1*predict.glm(q1, type = "response", newdata = wideSimdata %>% mutate(A_1 = 1))) + mean((1-wideSimdata$A_1)*predict.glm(q1, type = "response", newdata = wideSimdata %>% mutate(A_1 = 0)))
      wideSimdata$Y2 <- mean(wideSimdata$A_2*predict.glm(q2, type = "response", newdata = wideSimdata %>% mutate(A_2 = 1))) + mean((1-wideSimdata$A_2)*predict.glm(q2, type = "response", newdata = wideSimdata %>% mutate(A_2 = 0)))
      wideSimdata$Y3 <- mean(wideSimdata$A_3*predict.glm(q3, type = "response", newdata = wideSimdata %>% mutate(A_3 = 1))) + mean((1-wideSimdata$A_3)*predict.glm(q3, type = "response", newdata = wideSimdata %>% mutate(A_3 = 0)))
      wideSimdata$Y4 <- mean(wideSimdata$A_4*predict.glm(q4, type = "response", newdata = wideSimdata %>% mutate(A_4 = 1))) + mean((1-wideSimdata$A_4)*predict.glm(q4, type = "response", newdata = wideSimdata %>% mutate(A_4 = 0)))
      
      
      wideSimdata$missY0 <- mean(wideSimdata$A_0*predict.glm(q0miss, type = "response", newdata = wideSimdata %>% mutate(A_0 = 1))) + mean((1-wideSimdata$A_0)*predict.glm(q0miss, type = "response", newdata = wideSimdata %>% mutate(A_0 = 0)))
      wideSimdata$missY1 <- mean(wideSimdata$A_1*predict.glm(q1miss, type = "response", newdata = wideSimdata %>% mutate(A_1 = 1))) + mean((1-wideSimdata$A_1)*predict.glm(q1miss, type = "response", newdata = wideSimdata %>% mutate(A_1 = 0)))
      wideSimdata$missY2 <- mean(wideSimdata$A_2*predict.glm(q2miss, type = "response", newdata = wideSimdata %>% mutate(A_2 = 1))) + mean((1-wideSimdata$A_2)*predict.glm(q2miss, type = "response", newdata = wideSimdata %>% mutate(A_2 = 0)))
      wideSimdata$missY3 <- mean(wideSimdata$A_3*predict.glm(q3miss, type = "response", newdata = wideSimdata %>% mutate(A_3 = 1))) + mean((1-wideSimdata$A_3)*predict.glm(q3miss, type = "response", newdata = wideSimdata %>% mutate(A_3 = 0)))
      wideSimdata$missY4 <- mean(wideSimdata$A_4*predict.glm(q4miss, type = "response", newdata = wideSimdata %>% mutate(A_4 = 1))) + mean((1-wideSimdata$A_4)*predict.glm(q4miss, type = "response", newdata = wideSimdata %>% mutate(A_4 = 0)))
      
      
      simdata_imputed <- pivot_longer(wideSimdata %>% dplyr::select(ID, Y0, Y1, Y2, Y3, Y4),
                                   cols = Y0:Y4, names_to = 't', names_prefix = "Y", values_to = "Yhat") %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, missY0, missY1, missY2, missY3, missY4),
                           cols = missY0:missY4, names_to = 't', names_prefix = "missY", values_to = "missYhat"),by = c("ID", "t")) %>% 
        dplyr::mutate(t = as.numeric(t))
      
    
      ########### Data prep for calibration ###############
      outcome_model <- glm(Y ~  X1 + X2 + A,
                           data = simdata, family = 'gaussian')
      outcome_model_miss <- glm(Y ~  missX1 + missX2 + A,
                           data = simdata, family = 'gaussian')
      
      simdata$Ythat <- mean(simdata$A*predict.glm(outcome_model, type = "response", newdata = simdata %>% mutate(A =1))) + mean((1-simdata$A)*predict.glm(outcome_model, type = "response", newdata = simdata %>% mutate(A =0)))
      simdata$Ytmisshat <- mean(simdata$A*predict.glm(outcome_model_miss, type = "response", newdata = simdata %>% mutate(A =1))) + mean((1-simdata$A)*predict.glm(outcome_model_miss, type = "response", newdata = simdata %>% mutate(A =0)))
      
      data_restric <- simdata %>% 
        merge(simdata_imputed, by = c("ID", "t")) %>% 
        dplyr::group_by(ID) %>% 
        dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1)) %>% 
        dplyr::mutate(CRA = cumsum(RA),
                      RA = ifelse(CRA == t+1,1,0),
                      RC = ifelse(lag(C) == 0,1,0),
                      A1X2 = A_0*X2,
                      A0X2 = (1-A_0)*X2,
                      A1X1 = A_0*X1,
                      A0X1 = (1-A_0)*X1,
                      A1 = A_0,
                      A0 = 1-A_0,
                      A1Yhat = A_0*Yhat,
                      A0Yhat = (1-A_0)*Yhat,
                      A1missYhat = A_0*missYhat,
                      A0missYhat = (1-A_0)*missYhat,
                      A1Ythat = A_0*Ythat,
                      A0Ythat = (1-A_0)*Ythat,
                      A1Ytmisshat = A_0*Ytmisshat,
                      A0Ytmisshat = (1-A_0)*Ytmisshat,
                      sub = ID,
                      tall = t,
                      One = 1.0) %>% 
        merge(dplyr::select(switch_data,id, followup_time, weight), 
              by.x = c('ID', 't'), by.y = c('id', 'followup_time'), all.x = T) %>% 
        dplyr::mutate(weights = ifelse(!is.na(weight), weight, 0)) %>% 
        dplyr::arrange(ID, t) 
    
      ############# IPW estimation ################
      weight_model <- glm(data = simdata, 
                          formula = A ~ Ap + X1 + X2, family = 'binomial')
      weight_model_miss <- glm(data = simdata, 
                          formula = A ~ Ap + X1, family = 'binomial')
      summary(weight_model)
      summary(weight_model_miss)
      data_restric$p_1 <- 1.0
      data_restric$p_1miss <- 1.0
      
      data_restric[data_restric$t != 0,]$p_1 <- predict.glm(weight_model, data_restric[data_restric$t != 0,], type = 'response')
      data_restric[data_restric$t != 0,]$p_1miss <- predict.glm(weight_model_miss, data_restric[data_restric$t != 0,], type = 'response')
      
      data_restric <- data_restric %>% 
        arrange(ID, t) %>% 
        group_by(ID) %>% 
        dplyr::mutate(
          truep_1 = 1/(1+exp(-(2.5*Ap -2.5*(1-Ap) + 1.5*(X1+X2)))),
          truewt = ifelse( t == 0, 1.0,ifelse(A == 1, 1/truep_1, 1/(1-truep_1))),
          truewtprod = cumprod(truewt),
          true_weights = ifelse(weights !=0.0,truewtprod,0.0),
          wt = ifelse( t == 0, 1.0,ifelse(A == 1, 1/p_1, 1/(1-p_1))),
          wtprod = cumprod(wt),
          weights = ifelse(weights !=0.0,wtprod,0.0),
          wtmiss = ifelse( t == 0, 1.0,ifelse(A == 1, 1/p_1miss, 1/(1-p_1miss))),
          wtprodmiss = cumprod(wtmiss),
          weights_miss = ifelse(weights !=0.0,wtprodmiss,0.0))
    
      ################### Calibration by time Lk #######################
      simdatafinal1 <- calibration(simdatafinal = data_restric, 
                                           var = c('A1','A1X1', 'A1X2',
                                                   'A0','A0X1','A0X2'))
      
      

      ################### Calibration by time g(Lk) gimp #######################
      simdatafinal2 <- calibration(simdatafinal = data_restric, 
                                           var = c('A1','A1Yhat',
                                                   'A0','A0Yhat'))
      
      
      ################## Calibration by time Lk, g(Lk) gimp###########################
      simdatafinal3 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1X1', 'A1X2','A1Yhat',
                                                   'A0','A0X1','A0X2','A0Yhat'))
      
      
      
      ################## Calibration by time g(Lk) gimp miss #######################
      simdatafinal4 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1missYhat',
                                                   'A0','A0missYhat'))
      
      ################## Calibration by time Lk, g(Lk) gimp ###########################
      simdatafinal5 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1X1', 'A1X2','A1missYhat',
                                                   'A0','A0X1','A0X2','A0missYhat'))
      
      ################### Calibration by time Lk miss #######################
      simdatafinal6 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1X1',
                                                   'A0','A0X1'))
      ################### Calibration by time miss g(Lk) gimp correct #######################
      simdatafinal7 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1Yhat',
                                                   'A0','A0Yhat'))
      ################### Calibration by time miss Lk g(Lk)gimp correct #######################
      simdatafinal8 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1X1', 'A1Yhat',
                                                   'A0','A0X1', 'A0Yhat'))
      
      ################### Calibration by time all gimp miss #######################
      simdatafinal9 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1X1', 'A1missYhat',
                                                   'A0','A0X1', 'A0missYhat'))
      
    
      
      ################### Calibration by time g(Lk) basic #######################
      simdatafinal10 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1Ythat',
                                                   'A0','A0Ythat'))
      
      
      ################## Calibration by time Lk, g(Lk) basic###########################
      simdatafinal11 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1X1', 'A1X2','A1Ythat',
                                                   'A0','A0X1','A0X2','A0Ythat'))
      
      
      
      ################## Calibration by time g(Lk) basic miss #######################
      simdatafinal12 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1Ytmisshat',
                                                   'A0','A0Ytmisshat'))
      
      ################## Calibration by time Lk, g(Lk) basic miss  ###########################
      simdatafinal13 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1X1', 'A1X2','A1Ytmisshat',
                                                   'A0','A0X1','A0X2','A0Ytmisshat'))
      
  
      ################### Calibration by time miss g(Lk) gimp correct #######################
      simdatafinal14 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1Ythat',
                                                   'A0','A0Ythat'))
     
      
      ################### Calibration by time miss Lk g(Lk)gimp correct #######################
      simdatafinal15 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1X1', 'A1Ythat',
                                                   'A0','A0X1', 'A0Ythat'))
      ################### Calibration by time all gimp miss #######################
      simdatafinal16 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1X1', 'A1Ytmisshat',
                                                   'A0','A0X1', 'A0Ytmisshat'))
      
      
      
      switch_data$IPWcorrect <- simdatafinal1$data[simdatafinal1$data$RA == 1,]$weights
      switch_data$CaliLkcorrect <- simdatafinal1$data[simdatafinal1$data$RA == 1,]$Cweights
      switch_data$CaliGLkcorrect <- simdatafinal2$data[simdatafinal2$data$RA == 1,]$Cweights
      switch_data$CaliAllCorrect <- simdatafinal3$data[simdatafinal3$data$RA == 1,]$Cweights
      switch_data$CaliGLkmiss <- simdatafinal4$data[simdatafinal4$data$RA == 1,]$Cweights
      switch_data$CaliAllGLkmiss <- simdatafinal5$data[simdatafinal5$data$RA == 1,]$Cweights
      switch_data$CaliLkmiss <- simdatafinal6$data[simdatafinal6$data$RA == 1,]$Cweights
      switch_data$CaliMissGLkcorrect <- simdatafinal7$data[simdatafinal7$data$RA == 1,]$Cweights
      switch_data$CaliMissAllGLkcorrect <- simdatafinal8$data[simdatafinal8$data$RA == 1,]$Cweights
      switch_data$CaliAllMiss <- simdatafinal9$data[simdatafinal9$data$RA == 1,]$Cweights
      switch_data$IPWmiss <- simdatafinal9$data[simdatafinal9$data$RA == 1,]$weights
      switch_data$CaliGLkcorrectbasic <- simdatafinal10$data[simdatafinal10$data$RA == 1,]$Cweights
      switch_data$CaliAllCorrectbasic <- simdatafinal11$data[simdatafinal11$data$RA == 1,]$Cweights
      switch_data$CaliGLkmissbasic <- simdatafinal12$data[simdatafinal12$data$RA == 1,]$Cweights
      switch_data$CaliAllGLkmissbasic <- simdatafinal13$data[simdatafinal13$data$RA == 1,]$Cweights
      switch_data$CaliMissGLkcorrectbasic <- simdatafinal14$data[simdatafinal14$data$RA == 1,]$Cweights
      switch_data$CaliMissAllGLkcorrectbasic <- simdatafinal15$data[simdatafinal15$data$RA == 1,]$Cweights
      switch_data$CaliAllMissbasic <- simdatafinal16$data[simdatafinal16$data$RA == 1,]$Cweights

      MSM_data  <- switch_data %>% 
        filter(t == 4)
      
      PP_naive <- glm(data = MSM_data,
                      formula = Y ~ CA,
                      weights = NULL, family = 'gaussian')
      summary(PP_naive)
      
      PP_IPWcorrect <- glm(data = MSM_data,
                    formula = Y ~ CA,
                    weights = IPWcorrect, family = 'gaussian')
      summary(PP_IPWcorrect)
      
      PP_CaliLkcorrect <- glm(data = MSM_data,
                           formula = Y ~ CA,
                           weights = CaliLkcorrect, family = 'gaussian')
      summary(PP_CaliLkcorrect)
      
      PP_CaliGLkcorrect <- glm(data = MSM_data,
                                      formula = Y ~ CA,
                                      weights = CaliGLkcorrect, family = 'gaussian')
      summary(PP_CaliGLkcorrect)
      
      PP_CaliAllCorrect <- glm(data = MSM_data,
                             formula = Y ~ CA,
                             weights = CaliAllCorrect, family = 'gaussian')
      summary(PP_CaliAllCorrect)
      
      PP_CaliGLkmiss <- glm(data = MSM_data,
                               formula = Y ~ CA,
                               weights = CaliGLkmiss, family = 'gaussian')
      summary(PP_CaliGLkmiss)
      
      PP_CaliAllGLkmiss <- glm(data = MSM_data,
                 formula = Y ~ CA,
                 weights = CaliAllGLkmiss, family = 'gaussian')
      summary(PP_CaliAllGLkmiss)
      
      PP_IPWmiss <- glm(data = MSM_data,
                 formula = Y ~ CA,
                 weights =IPWmiss , family = 'gaussian')
      summary(PP_IPWmiss)
      
      PP_CaliLkmiss <- glm(data = MSM_data,
                 formula = Y ~ CA,
                 weights = CaliLkmiss, family = 'gaussian')
      summary(PP_CaliLkmiss)
      
      PP_CaliMissGLkcorrect <- glm(data = MSM_data,
                 formula = Y ~ CA,
                 weights = CaliMissGLkcorrect, family = 'gaussian')
      summary(PP_CaliMissGLkcorrect)
      
      PP_CaliMissAllGLkcorrect <- glm(data = MSM_data,
                 formula = Y ~ CA,
                 weights = CaliMissAllGLkcorrect, family = 'gaussian')
      summary(PP_CaliMissAllGLkcorrect)
      
      PP_CaliAllMiss <- glm(data = MSM_data,
                 formula = Y ~ CA,
                 weights = CaliAllMiss, family = 'gaussian')
      summary(PP_CaliAllMiss)
      
      
     
      PP_CaliGLkcorrectbasic <- glm(data = MSM_data,
                               formula = Y ~ CA,
                               weights = CaliGLkcorrectbasic, family = 'gaussian')
      summary(PP_CaliGLkcorrectbasic)
      
      PP_CaliAllCorrectbasic <- glm(data = MSM_data,
                               formula = Y ~ CA,
                               weights = CaliAllCorrectbasic, family = 'gaussian')
      summary(PP_CaliAllCorrectbasic)
      
      PP_CaliGLkmissbasic <- glm(data = MSM_data,
                            formula = Y ~ CA,
                            weights = CaliGLkmissbasic, family = 'gaussian')
      summary(PP_CaliGLkmissbasic)
      
      PP_CaliAllGLkmissbasic <- glm(data = MSM_data,
                               formula = Y ~ CA,
                               weights = CaliAllGLkmissbasic, family = 'gaussian')
      summary(PP_CaliAllGLkmissbasic)
      
      
      PP_CaliMissGLkcorrectbasic <- glm(data = MSM_data,
                                   formula = Y ~ CA,
                                   weights = CaliMissGLkcorrectbasic, family = 'gaussian')
      summary(PP_CaliMissGLkcorrectbasic)
      
      PP_CaliMissAllGLkcorrectbasic <- glm(data = MSM_data,
                                      formula = Y ~ CA,
                                      weights = CaliMissAllGLkcorrectbasic, family = 'gaussian')
      summary(PP_CaliMissAllGLkcorrectbasic)
      
      PP_CaliAllMissbasic <- glm(data = MSM_data,
                            formula = Y ~ CA,
                            weights = CaliAllMissbasic, family = 'gaussian')
      summary(PP_CaliAllMissbasic)
      

      result$hr_estimates<- array(,dim = c(19,2))
      result$hr_estimates[1,] <- PP_naive$coefficients
      result$hr_estimates[2,] <- PP_IPWcorrect$coefficients
      result$hr_estimates[3,] <- PP_CaliLkcorrect$coefficients
      result$hr_estimates[4,] <- PP_CaliGLkcorrect$coefficients
      result$hr_estimates[5,] <- PP_CaliAllCorrect$coefficients
      result$hr_estimates[6,] <- PP_CaliGLkmiss$coefficients
      result$hr_estimates[7,] <- PP_CaliAllGLkmiss$coefficients
      result$hr_estimates[8,] <- PP_IPWmiss$coefficients
      result$hr_estimates[9,] <- PP_CaliLkmiss$coefficients
      result$hr_estimates[10,] <- PP_CaliMissGLkcorrect$coefficients
      result$hr_estimates[11,] <- PP_CaliMissAllGLkcorrect$coefficients
      result$hr_estimates[12,] <- PP_CaliAllMiss$coefficients
      result$hr_estimates[13,] <- PP_CaliGLkcorrectbasic$coefficients
      result$hr_estimates[14,] <- PP_CaliAllCorrectbasic$coefficients
      result$hr_estimates[15,] <- PP_CaliGLkmissbasic$coefficients
      result$hr_estimates[16,] <- PP_CaliAllGLkmissbasic$coefficients
      result$hr_estimates[17,] <- PP_CaliMissGLkcorrectbasic$coefficients
      result$hr_estimates[18,] <- PP_CaliMissAllGLkcorrectbasic$coefficients
      result$hr_estimates[19,] <- PP_CaliAllMissbasic$coefficients

      
      ATE <- data.frame(k = c(4)) %>% 
        mutate(ATE_naive = predict.glm(PP_naive, data.frame(CA = 5))- predict.glm(PP_naive, data.frame(CA = 0)),
               ATE_IPWcorrect = predict.glm(PP_IPWcorrect, data.frame(CA = 5))- predict.glm(PP_IPWcorrect, data.frame(CA = 0)),
               ATE_CaliLkcorrect = predict.glm(PP_CaliLkcorrect, data.frame(CA = 5))- predict.glm(PP_CaliLkcorrect, data.frame(CA = 0)),
               ATE_CaliGLkcorrect = predict.glm(PP_CaliGLkcorrect, data.frame(CA = 5))- predict.glm(PP_CaliGLkcorrect, data.frame(CA = 0)),
               ATE_CaliAllCorrect = predict.glm(PP_CaliAllCorrect, data.frame(CA = 5))- predict.glm(PP_CaliAllCorrect, data.frame(CA = 0)),
               ATE_CaliGLkmiss = predict.glm(PP_CaliGLkmiss, data.frame(CA = 5))- predict.glm(PP_CaliGLkmiss, data.frame(CA = 0)),
               ATE_CaliAllGLkmiss = predict.glm(PP_CaliAllGLkmiss, data.frame(CA = 5))- predict.glm(PP_CaliAllGLkmiss, data.frame(CA = 0)),
               ATE_IPWmiss = predict.glm(PP_IPWmiss, data.frame(CA = 5))- predict.glm(PP_IPWmiss, data.frame(CA = 0)),
               ATE_CaliLkmiss = predict.glm(PP_CaliLkmiss, data.frame(CA = 5))- predict.glm(PP_CaliLkmiss, data.frame(CA = 0)),
               ATE_CaliMissGLkcorrect = predict.glm(PP_CaliMissGLkcorrect, data.frame(CA = 5))- predict.glm(PP_CaliMissGLkcorrect, data.frame(CA = 0)),
               ATE_CaliMissAllGLkcorrect = predict.glm(PP_CaliMissAllGLkcorrect, data.frame(CA = 5))- predict.glm(PP_CaliMissAllGLkcorrect, data.frame(CA = 0)),
               ATE_CaliAllMiss = predict.glm(PP_CaliAllMiss, data.frame(CA = 5))- predict.glm(PP_CaliAllMiss, data.frame(CA = 0)),
               ATE_CaliGLkcorrectbasic = predict.glm(PP_CaliGLkcorrectbasic, data.frame(CA = 5))- predict.glm(PP_CaliGLkcorrectbasic, data.frame(CA = 0)),
               ATE_CaliAllCorrectbasic = predict.glm(PP_CaliAllCorrectbasic, data.frame(CA = 5))- predict.glm(PP_CaliAllCorrectbasic, data.frame(CA = 0)),
               ATE_CaliGLkmissbasic = predict.glm(PP_CaliGLkmissbasic, data.frame(CA = 5))- predict.glm(PP_CaliGLkmissbasic, data.frame(CA = 0)),
               ATE_CaliAllGLkmissbasic = predict.glm(PP_CaliAllGLkmissbasic, data.frame(CA = 5))- predict.glm(PP_CaliAllGLkmissbasic, data.frame(CA = 0)),
               ATE_CaliMissGLkcorrectbasic = predict.glm(PP_CaliMissGLkcorrectbasic, data.frame(CA = 5))- predict.glm(PP_CaliMissGLkcorrectbasic, data.frame(CA = 0)),
               ATE_CaliMissAllGLkcorrectbasic = predict.glm(PP_CaliMissAllGLkcorrectbasic, data.frame(CA = 5))- predict.glm(PP_CaliMissAllGLkcorrectbasic, data.frame(CA = 0)),
               ATE_CaliAllMissbasic = predict.glm(PP_CaliAllMissbasic, data.frame(CA = 5))- predict.glm(PP_CaliAllMissbasic, data.frame(CA = 0)))
      
               
      result$predict_estimates <- ATE
     
      return(result)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  save(oper, file = paste("Simulation results/result_simu_linear_outcome_cali_test_",as.character(l),".rda", sep = ""))
}

scenarios <- tidyr::crossing(size, treat,conf)
iters <- 3

bias_ate <- array(,dim = c(19,5))
bias_hr <- array(, dim = c(19,2,4))
sd_ate <- array(,dim = c(19,5))
sd_hr <- array(, dim = c(19,2,4))
mae_ate <- array(,dim = c(19,5))
mae_hr <- array(, dim = c(19,2,4))
medae_ate <- array(,dim = c(19,5))
medae_hr <- array(, dim = c(19,2,4))
rootmse_ate <- array(,dim = c(19,5))
rootmse_hr <- array(, dim = c(19,2,4))


for (l in 1:4){
  load(paste0("Simulation results/result_simu_linear_outcome_cali_test_",as.character(l),".rda"))
  simu.t <- as.data.frame(1:iters)
  ate_all <- as.data.frame(list.rbind(oper[3,]))
  ate_all$simu <- simu.t[,1]
  ate_all <- ate_all[,-1]
  
  bias_ate[,1] <- colnames(ate_all[,-20])
  sd_ate[,1] <- colnames(ate_all[,-20])
  mae_ate[,1] <- colnames(ate_all[,-20])
  medae_ate[,1] <- colnames(ate_all[,-20])
  rootmse_ate[,1] <- colnames(ate_all[,-20])
  
  simu.scenario <- as.data.frame(tidyr::crossing(1:iters, 1:19))
  hr_all <- as.data.frame(list.rbind(oper[2,])) %>% 
    dplyr::mutate(simu = simu.scenario[,1], scenario =  simu.scenario[,2])
  
  bias_ate[,i+1] <- colMeans(ate_all[,-20]) - (-15)
  sd_ate[,i+1] <- colSds(as.matrix(ate_all[,-20]))
  mae_ate[,i+1] <- colMeans(abs(ate_all[,-20] - (-15)))
  medae_ate[,i+1] <- colMedians(as.matrix(abs(ate_all[,-20] - (-15))))
  rootmse_ate[,i+1] <- sqrt(as.numeric(bias_ate[,i+1])^2 +as.numeric(sd_ate[,i+1])^2)

  
}

colnames(bias_ate) <- c('Method', 'N = 200', 'N = 500', 'N = 1000', 'N = 5000')
colnames(sd_ate) <- c('Method', 'N = 200', 'N = 500', 'N = 1000', 'N = 5000')
colnames(mae_ate) <- c('Method', 'N = 200', 'N = 500', 'N = 1000', 'N = 5000')
colnames(medae_ate) <- c('Method', 'N = 200', 'N = 500', 'N = 1000', 'N = 5000')
colnames(rootmse_ate) <- c('Method', 'N = 200', 'N = 500', 'N = 1000', 'N = 5000')

bias_hr_table <- bias_hr_all[,,1]
sd_hr_table <- sd_hr_all[,,1]
rootmse_hr_table <- rootmse_hr_all[,,1]
mae_hr_table <- mae_hr_all[,,1]
medae_hr_table <- medae_hr_all[,,1]
for (i in 2:27){
  bias_hr_table <- rbind(bias_hr_table,bias_hr_all[,,i])
  sd_hr_table <- rbind(sd_hr_table,sd_hr_all[,,i])
  rootmse_hr_table <- rbind(rootmse_hr_table,rootmse_hr_all[,,i])
  mae_hr_table <- rbind(mae_hr_table, mae_hr_all[,,i])
  medae_hr_table <- rbind(medae_hr_table, medae_hr_all[,,i])
}

colnames(bias_hr_table) <- c('Bias 1', 'Bias 2', 'Bias 3')
colnames(sd_hr_table) <- c('SD 1', 'SD 2', 'SD 3')
colnames(rootmse_hr_table) <- c('rootMSE 1', 'rootMSE 2', 'rootMSE 3')
colnames(mae_hr_table) <- c('MAE 1', 'MAE 2', 'MAE 3')
colnames(medae_hr_table) <- c('MedAE 1', 'MedAE 2', 'MedAE 3')
size <- c(200,500,1000)
treat <- c(0.25,0.5,0.75)
conf <- c(1,1.5,2)
var <- c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
         'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
scenarios <- tidyr::crossing(size, treat,conf,var) %>% 
  dplyr::arrange(size, treat, conf,
                 match(var, c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
                              'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
                 ))

na_list <- c(6)
for (r in 2:54){
  na_list <- cbind(na_list, r*5 + r)
}
hr_table <- cbind(scenarios,bias_hr_table, sd_hr_table,rootmse_hr_table, mae_hr_table,medae_hr_table) %>% 
  dplyr::select(size, treat, conf,var, 'Bias 1', 'SD 1', 'rootMSE 1', 'MAE 1','MedAE 1','Bias 2', 'SD 2', 'rootMSE 2','MAE 2','MedAE 2','Bias 3', 'SD 3', 'rootMSE 3', 'MAE 3','MedAE 3')%>% 
  insertRows(r = na_list, new = NA)


print(xtable(hr_table[hr_table$size == 200,],digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')

print(xtable(hr_table[hr_table$size == 500,],digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')

print(xtable(hr_table[hr_table$size == 1000,],digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')

write_csv(hr_table, file = 'simulation_results_hr_table.csv')

bias_ate_table <- bias_ate_all[,,1]
sd_ate_table <- sd_ate_all[,,1]
rootmse_ate_table <- rootmse_ate_all[,,1]
mae_ate_table <- mae_ate_all[,,1]
medae_ate_table <- medae_ate_all[,,1]
for (i in 2:27){
  bias_ate_table <- rbind(bias_ate_table,bias_ate_all[,,i])
  sd_ate_table <- rbind(sd_ate_table,sd_ate_all[,,i])
  rootmse_ate_table <- rbind(rootmse_ate_table,rootmse_ate_all[,,i])
  mae_ate_table <- rbind(mae_ate_table, mae_ate_all[,,i])
  medae_ate_table <- rbind(medae_ate_table, medae_ate_all[,,i])
}

colnames(bias_ate_table) <- c('Bias 0', 'Bias 1','Bias 2', 'Bias 3','Bias 4')
colnames(sd_ate_table) <- c('SD 0', 'SD 1','SD 2', 'SD 3','SD 4')
colnames(rootmse_ate_table) <- c('rootMSE 0', 'rootMSE 1','rootMSE 2', 'rootMSE 3','rootMSE 4')
colnames(mae_ate_table) <- c('MAE 0','MAE 1', 'MAE 2', 'MAE 3','MAE 4')
colnames(medae_ate_table) <- c('MedAE 0','MedAE 1', 'MedAE 2', 'MedAE 3','MedAE 4')
size <- c(200,500,1000)
treat <- c(0.25,0.50,0.75)
conf <- c(1,1.5,2)
var <- c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
         'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
scenarios <- tidyr::crossing(size, treat,conf,var)%>% 
  dplyr::arrange(size, treat, conf,
                 match(var, c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
                              'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')))


ate_table <- cbind(scenarios,bias_ate_table, sd_ate_table,rootmse_ate_table,mae_ate_table,medae_ate_table) %>% 
  dplyr::select(size, treat, conf,var,
                'Bias 0', 'SD 0', 'rootMSE 0', 'MAE 0','MedAE 0',
                'Bias 1', 'SD 1', 'rootMSE 1', 'MAE 1','MedAE 1',
                'Bias 2', 'SD 2', 'rootMSE 2','MAE 2','MedAE 2',
                'Bias 3', 'SD 3', 'rootMSE 3', 'MAE 3','MedAE 3',
                'Bias 4', 'SD 4', 'rootMSE 4', 'MAE 4','MedAE 4')%>% 
  insertRows(r = na_list, new = NA)

write_csv(ate_table, file = 'simulation_results_ate_table.csv')
print(xtable(ate_table[ate_table$size == 200,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')
print(xtable(ate_table[ate_table$size == 500,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')
print(xtable(ate_table[ate_table$size == 1000,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')


bias_mean1_table <- bias_mean1_all[,,1]
sd_mean1_table <- sd_mean1_all[,,1]
rootmse_mean1_table <- rootmse_mean1_all[,,1]
mae_mean1_table <- mae_mean1_all[,,1]
medae_mean1_table <- medae_mean1_all[,,1]
for (i in 2:27){
  bias_mean1_table <- rbind(bias_mean1_table,bias_mean1_all[,,i])
  sd_mean1_table <- rbind(sd_mean1_table,sd_mean1_all[,,i])
  rootmse_mean1_table <- rbind(rootmse_mean1_table,rootmse_mean1_all[,,i])
  mae_mean1_table <- rbind(mae_mean1_table, mae_mean1_all[,,i])
  medae_mean1_table <- rbind(medae_mean1_table, medae_mean1_all[,,i])
}

colnames(bias_mean1_table) <- c('Bias 0', 'Bias 1','Bias 2', 'Bias 3','Bias 4')
colnames(sd_mean1_table) <- c('SD 0', 'SD 1','SD 2', 'SD 3','SD 4')
colnames(rootmse_mean1_table) <- c('rootMSE 0', 'rootMSE 1','rootMSE 2', 'rootMSE 3','rootMSE 4')
colnames(mae_mean1_table) <- c('MAE 0','MAE 1', 'MAE 2', 'MAE 3','MAE 4')
colnames(medae_mean1_table) <- c('MedAE 0','MedAE 1', 'MedAE 2', 'MedAE 3','MedAE 4')
size <- c(200,500,1000)
treat <- c(0.25,0.50,0.75)
conf <- c(1,1.5,2)
var <- c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
         'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
scenarios <- tidyr::crossing(size, treat,conf,var)%>% 
  dplyr::arrange(size, treat, conf,
                 match(var, c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
                              'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')))


mean1_table <- cbind(scenarios,bias_mean1_table, sd_mean1_table,rootmse_mean1_table,mae_mean1_table,medae_mean1_table) %>% 
  dplyr::select(size, treat, conf,var,
                'Bias 0', 'SD 0', 'rootMSE 0', 'MAE 0','MedAE 0',
                'Bias 1', 'SD 1', 'rootMSE 1', 'MAE 1','MedAE 1',
                'Bias 2', 'SD 2', 'rootMSE 2','MAE 2','MedAE 2',
                'Bias 3', 'SD 3', 'rootMSE 3', 'MAE 3','MedAE 3',
                'Bias 4', 'SD 4', 'rootMSE 4', 'MAE 4','MedAE 4')%>% 
  insertRows(r = na_list, new = NA)

write_csv(mean1_table, file = 'simulation_results_mean1_table.csv')
print(xtable(mean1_table[mean1_table$size == 200,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')
print(xtable(mean1_table[mean1_table$size == 500,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')
print(xtable(mean1_table[mean1_table$size == 1000,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')

