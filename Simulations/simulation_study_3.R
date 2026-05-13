#!/usr/bin R
###############################
# This script runs simulations with Q functions on logit form,
# and then predicted values are scaled before targeting step.
# Number of visits: 10
#
# For k = T,..., 0:
#   Fit logit model of Q_{k+1}^d regressed on \bar A_k, \bar L_k
#   Obtain predicted values Q_{k}^d
#   Use logistic targeting step on scaled Q_{k}^d
###############################
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/utils/dgm_10_visits.R")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/utils/calibration_func_trials.R")
library(MASS)
library(foreach)
library(doParallel)
library(doRNG)
library(nleqslv)
library(data.table)
library(matrixStats)
library(xtable)

set.seed(160625)
seeds <- floor(runif(1000)*10^8)

simulation_code <- function(iters, transformed = FALSE, sample_size = 300,seeds,conf = 0.2, 
                            treat_prev_0 = 0, 
                            treat_prev_d1 = c(1.85,1.65,1.45,1.25,1.05,0.85,0.65,0.45,0.25), 
                            treat_prev_d0 = c(-2.15, -2.15, -2.15, -2.15, -2.15, -2.15, -2.15, -2.15, -2.15)){
  time <- proc.time()
  if(transformed){
    print("Transformed covariates")
  } else {print("Correct covariates")}
  simulation <- foreach(i = 1:iters, .combine=cbind) %dopar% {
    set.seed(seeds[i])
    suppressMessages(suppressWarnings({
      simdata<-DATA_GEN_TEN(ns = sample_size, conf = conf,treat_prev_0 = treat_prev_0, 
                            treat_prev_d1 = treat_prev_d1, 
                            treat_prev_d0 = treat_prev_d0)
      if(transformed){
        simdata$X1 <- simdata$TX1
        simdata$X2 <- simdata$TX2
        simdata$X3 <- simdata$TX3
        simdata$X4 <- simdata$TX4
      }
      
      
      treat_model_A_0 <- glm(A~X1 + X2 + X3 + X4-1, data = simdata[simdata$t == 0,], family = 'binomial')
      treat_model_A_1 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 1,], family = 'binomial')
      treat_model_A_2 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 2,], family = 'binomial')
      treat_model_A_3 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 3,], family = 'binomial')
      treat_model_A_4 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 4,], family = 'binomial')
      treat_model_A_5 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 5,], family = 'binomial')
      treat_model_A_6 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 6,], family = 'binomial')
      treat_model_A_7 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 7,], family = 'binomial')
      treat_model_A_8 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 8,], family = 'binomial')
      treat_model_A_9 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 9,], family = 'binomial')
      
      simdata$ps <- 1.0
      simdata[simdata$t == 0& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 1& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 2& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 3& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_3, newdata = simdata[simdata$t == 3& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 4& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_4, newdata = simdata[simdata$t == 4& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 5& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_5, newdata = simdata[simdata$t == 5& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 6& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_6, newdata = simdata[simdata$t == 6& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 7& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_7, newdata = simdata[simdata$t == 7& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 8& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_8, newdata = simdata[simdata$t == 8& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 9& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_9, newdata = simdata[simdata$t == 9& simdata$A == 1,], type = 'response'))
      
      simdata[simdata$t == 0& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 1& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 2& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 3& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_3, newdata = simdata[simdata$t == 3& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 4& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_4, newdata = simdata[simdata$t == 4& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 5& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_5, newdata = simdata[simdata$t == 5& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 6& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_6, newdata = simdata[simdata$t == 6& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 7& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_7, newdata = simdata[simdata$t == 7& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 8& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_8, newdata = simdata[simdata$t == 8& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 9& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_9, newdata = simdata[simdata$t == 9& simdata$A == 0,], type = 'response'))
      
      simdata$weight <- ave(simdata$ps, simdata$ID, FUN = function(X) 1/cumprod(X))
      simdata$tall <- simdata$t
      simdata$sub <- simdata$ID
      simdata$RA <- 1
      simdata[simdata$t == 0 & !(simdata$CA == 1),]$RA <- 0
      simdata[simdata$t == 1 & !(simdata$CA == 2),]$RA <- 0
      simdata[simdata$t == 2 & !(simdata$CA == 3),]$RA <- 0
      simdata[simdata$t == 3 & !(simdata$CA == 4),]$RA <- 0
      simdata[simdata$t == 4 & !(simdata$CA == 5),]$RA <- 0
      simdata[simdata$t == 5 & !(simdata$CA == 6),]$RA <- 0
      simdata[simdata$t == 6 & !(simdata$CA == 7),]$RA <- 0
      simdata[simdata$t == 7 & !(simdata$CA == 8),]$RA <- 0
      simdata[simdata$t == 8 & !(simdata$CA == 9),]$RA <- 0
      simdata[simdata$t == 9 & !(simdata$CA == 10),]$RA <- 0
      simdata$tX1 <- simdata$t*simdata$X1
      simdata$tX2 <- simdata$t*simdata$X2
      simdata$tX3 <- simdata$t*simdata$X3
      simdata$tX4 <- simdata$t*simdata$X4
      
      calibrate_always_treated <- calibration_by_time_from_baseline(simdata, var = c("X1", "X2", "X3", "X4"), weights_var = 'weight')
      calibrate_always_treated_aggr <- aggregated_calibration_from_baseline(simdata, var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), weights_var = 'weight')
      
      simdata$Cweights <- calibrate_always_treated$data$Cweights
      simdata$Cweights_aggr <- calibrate_always_treated_aggr$data$Cweights
      
      simdata$RA <- 1
      simdata[simdata$t == 0 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 1 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 2 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 3 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 4 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 5 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 6 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 7 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 8 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 9 & !(simdata$CA == 0),]$RA <- 0
      
      calibrate_never_treated <- calibration_by_time_from_baseline(simdata, var = c("X1", "X2", "X3", "X4"), weights_var = 'Cweights')
      calibrate_never_treated_aggr <- aggregated_calibration_from_baseline(simdata, var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), weights_var = 'Cweights_aggr')
      
      simdata$Cweights <- calibrate_never_treated$data$Cweights
      simdata$Cweights_aggr <- calibrate_never_treated_aggr$data$Cweights
      
      simdata$weights <- simdata$weight
      
      
      #con4<-xtabs(~t + switch + A, data=simdata)
      #ftable(con4)
      
      # library(cobalt)
      #
      # balance_list <- lapply(unique(simdata$t), function(time) {
      #   dat <- subset(simdata, t == time)
      #   bal <- bal.tab(A ~ X1 + X2,
      #                  data = dat,
      #                  estimand = "ATE")
      #   bal
      # })
      #
      # names(balance_list) <- unique(simdata$t)
      # balance_list[[3]]
      
      
      #con4<-xtabs(~followup_time + assigned_treatment, data=switch_data)
      #ftable(con4)
      
      ########## Manual LTMLE MSM ############
      wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, value.var = c("A", "X1", "X2","X3", "X4", "CA", "Y", "weights", "Cweights", "Cweights_aggr"))
      
      a = min(simdata$Y,na.rm = TRUE); b = max(simdata$Y, na.rm = TRUE)
      
      wideSimdata$Y_9_scaled = (wideSimdata$Y_9-a)/(b-a)
      wideSimdata$Y_8_scaled = (wideSimdata$Y_8-a)/(b-a)
      wideSimdata$Y_7_scaled = (wideSimdata$Y_7-a)/(b-a)
      wideSimdata$Y_6_scaled = (wideSimdata$Y_6-a)/(b-a)
      wideSimdata$Y_5_scaled = (wideSimdata$Y_5-a)/(b-a)
      wideSimdata$Y_4_scaled = (wideSimdata$Y_4-a)/(b-a)
      wideSimdata$Y_3_scaled = (wideSimdata$Y_3-a)/(b-a)
      wideSimdata$Y_2_scaled = (wideSimdata$Y_2-a)/(b-a)
      wideSimdata$Y_1_scaled = (wideSimdata$Y_1-a)/(b-a)
      wideSimdata$Y_0_scaled = (wideSimdata$Y_0-a)/(b-a)
      
      
      logitQforms <- vector("list", 10 * 10)
      dim(logitQforms) <- c(10, 10)
      
      logitQforms_cali <- vector("list", 10 * 10)
      dim(logitQforms_cali) <- c(10, 10)
      
      logitQforms_cali_aggr <- vector("list", 10 * 10)
      dim(logitQforms_cali_aggr) <- c(10, 10)
      
      Qstars <- vector("list", 10 * 10)
      dim(Qstars) <- c(10, 10)
      
      Qstars_cali <- vector("list", 10 * 10)
      dim(Qstars_cali) <- c(10, 10)
      
      Qstars_cali_aggr <- vector("list", 10 * 10)
      dim(Qstars_cali_aggr) <- c(10, 10)
      
      ########### t = 9 ##########
      
      for (k in 0:9){
        for (j in k:0){
          regimen_ind <- c(as.numeric(wideSimdata[[paste0('CA_',j)]] ==j+1), as.numeric(wideSimdata[[paste0('CA_',j)]] ==0))
          weight <- rep(wideSimdata[[paste0('weights_',j)]],2)
          Cweight <- rep(wideSimdata[[paste0('Cweights_',j)]],2)
          Cweight_aggr <- rep(wideSimdata[[paste0('Cweights_aggr_',j)]],2)
          
          if (k == j){
            if(j > 0){ 
              formula_string <- paste0('Y_',j,'_scaled ~ A_',j,' + A_',j-1, '+ X1_',j,' + X2_',j,' + X3_',j,' +X4_',j,' + X1_',j-1, ' + X2_',j-1, ' + X3_',j-1, ' +X4_',j-1)
              predict_data <- data.frame(c(rep(1,sample_size), rep(0,sample_size)),c(rep(1, sample_size), rep(0, sample_size)))
              colnames(predict_data) <- paste0("A_", j:(j-1))
              for (l in j:(j-1)){
                predict_data[, paste0('X1_',l)] <- rep(wideSimdata[[paste0('X1_',l)]], 2)
                predict_data[, paste0('X2_',l)] <- rep(wideSimdata[[paste0('X2_',l)]], 2)
                predict_data[, paste0('X3_',l)] <- rep(wideSimdata[[paste0('X3_',l)]], 2)
                predict_data[, paste0('X4_',l)] <- rep(wideSimdata[[paste0('X4_',l)]], 2)
              }
              
            } else {
              formula_string <- paste0('Y_',j,'_scaled ~ A_',j, '+ X1_',j,' + X2_',j,' + X3_',j,' +X4_',j)
              predict_data <- data.frame(c(rep(1,sample_size), rep(0,sample_size)))
              colnames(predict_data) <- paste0("A_", j)
              predict_data[, paste0('X1_',j)] <- rep(wideSimdata[[paste0('X1_',j)]], 2)
              predict_data[, paste0('X2_',j)] <- rep(wideSimdata[[paste0('X2_',j)]], 2)
              predict_data[, paste0('X3_',j)] <- rep(wideSimdata[[paste0('X3_',j)]], 2)
              predict_data[, paste0('X4_',j)] <- rep(wideSimdata[[paste0('X4_',j)]], 2)
              
            }
            Qform <- glm(data = wideSimdata, formula = as.formula(formula_string), family = 'quasibinomial')
            
            logitQforms[[(k+1),(j+1)]] <- predict.glm(Qform, newdata = predict_data, type = 'link')
            
            
            
            update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                      off = logitQforms[[(k+1),(j+1)]],
                                      intercept = rep(1,sample_size*2), 
                                      cumA = c(rep(k+1,sample_size), rep(0,sample_size)))
            update_data[[paste0('Y_',j)]] <-rep(wideSimdata[[paste0('Y_',j,'_scaled')]], 2)
            
            formula_string <- paste0('Y_',j,' ~ intercept + cumA + offset(off) - 1')
            
            Qstar_fit <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(weight[regimen_ind == 1]))
            
            Qstars[[(k+1),(j+1)]] <- predict.glm(Qstar_fit, newdata = update_data, type = 'response')
            
            Qstar_fit_cali <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                  data = update_data[regimen_ind == 1,],
                                  weights = as.vector(Cweight[regimen_ind == 1]))
            
            Qstars_cali[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali, newdata = update_data, type = 'response')
            
            Qstar_fit_cali_aggr <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                       data = update_data[regimen_ind == 1,],
                                       weights = as.vector(Cweight_aggr[regimen_ind == 1]))
            
            Qstars_cali_aggr[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_aggr, newdata = update_data, type = 'response')
            
          } else{ 
            fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                                       Qk_j1star = Qstars[[(k+1),(j+2)]], 
                                       Qk_j1star_cali = Qstars_cali[[(k+1),(j+2)]],
                                       Qk_j1star_cali_aggr = Qstars_cali_aggr[[(k+1),(j+2)]],
                                       CA_j = rep(wideSimdata[[paste0('CA_',j)]], 2),
                                       X1_j = rep(wideSimdata[[paste0('X1_',j)]], 2), 
                                       X2_j = rep(wideSimdata[[paste0('X2_',j)]], 2),
                                       X3_j = rep(wideSimdata[[paste0('X3_',j)]], 2), 
                                       X4_j = rep(wideSimdata[[paste0('X4_',j)]], 2))
            
            Q_predict_data_d1 <- fitting_data[1:sample_size,]
            Q_predict_data_d1$CA_j <- rep(j+1,sample_size)
            
            Q_predict_data_d0 <- fitting_data[1:sample_size,]
            Q_predict_data_d0$CA_j <- rep(0,sample_size)
            if (j == (k-1)){
              formula_string <- paste0('Qk_j1star ~ CA_j + X1_j + X2_j + X3_j + X4_j')
              formula_string_cali <- paste0('Qk_j1star_cali ~ CA_j + X1_j + X2_j + X3_j + X4_j')
              formula_string_cali_aggr <- paste0('Qk_j1star_cali_aggr ~ CA_j + X1_j + X2_j + X3_j + X4_j')
            } else if (j == (k-2)){
              formula_string <- paste0('Qk_j1star ~ CA_j')
              formula_string_cali <- paste0('Qk_j1star_cali ~ CA_j ')
              formula_string_cali_aggr <- paste0('Qk_j1star_cali_aggr ~ CA_j ')
            } else if (j < (k-2)){
              formula_string <- paste0('Qk_j1star ~ 1')
              formula_string_cali <- paste0('Qk_j1star_cali ~ 1 ')
              formula_string_cali_aggr <- paste0('Qk_j1star_cali_aggr ~ 1')
            }
            Qform_d1 <- glm(data = fitting_data[1:sample_size,], formula = as.formula(formula_string),
                            family = 'quasibinomial')
            Qform_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = as.formula(formula_string),
                            family = 'quasibinomial')
            
            logitQforms[[(k+1), (j+1)]] <- as.matrix(c(predict.glm(Qform_d1, 
                                                                           newdata = Q_predict_data_d1, type = 'link'),
                                                       predict.glm(Qform_d0, 
                                                                           newdata = Q_predict_data_d0, type = 'link')
            ))
            
            Qform_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = as.formula(formula_string_cali),
                                 family = 'quasibinomial')
            Qform_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = as.formula(formula_string_cali),
                                 family = 'quasibinomial')
            
            logitQforms_cali[[(k+1), (j+1)]] <- as.matrix(c(predict.glm(Qform_d1_cali, 
                                                                                newdata = Q_predict_data_d1, type = 'link'),
                                                            predict.glm(Qform_d0_cali, 
                                                                                newdata = Q_predict_data_d0, type = 'link')
            ))
            
            Qform_d1_cali_aggr <- glm(data = fitting_data[1:sample_size,], formula = as.formula(formula_string_cali_aggr),
                                      family = 'quasibinomial')
            Qform_d0_cali_aggr <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = as.formula(formula_string_cali_aggr),
                                      family = 'quasibinomial')
            
            logitQforms_cali_aggr[[(k+1), (j+1)]] <- as.matrix(c(predict.glm(Qform_d1_cali_aggr, 
                                                                                     newdata = Q_predict_data_d1, type = 'link'),
                                                                 predict.glm(Qform_d0_cali_aggr, 
                                                                                     newdata = Q_predict_data_d0, type = 'link')
            ))
            
            
            update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                      Qk_j1star = Qstars[[(k+1),(j+2)]], 
                                      Qk_j1star_cali = Qstars_cali[[(k+1),(j+2)]],
                                      Qk_j1star_cali_aggr = Qstars_cali_aggr[[(k+1),(j+2)]],
                                      off = logitQforms[[(k+1),(j+1)]],
                                      off_cali = logitQforms_cali[[(k+1),(j+1)]],
                                      off_cali_aggr = logitQforms_cali_aggr[[(k+1),(j+1)]],
                                      intercept = rep(1,sample_size*2), 
                                      cumA = c(rep(k+1,sample_size), rep(0,sample_size)))
            
            Qk_j_star_fit <- glm(Qk_j1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                                 data = update_data[regimen_ind == 1,],
                                 weights = as.vector(weight[regimen_ind == 1]))
            
            Qstars[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit, newdata = update_data, type = 'response')
            
            Qk_j_star_fit_cali <- glm(Qk_j1star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                                      data = update_data[regimen_ind == 1,],
                                      weights = as.vector(Cweight[regimen_ind == 1]))
            
            Qstars_cali[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali, newdata = update_data, type = 'response')
            
            Qk_j_star_fit_cali_aggr <- glm(Qk_j1star_cali_aggr ~ intercept + cumA + offset(off_cali_aggr) - 1, family = 'quasibinomial', 
                                           data = update_data[regimen_ind == 1,],
                                           weights = as.vector(Cweight_aggr[regimen_ind == 1]))
            
            Qstars_cali_aggr[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_aggr, newdata = update_data, type = 'response')
            
          }
        }
      }
      ################# Fit MSM ##########################
      
      msm_fitting_data <- data.frame(id = rep(1:sample_size,20), 
                                     t = c(rep(0,sample_size*2), 
                                           rep(1,sample_size*2), 
                                           rep(2,sample_size*2),
                                           rep(3,sample_size*2),
                                           rep(4,sample_size*2),
                                           rep(5,sample_size*2),
                                           rep(6,sample_size*2),
                                           rep(7,sample_size*2),
                                           rep(8,sample_size*2),
                                           rep(9,sample_size*2)),
                                     Y =  unlist(Qstars[, 1])*(b-a) + a, 
                                     Y_cali =  unlist(Qstars_cali[, 1])*(b-a) + a,
                                     Y_cali_aggr =  unlist(Qstars_cali_aggr[, 1])*(b-a) + a,
                                     cumA = c(rep(1,sample_size), rep(0,sample_size), 
                                              rep(2,sample_size), rep(0,sample_size), 
                                              rep(3,sample_size), rep(0,sample_size),
                                              rep(4,sample_size), rep(0,sample_size),
                                              rep(5,sample_size), rep(0,sample_size),
                                              rep(6,sample_size), rep(0,sample_size),
                                              rep(7,sample_size), rep(0,sample_size),
                                              rep(8,sample_size), rep(0,sample_size),
                                              rep(9,sample_size), rep(0,sample_size),
                                              rep(10,sample_size), rep(0,sample_size)))
      
      msm <- glm(Y ~ cumA, data = msm_fitting_data)
      msm_cali <- glm(Y_cali ~ cumA, data = msm_fitting_data)
      msm_cali_aggr <- glm(Y_cali_aggr ~ cumA, data = msm_fitting_data)
      
      
      ##################### IPW - MSM ##################################
      
      simdata$RA <- 1
      simdata[simdata$t == 0 & !(simdata$CA == 1 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 1 & !(simdata$CA == 2 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 2 & !(simdata$CA == 3 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 3 & !(simdata$CA == 4 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 4 & !(simdata$CA == 5 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 5 & !(simdata$CA == 6 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 6 & !(simdata$CA == 7 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 7 & !(simdata$CA == 8 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 8 & !(simdata$CA == 9 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 9 & !(simdata$CA == 10 | simdata$CA == 0),]$RA <- 0
      
      switch_data <- simdata[simdata$RA == 1,]
      
      switch_data$X1 <- ave(switch_data$X1, switch_data$ID, FUN = first)
      switch_data$X2 <- ave(switch_data$X2, switch_data$ID, FUN = first)
      switch_data$X3 <- ave(switch_data$X3, switch_data$ID, FUN = first)
      switch_data$X4 <- ave(switch_data$X4, switch_data$ID, FUN = first)
      
      
      PP_strat <- glm(Y ~ CA, data = switch_data, weights = weights, family = 'gaussian')
      PP_cali <- glm(Y ~ CA, data = switch_data, weights = Cweights, family = 'gaussian')
      PP_cali_aggr <- glm(Y ~ CA, data = switch_data, weights = Cweights_aggr, family = 'gaussian')
      
    }))
    
    
    c(PP_strat$coefficients[1],
      PP_cali$coefficients[1], 
      PP_cali_aggr$coefficients[1], 
      msm$coefficients[1],
      msm_cali$coefficients[1],
      msm_cali_aggr$coefficients[1],
      PP_strat$coefficients[2],
      PP_cali$coefficients[2],
      PP_cali_aggr$coefficients[2],
      msm$coefficients[2],
      msm_cali$coefficients[2],
      msm_cali_aggr$coefficients[2]
    )
  }
  cat(paste('Transformed? ', transformed, ' \n'))
  cat(paste('Sample size = ', sample_size, ' \n'))
  cat(paste('Confounding = ', conf, ' \n'))
  
  cat(paste('Time taken for iters = ', iters, ' \n'))
  print(proc.time() - time)
  
  cat('\n')
  rownames(simulation) <- c('MLE',
                            'CMLE',
                            'Aggr. CMLE',
                            'MLE LTMLE',
                            'CMLE LTMLE',
                            'Aggr. CMLE LTMLE',
                            'MLE',
                            'CMLE',
                            'Aggr. CMLE',
                            'MLE LTMLE',
                            'CMLE LTMLE',
                            'Aggr. CMLE LTMLE'
  )
  
  results <- cbind(rowMeans(simulation, na.rm = TRUE)[1:6]-200, rowMeans(simulation, na.rm = TRUE)[7:12]-10,
                   rowSds(simulation, na.rm = TRUE)[1:6], rowSds(simulation, na.rm = TRUE)[7:12],
                   sqrt((rowMeans(simulation, na.rm = TRUE)[1:6]-200)^2+ rowSds(simulation, na.rm = TRUE)[1:6]^2),
                   sqrt((rowMeans(simulation, na.rm = TRUE)[7:12]-10)^2+ rowSds(simulation, na.rm = TRUE)[7:12]^2))
  
  always_treat_Y <- cbind(rowMeans(simulation[1:6,] + 3*simulation[7:12,], na.rm = TRUE)-230,
                          rowSds(simulation[1:6,] + 3*simulation[7:12,], na.rm = TRUE),
                          sqrt((rowMeans(simulation[1:6,] + 3*simulation[7:12,], na.rm = TRUE)-230)^2 + rowSds(simulation[1:6,] + 3*simulation[7:12,], na.rm = TRUE)^2))
  
  return(list(MSM = results,
              EYT = always_treat_Y))
}


iters = 1000
registerDoParallel(cores = 10)

sink("10_visits_simu_logitmodels.txt")

low_300_correct <- simulation_code(iters = iters, sample_size = 300, conf = 0.2, seeds = seeds)
low_300_trans <- simulation_code(iters = iters, transformed = TRUE, sample_size = 300, conf = 0.2, seeds = seeds)
low_500_correct <- simulation_code(iters = iters, sample_size = 500, conf = 0.2, seeds = seeds)
low_500_trans <- simulation_code(iters = iters, transformed = TRUE, sample_size = 500, conf = 0.2, seeds = seeds)
low_1000_correct <- simulation_code(iters = iters, sample_size = 1000, conf = 0.2, seeds = seeds)
low_1000_trans <- simulation_code(iters = iters, transformed = TRUE, sample_size = 1000, conf = 0.2, seeds = seeds)
low_2500_correct <- simulation_code(iters = iters, sample_size = 2500, conf = 0.2, seeds = seeds)
low_2500_trans <- simulation_code(iters = iters, transformed = TRUE, sample_size = 2500, conf = 0.2, seeds = seeds)

high_300_correct <- simulation_code(iters = iters, sample_size = 300,  seeds = seeds,conf = 2.5, 
                                    treat_prev_d1 = c(2,-0.5,-3,-5.5,-8,-10.5,-13,-15.5,-18), 
                                    treat_prev_d0 = c(-4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55))
high_300_trans <- simulation_code(iters = iters, transformed = TRUE, sample_size = 300,  seeds = seeds,conf = 2.5, 
                                  treat_prev_d1 = c(2,-0.5,-3,-5.5,-8,-10.5,-13,-15.5,-18), 
                                  treat_prev_d0 = c(-4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55))
high_500_correct <- simulation_code(iters = iters, sample_size = 500,  seeds = seeds,conf = 2.5,
                                    treat_prev_d1 = c(2,-0.5,-3,-5.5,-8,-10.5,-13,-15.5,-18),
                                    treat_prev_d0 = c(-4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55))
high_500_trans <- simulation_code(iters = iters, transformed = TRUE, sample_size = 500,  seeds = seeds,conf = 2.5,
                                  treat_prev_d1 = c(2,-0.5,-3,-5.5,-8,-10.5,-13,-15.5,-18),
                                  treat_prev_d0 = c(-4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55))
high_1000_correct <- simulation_code(iters = iters, sample_size = 1000,  seeds = seeds,conf = 2.5,
                                     treat_prev_d1 = c(2,-0.5,-3,-5.5,-8,-10.5,-13,-15.5,-18),
                                     treat_prev_d0 = c(-4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55))
high_1000_trans <- simulation_code(iters = iters, transformed = TRUE, sample_size = 1000,  seeds = seeds,conf = 2.5,
                                   treat_prev_d1 = c(2,-0.5,-3,-5.5,-8,-10.5,-13,-15.5,-18),
                                   treat_prev_d0 = c(-4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55))
high_2500_correct <- simulation_code(iters = iters, sample_size = 2500,  seeds = seeds,conf = 2.5,
                                     treat_prev_d1 = c(2,-0.5,-3,-5.5,-8,-10.5,-13,-15.5,-18),
                                     treat_prev_d0 = c(-4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55))
high_2500_trans <- simulation_code(iters = iters, transformed = TRUE, sample_size = 2500,  seeds = seeds,conf = 2.5,
                                   treat_prev_d1 = c(2,-0.5,-3,-5.5,-8,-10.5,-13,-15.5,-18),
                                   treat_prev_d0 = c(-4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55, -4.55))

cat(paste('Table 4 results \n'))

print(xtable(cbind(low_300_correct$MSM,
                   low_300_trans$MSM),
             type = "latex"))
print(xtable(cbind(low_500_correct$MSM,
                   low_500_trans$MSM),
             type = "latex"))
print(xtable(cbind(low_1000_correct$MSM,
                   low_1000_trans$MSM),
             type = "latex"))
print(xtable(cbind(low_2500_correct$MSM,
                   low_2500_trans$MSM),
             type = "latex"))
cat(paste('%---------------------------------------------------------------------------------- \n'))
cat(paste('Table 5 results \n'))
print(xtable(cbind(high_300_correct$MSM,
                   high_300_trans$MSM),
             type = "latex"))
print(xtable(cbind(high_500_correct$MSM,
                   high_500_trans$MSM),
             type = "latex"))
print(xtable(cbind(high_1000_correct$MSM,
                   high_1000_trans$MSM),
             type = "latex"))
print(xtable(cbind(high_2500_correct$MSM,
                   high_2500_trans$MSM),
             type = "latex"))

cat(paste('%---------------------------------------------------------------------------------- \n'))

cat(paste('Always-treated results \n'))
print(xtable(cbind(low_300_correct$EYT,
                   low_300_trans$EYT),
             type = "latex"))
print(xtable(cbind(low_500_correct$EYT,
                   low_500_trans$EYT),
             type = "latex"))
print(xtable(cbind(low_1000_correct$EYT,
                   low_1000_trans$EYT),
             type = "latex"))
print(xtable(cbind(low_2500_correct$EYT,
                   low_2500_trans$EYT),
             type = "latex"))
cat(paste('%---------------------------------------------------------------------------------- \n'))

print(xtable(cbind(high_300_correct$EYT,
                   high_300_trans$EYT),
             type = "latex"))
print(xtable(cbind(high_500_correct$EYT,
                   high_500_trans$EYT),
             type = "latex"))
print(xtable(cbind(high_1000_correct$EYT,
                   high_1000_trans$EYT),
             type = "latex"))
print(xtable(cbind(high_2500_correct$EYT,
                   high_2500_trans$EYT),
             type = "latex"))

