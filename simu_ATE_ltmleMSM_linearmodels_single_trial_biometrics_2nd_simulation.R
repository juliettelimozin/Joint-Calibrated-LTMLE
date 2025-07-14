#!/usr/bin R
###############################
# This script runs simulations with Q functions on linear form,
# and then predicted values are scaled before targeting step.
# For k = T,..., 0:
#   Fit linear model of unscaled Q_{k+1}^d regressed on \bar A_k, \bar L_k
#   Obtain predicted values Q_{k}^d and scale them to (0,1) using minimum and maximum values of Y
#   Use logistic targeting step on scaled Q_{k}^d
###############################
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/dgm_2nd_simulation_biometrics.R")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/calibration_func_trials.R")
library(MASS)
library(foreach)
library(doParallel)
library(doRNG)
library(nleqslv)
library(ltmle)
library(data.table)
library(matrixStats)
library(xtable)

set.seed(160625)
seeds <- floor(runif(1000)*10^8)

simulation_code <- function(iters, transformed = FALSE, sample_size,seeds,conf = 0.2, 
                            treat_prev_0 = 0, 
                            treat_prev_d1_1 = 1, 
                            treat_prev_d0_1 = -1.25, 
                            treat_prev_d1_2 =0.8, 
                            treat_prev_d0_2 = -1.25){
  time <- proc.time()
  if(transformed){
    print("Transformed covariates")
  } else {print("Correct covariates")}
  simulation <- foreach(i = 1:iters, .combine=cbind) %dopar% {
    set.seed(seeds[i])
    suppressMessages(suppressWarnings({
      simdata<-DATA_GEN(ns = sample_size, conf = conf,treat_prev_0 = treat_prev_0, 
                        treat_prev_d1_1 = treat_prev_d1_1, 
                        treat_prev_d0_1 = treat_prev_d0_1, 
                        treat_prev_d1_2 =treat_prev_d1_2, 
                        treat_prev_d0_2 = treat_prev_d0_2)
      if(transformed){
        simdata$X1 <- simdata$TX1
        simdata$X2 <- simdata$TX2
        simdata$X3 <- simdata$TX3
        simdata$X4 <- simdata$TX4
      }
      
      treatment_model_pooled <- glm(A~Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t !=0,], family = 'binomial')
      
      treat_model_A_0 <- glm(A~X1 + X2 + X3 + X4-1, data = simdata[simdata$t == 0,], family = 'binomial')
      treat_model_A_1 <- glm(A~ Ap+X1 + X2 + X3 + X4, data = simdata[simdata$t == 1,], family = 'binomial')
      treat_model_A_2 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 2,], family = 'binomial')
      
      
      simdata$ps <- 1.0
      simdata[simdata$t == 0& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 1& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 1,], type = 'response'))
      simdata[simdata$t == 2& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 1,], type = 'response'))
      
      simdata[simdata$t == 0& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 1& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 0,], type = 'response'))
      simdata[simdata$t == 2& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 0,], type = 'response'))
      
      simdata$weight <- ave(simdata$ps, simdata$ID, FUN = function(X) 1/cumprod(X))
      simdata$tall <- simdata$t
      simdata$sub <- simdata$ID
      simdata$RA <- 1
      simdata[simdata$t == 0 & !(simdata$CA == 1),]$RA <- 0
      simdata[simdata$t == 1 & !(simdata$CA == 2),]$RA <- 0
      simdata[simdata$t == 2 & !(simdata$CA == 3),]$RA <- 0
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
      
      wideSimdata$g_treat_1_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_1, wideSimdata$X2_1,wideSimdata$X3_1, wideSimdata$X4_1)) %*% treatment_model_pooled$coefficients)
      wideSimdata$g_treat_2_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_2, wideSimdata$X2_2,wideSimdata$X3_2, wideSimdata$X4_2)) %*% treatment_model_pooled$coefficients)
      
      wideSimdata$g_control_1_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_1, wideSimdata$X2_1,wideSimdata$X3_1, wideSimdata$X4_1)) %*% treatment_model_pooled$coefficients)
      wideSimdata$g_control_2_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_2, wideSimdata$X2_2,wideSimdata$X3_2, wideSimdata$X4_2)) %*% treatment_model_pooled$coefficients)
      
      a = min(simdata$Y,na.rm = TRUE); b = max(simdata$Y, na.rm = TRUE)
      
      wideSimdata$Y_2_scaled = (wideSimdata$Y_2-a)/(b-a)
      wideSimdata$Y_1_scaled = (wideSimdata$Y_1-a)/(b-a)
      wideSimdata$Y_0_scaled = (wideSimdata$Y_0-a)/(b-a)
      
      
      
      ########### t = 2 ##########
      
      Q2_2_fit <- glm(data = wideSimdata, formula = Y_2 ~ A_2 + A_1 + X1_2 + X2_2 + X3_2 +X4_2 + X1_1 + X2_1 + X3_1 +X4_1)
      
      logitQ2_2 <- qlogis((predict.glm(Q2_2_fit, newdata = data.frame(A_2 = c(rep(1,sample_size), rep(0,sample_size)),
                                                                      A_1 = c(rep(1, sample_size), rep(0, sample_size)),
                                                                      X1_2 = rep(wideSimdata$X1_2, 2), 
                                                                      X2_2 = rep(wideSimdata$X2_2, 2),
                                                                      X3_2 = rep(wideSimdata$X3_2, 2), 
                                                                      X4_2 = rep(wideSimdata$X4_2, 2),
                                                                      X1_1 = rep(wideSimdata$X1_1, 2), 
                                                                      X2_1 = rep(wideSimdata$X2_1, 2),
                                                                      X3_1 = rep(wideSimdata$X3_1, 2), 
                                                                      X4_1 = rep(wideSimdata$X4_1, 2)))-a)/(b-a))
      
      #------------- update Q2_2-----------------------
      regimen_ind <- c(as.numeric(wideSimdata$CA_2 ==3), as.numeric(wideSimdata$CA_2 ==0))
      weight <- rep(wideSimdata$weights_2,2)
      Cweight <- rep(wideSimdata$Cweights_2,2)
      Cweight_aggr <- rep(wideSimdata$Cweights_aggr_2,2)
      
      update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                Y_2 = rep(wideSimdata$Y_2_scaled, 2), 
                                off = logitQ2_2,
                                intercept = rep(1,sample_size*2), 
                                cumA = c(rep(3,sample_size), rep(0,sample_size)))
      
      Q2_2star_fit <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                          data = update_data[regimen_ind == 1,],
                          weights = as.vector(weight[regimen_ind == 1]))
      
      Q2_2star <- predict.glm(Q2_2star_fit, newdata = update_data, type = 'response')*(b-a) +a
      
      Q2_2star_fit_cali <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q2_2star_cali <- predict.glm(Q2_2star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
      
      Q2_2star_fit_cali_aggr <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q2_2star_cali_aggr <- predict.glm(Q2_2star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
      
      ########### t = 1 ##########
      #------------- Get Q2_1, Q1_1 ---------------------
      fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                                 Q2_2star = Q2_2star, 
                                 Q2_2star_cali = Q2_2star_cali,
                                 Q2_2star_cali_aggr = Q2_2star_cali_aggr,
                                 CA_1 = rep(wideSimdata$CA_1, 2),
                                 A_1 = rep(wideSimdata$A_1, 2),
                                 A_0 = rep(wideSimdata$A_0, 2),
                                 X1_1 = rep(wideSimdata$X1_1, 2), 
                                 X2_1 = rep(wideSimdata$X2_1, 2),
                                 X3_1 = rep(wideSimdata$X3_1, 2), 
                                 X4_1 = rep(wideSimdata$X4_1, 2))
      
      Q2_1_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1)
      Q2_1_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1)
      
      logitQ2_1 <- as.matrix(c(qlogis((predict.glm(Q2_1_fit_d1, 
                                                   newdata = data.frame(CA_1 = rep(2,sample_size),
                                                                        X1_1 = wideSimdata$X1_1, 
                                                                        X2_1 = wideSimdata$X2_1,
                                                                        X3_1 = wideSimdata$X3_1, 
                                                                        X4_1 = wideSimdata$X4_1))-a)/(b-a)),
                               qlogis((predict.glm(Q2_1_fit_d0, 
                                                   newdata = data.frame(CA_1 = rep(0,sample_size),
                                                                        X1_1 = wideSimdata$X1_1, 
                                                                        X2_1 = wideSimdata$X2_1,
                                                                        X3_1 = wideSimdata$X3_1, 
                                                                        X4_1 = wideSimdata$X4_1))-a)/(b-a))
      ))
      
      Q2_1_fit_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star_cali ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1)
      Q2_1_fit_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star_cali ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1)
      
      logitQ2_1_cali <- as.matrix(c(qlogis((predict.glm(Q2_1_fit_d1_cali, 
                                                        newdata = data.frame(CA_1 = rep(2,sample_size),
                                                                             X1_1 = wideSimdata$X1_1, 
                                                                             X2_1 = wideSimdata$X2_1,
                                                                             X3_1 = wideSimdata$X3_1, 
                                                                             X4_1 = wideSimdata$X4_1))-a)/(b-a)),
                                    qlogis((predict.glm(Q2_1_fit_d0_cali, 
                                                        newdata = data.frame(CA_1 = rep(0,sample_size),
                                                                             X1_1 = wideSimdata$X1_1, 
                                                                             X2_1 = wideSimdata$X2_1,
                                                                             X3_1 = wideSimdata$X3_1, 
                                                                             X4_1 = wideSimdata$X4_1))-a)/(b-a))
      ))
      
      Q2_1_fit_d1_cali_aggr <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star_cali_aggr ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1)
      Q2_1_fit_d0_cali_aggr <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star_cali_aggr ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1)
      
      logitQ2_1_cali_aggr <- as.matrix(c(qlogis((predict.glm(Q2_1_fit_d1_cali_aggr, 
                                                             newdata = data.frame(CA_1 = rep(2,sample_size),
                                                                                  X1_1 = wideSimdata$X1_1, 
                                                                                  X2_1 = wideSimdata$X2_1,
                                                                                  X3_1 = wideSimdata$X3_1, 
                                                                                  X4_1 = wideSimdata$X4_1))-a)/(b-a)),
                                         qlogis((predict.glm(Q2_1_fit_d0_cali_aggr, 
                                                             newdata = data.frame(CA_1 = rep(0,sample_size),
                                                                                  X1_1 = wideSimdata$X1_1, 
                                                                                  X2_1 = wideSimdata$X2_1,
                                                                                  X3_1 = wideSimdata$X3_1, 
                                                                                  X4_1 = wideSimdata$X4_1))-a)/(b-a))
      ))
      
      Q1_1_fit <- glm(data = wideSimdata, formula = Y_1 ~ A_1 + A_0 + X1_1 + X2_1 + X3_1 + X4_1 + X1_0 + X2_0 + X3_0 + X4_0)
      
      logitQ1_1 <- qlogis((predict.glm(Q1_1_fit, newdata = data.frame(A_1 = c(rep(1,sample_size), rep(0,sample_size)),
                                                                      A_0 = c(rep(1,sample_size), rep(0,sample_size)),
                                                                      X1_1 = rep(wideSimdata$X1_1, 2), 
                                                                      X2_1 = rep(wideSimdata$X2_1, 2),
                                                                      X3_1 = rep(wideSimdata$X3_1, 2), 
                                                                      X4_1 = rep(wideSimdata$X4_1, 2),
                                                                      X1_0 = rep(wideSimdata$X1_0, 2), 
                                                                      X2_0 = rep(wideSimdata$X2_0, 2),
                                                                      X3_0 = rep(wideSimdata$X3_0, 2), 
                                                                      X4_0 = rep(wideSimdata$X4_0, 2)))-a)/(b-a))
      
      
      #------------- update Q1_1-----------------------
      regimen_ind = c(as.numeric(wideSimdata$CA_1 ==2), as.numeric(wideSimdata$CA_1 ==0))
      weight <- rep(wideSimdata$weights_1,2)
      Cweight <- rep(wideSimdata$Cweights_1,2)
      Cweight_aggr <- rep(wideSimdata$Cweights_aggr_1,2)
      
      update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                Y_1 = rep(wideSimdata$Y_1_scaled, 2), 
                                off = logitQ1_1,
                                intercept = rep(1,sample_size*2), 
                                cumA = c(rep(2,sample_size), rep(0,sample_size)))
      
      Q1_1star_fit <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                          data = update_data[regimen_ind == 1,],
                          weights = as.vector(weight[regimen_ind == 1]))
      
      Q1_1star <- predict.glm(Q1_1star_fit, newdata = update_data, type = 'response')*(b-a) +a
      
      Q1_1star_fit_cali <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q1_1star_cali <- predict.glm(Q1_1star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
      
      Q1_1star_fit_cali_aggr <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q1_1star_cali_aggr <- predict.glm(Q1_1star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
      
      #------------- update Q2_1-----------------------
      
      update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                Q2_2star = (Q2_2star-a)/(b-a),
                                Q2_2star_cali = (Q2_2star_cali-a)/(b-a),
                                Q2_2star_cali_aggr = (Q2_2star_cali_aggr-a)/(b-a),
                                off = logitQ2_1,
                                off_cali = logitQ2_1_cali,
                                off_cali_aggr = logitQ2_1_cali_aggr,
                                intercept = rep(1,sample_size*2), 
                                cumA = c(rep(3,sample_size), rep(0,sample_size)))
      
      Q2_1star_fit <- glm(Q2_2star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                          data = update_data[regimen_ind == 1,],
                          weights = as.vector(weight[regimen_ind == 1]))
      
      Q2_1star <- predict.glm(Q2_1star_fit, newdata = update_data, type = 'response')*(b-a) +a
      
      Q2_1star_fit_cali <- glm(Q2_2star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q2_1star_cali <- predict.glm(Q2_1star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
      
      Q2_1star_fit_cali_aggr <- glm(Q2_2star_cali_aggr ~ intercept + cumA + offset(off_cali_aggr) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q2_1star_cali_aggr <- predict.glm(Q2_1star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
      
      ########### t = 0 ##########
      #------------- Get Q4_0, Q3_0, Q2_0, Q1_0, Q0_0 ---------------------
      fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                                 Q2_1star = Q2_1star, 
                                 Q1_1star = Q1_1star, 
                                 Q2_1star_cali = Q2_1star_cali, 
                                 Q2_1star_cali_aggr = Q2_1star_cali_aggr,
                                 Q1_1star_cali = Q1_1star_cali, 
                                 Q1_1star_cali_aggr = Q1_1star_cali_aggr, 
                                 A_0 = rep(wideSimdata$A_0, 2),
                                 X1_0 = rep(wideSimdata$X1_0, 2), 
                                 X2_0 = rep(wideSimdata$X2_0, 2),
                                 X3_0 = rep(wideSimdata$X3_0, 2), 
                                 X4_0 = rep(wideSimdata$X4_0, 2))
      
      
      Q2_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star ~ A_0)
      Q2_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star ~ A_0)
      
      logitQ2_0 <- as.matrix(c(qlogis((predict.glm(Q2_0_fit_d1, 
                                                   newdata = data.frame(A_0 = rep(1,sample_size)))-a)/(b-a)),
                               qlogis((predict.glm(Q2_0_fit_d0, 
                                                   newdata = data.frame(A_0 = rep(0,sample_size)))-a)/(b-a))
                               ))
      
      Q2_0_fit_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star_cali ~ A_0)
      Q2_0_fit_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star_cali ~ A_0)
      
      logitQ2_0_cali <- as.matrix(c(qlogis((predict.glm(Q2_0_fit_d1_cali, 
                                                        newdata = data.frame(A_0 = rep(1,sample_size)))-a)/(b-a)),
                                    qlogis((predict.glm(Q2_0_fit_d0_cali, 
                                                        newdata = data.frame(A_0 = rep(0,sample_size)))-a)/(b-a))))
      
      Q2_0_fit_d1_cali_aggr <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star_cali_aggr ~ A_0)
      Q2_0_fit_d0_cali_aggr <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star_cali_aggr ~ A_0)
      
      logitQ2_0_cali_aggr <- as.matrix(c(qlogis((predict.glm(Q2_0_fit_d1_cali_aggr, 
                                                             newdata = data.frame(A_0 = rep(1,sample_size)))-a)/(b-a)),
                                         qlogis((predict.glm(Q2_0_fit_d0_cali_aggr, 
                                                             newdata = data.frame(A_0 = rep(0,sample_size)))-a)/(b-a))))
      
      Q1_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0)
      Q1_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0)
      
      logitQ1_0 <- as.matrix(c(qlogis((predict.glm(Q1_0_fit_d1, 
                                                   newdata = data.frame(A_0 = rep(1,sample_size),
                                                                        X1_0 = wideSimdata$X1_0, 
                                                                        X2_0 = wideSimdata$X2_0,
                                                                        X3_0 = wideSimdata$X3_0, 
                                                                        X4_0 = wideSimdata$X4_0))-a)/(b-a)),
                               qlogis((predict.glm(Q1_0_fit_d0, 
                                                   newdata = data.frame(A_0 = rep(0,sample_size),
                                                                        X1_0 = wideSimdata$X1_0, 
                                                                        X2_0 = wideSimdata$X2_0,
                                                                        X3_0 = wideSimdata$X3_0, 
                                                                        X4_0 = wideSimdata$X4_0))-a)/(b-a))))
      
      Q1_0_fit_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star_cali ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0 )
      Q1_0_fit_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star_cali ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0)
      
      logitQ1_0_cali <- as.matrix(c(qlogis((predict.glm(Q1_0_fit_d1_cali, 
                                                        newdata = data.frame(A_0 = rep(1,sample_size),
                                                                             X1_0 = wideSimdata$X1_0, 
                                                                             X2_0 = wideSimdata$X2_0,
                                                                             X3_0 = wideSimdata$X3_0, 
                                                                             X4_0 = wideSimdata$X4_0))-a)/(b-a)),
                                    qlogis((predict.glm(Q1_0_fit_d0_cali, 
                                                        newdata = data.frame(A_0 = rep(0,sample_size),
                                                                             X1_0 = wideSimdata$X1_0, 
                                                                             X2_0 = wideSimdata$X2_0,
                                                                             X3_0 = wideSimdata$X3_0, 
                                                                             X4_0 = wideSimdata$X4_0))-a)/(b-a))
      ))
      
      Q1_0_fit_d1_cali_aggr <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star_cali_aggr ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0 )
      Q1_0_fit_d0_cali_aggr <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star_cali_aggr ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0)
      
      logitQ1_0_cali_aggr <- as.matrix(c(qlogis((predict.glm(Q1_0_fit_d1_cali_aggr, 
                                                             newdata = data.frame(A_0 = rep(1,sample_size),
                                                                                  X1_0 = wideSimdata$X1_0, 
                                                                                  X2_0 = wideSimdata$X2_0,
                                                                                  X3_0 = wideSimdata$X3_0, 
                                                                                  X4_0 = wideSimdata$X4_0))-a)/(b-a)),
                                         qlogis((predict.glm(Q1_0_fit_d0_cali_aggr, 
                                                             newdata = data.frame(A_0 = rep(0,sample_size),
                                                                                  X1_0 = wideSimdata$X1_0, 
                                                                                  X2_0 = wideSimdata$X2_0,
                                                                                  X3_0 = wideSimdata$X3_0, 
                                                                                  X4_0 = wideSimdata$X4_0))-a)/(b-a))))
      
      Q0_0_fit <- glm(data = wideSimdata, formula = Y_0 ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0)
      
      logitQ0_0 <- qlogis((predict.glm(Q0_0_fit, newdata = data.frame(A_0 = c(rep(1,sample_size), rep(0,sample_size)), 
                                                                      X1_0 = rep(wideSimdata$X1_0, 2), 
                                                                      X2_0 = rep(wideSimdata$X2_0, 2),
                                                                      X3_0 = rep(wideSimdata$X3_0, 2), 
                                                                      X4_0 = rep(wideSimdata$X4_0, 2)))-a)/(b-a))
      
      
      #------------- update Q0_0-----------------------
      regimen_ind = c(as.numeric(wideSimdata$CA_0 ==1), as.numeric(wideSimdata$CA_0 ==0))
      weight <- rep(wideSimdata$weights_0,2)
      Cweight <- rep(wideSimdata$Cweights_0,2)
      Cweight_aggr <- rep(wideSimdata$Cweights_aggr_0,2)
      
      update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                Y_0 = rep(wideSimdata$Y_0_scaled, 2), 
                                off = logitQ0_0,
                                intercept = rep(1,sample_size*2), 
                                cumA = c(rep(1,sample_size), rep(0,sample_size)))
      
      Q0_0star_fit <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                          data = update_data[regimen_ind == 1,],
                          weights = as.vector(weight[regimen_ind == 1]))
      
      Q0_0star <- predict.glm(Q0_0star_fit, newdata = update_data, type = 'response')*(b-a) +a
      
      Q0_0star_fit_cali <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q0_0star_cali <- predict.glm(Q0_0star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
      
      Q0_0star_fit_cali_aggr <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q0_0star_cali_aggr <- predict.glm(Q0_0star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
      
      #------------- update Q2_0-----------------------
      
      update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                Q2_1star = (Q2_1star-a)/(b-a),
                                Q2_1star_cali = (Q2_1star_cali-a)/(b-a),
                                Q2_1star_cali_aggr = (Q2_1star_cali_aggr-a)/(b-a),
                                off = logitQ2_0,
                                off_cali = logitQ2_0_cali,
                                off_cali_aggr = logitQ2_0_cali_aggr,
                                intercept = rep(1,sample_size*2), 
                                cumA = c(rep(3,sample_size), rep(0,sample_size)))
      
      Q2_0star_fit <- glm(Q2_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                          data = update_data[regimen_ind == 1,],
                          weights = as.vector(weight[regimen_ind == 1]))
      
      Q2_0star <- predict.glm(Q2_0star_fit, newdata = update_data, type = 'response')*(b-a) +a
      
      Q2_0star_fit_cali <- glm(Q2_1star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q2_0star_cali <- predict.glm(Q2_0star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
      
      Q2_0star_fit_cali_aggr <- glm(Q2_1star_cali_aggr ~ intercept + cumA + offset(off_cali_aggr) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q2_0star_cali_aggr <- predict.glm(Q2_0star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
      
      #------------- update Q1_0-----------------------
      
      update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                Q1_1star = (Q1_1star-a)/(b-a),
                                Q1_1star_cali = (Q1_1star_cali-a)/(b-a),
                                Q1_1star_cali_aggr = (Q1_1star_cali_aggr-a)/(b-a),
                                off = logitQ1_0,
                                off_cali = logitQ1_0_cali,
                                off_cali_aggr = logitQ1_0_cali_aggr,
                                intercept = rep(1,sample_size*2), 
                                cumA = c(rep(2,sample_size), rep(0,sample_size)))
      
      Q1_0star_fit <- glm(Q1_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                          data = update_data[regimen_ind == 1,],
                          weights = as.vector(weight[regimen_ind == 1]))
      
      Q1_0star <- predict.glm(Q1_0star_fit, newdata = update_data, type = 'response')*(b-a) +a
      
      Q1_0star_fit_cali <- glm(Q1_1star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q1_0star_cali <- predict.glm(Q1_0star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
      
      Q1_0star_fit_cali_aggr <- glm(Q1_1star_cali_aggr ~ intercept + cumA + offset(off_cali_aggr) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q1_0star_cali_aggr <- predict.glm(Q1_0star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
      
      ################# Fit MSM ##########################
      
      msm_fitting_data <- data.frame(id = rep(1:sample_size,6), 
                                     t = c(rep(0,sample_size*2), 
                                           rep(1,sample_size*2), 
                                           rep(2,sample_size*2)),
                                     Y = c(Q0_0star, Q1_0star, Q2_0star), 
                                     Y_cali = c(Q0_0star_cali, Q1_0star_cali, Q2_0star_cali),
                                     Y_cali_aggr = c(Q0_0star_cali_aggr, Q1_0star_cali_aggr, Q2_0star_cali_aggr),
                                     cumA = c(rep(1,sample_size), rep(0,sample_size), 
                                              rep(2,sample_size), rep(0,sample_size), 
                                              rep(3,sample_size), rep(0,sample_size)))
      
      msm <- glm(Y ~ cumA, data = msm_fitting_data)
      msm_cali <- glm(Y_cali ~ cumA, data = msm_fitting_data)
      msm_cali_aggr <- glm(Y_cali_aggr ~ cumA, data = msm_fitting_data)
      
      
      ##################### IPW - MSM ##################################
      
      simdata$RA <- 1
      simdata[simdata$t == 0 & !(simdata$CA == 1 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 1 & !(simdata$CA == 2 | simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 2 & !(simdata$CA == 3 | simdata$CA == 0),]$RA <- 0
      
      switch_data <- simdata[simdata$RA == 1,]
      
      switch_data$X1 <- ave(switch_data$X1, switch_data$ID, FUN = first)
      switch_data$X2 <- ave(switch_data$X2, switch_data$ID, FUN = first)
      switch_data$X3 <- ave(switch_data$X3, switch_data$ID, FUN = first)
      switch_data$X4 <- ave(switch_data$X4, switch_data$ID, FUN = first)
      
      switch_data$weight <- 1.0
      switch_data[switch_data$t == 1& switch_data$A == 1,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 1,],]$g_treat_1_pooled)
      switch_data[switch_data$t == 2& switch_data$A == 1,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 1,],]$g_treat_1_pooled*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 1,],]$g_treat_2_pooled)
      
      switch_data[switch_data$t == 1& switch_data$A == 0,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 0,],]$g_control_1_pooled)
      switch_data[switch_data$t == 2& switch_data$A == 0,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 0,],]$g_control_1_pooled*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 0,],]$g_control_2_pooled)
      
      PP_pooled <- glm(Y ~ CA + X1 + X2 + X3 + X4, data = switch_data, weights = weight, family = 'gaussian')
      PP_strat <- glm(Y ~ CA, data = switch_data, weights = weights, family = 'gaussian')
      PP_cali <- glm(Y ~ CA, data = switch_data, weights = Cweights, family = 'gaussian')
      PP_cali_aggr <- glm(Y ~ CA, data = switch_data, weights = Cweights_aggr, family = 'gaussian')
      
    }))
    
    fitting_data_CA3 <- switch_data[switch_data$t == 0,]
    fitting_data_CA3$CA <- 3
    
    fitting_data_CA0 <- switch_data[switch_data$t == 0,]
    fitting_data_CA0$CA <- 0
    
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
  
  return(results)
}


iters = 1000
registerDoParallel(cores = 10)


cat(paste('Table 4 results \n'))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 300, conf = 0.2, seeds = seeds),
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 300, conf = 0.2, seeds = seeds)),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 500, conf = 0.2, seeds = seeds),
             simulation_code(iters = iters, transformed = TRUE, sample_size = 500, conf = 0.2, seeds = seeds)),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 1000, conf = 0.2, seeds = seeds),
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 1000, conf = 0.2, seeds = seeds)),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 2500, conf = 0.2, seeds = seeds),
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 2500, conf = 0.2, seeds = seeds)),
             type = "latex"))
cat(paste('%-------------------------------------------- \n'))
cat(paste('Table 5 results \n'))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 300,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1),
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 300,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1),
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 1000,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1),
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 1000,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 2500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1),
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 2500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)),
             type = "latex"))


