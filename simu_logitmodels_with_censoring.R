#!/usr/bin R
###############################
# This script runs simulations with Q functions on logit form
# Jointly calibrated weights for treatment switching AND dropout censoring
###############################
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/dgm_biometrics_censoring.R")
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
      simdata<-DATA_GEN_censored(ns = sample_size, conf = conf,treat_prev_0 = treat_prev_0, 
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
      
      censor_model_0 <- glm(C ~ X1 + X2 + X3 + X4, data = simdata[simdata$t == 0,], family = 'binomial')
      censor_model_1 <- glm(C ~ X1 + X2 + X3 + X4, data = simdata[simdata$t == 1,], family = 'binomial')
      
      simdata$lagX1 <- ave(simdata$X1, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
      simdata$lagX2 <- ave(simdata$X2, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
      simdata$lagX3 <- ave(simdata$X3, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
      simdata$lagX4 <- ave(simdata$X4, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
      
      simdata$pr_c <- 1.0
      simdata[simdata$t == 1,]$pr_c <- as.vector(1-predict.glm(censor_model_0, 
                                                               newdata = data.frame(X1 = simdata[simdata$t == 1,]$lagX1,
                                                                                    X2 = simdata[simdata$t == 1,]$lagX2,
                                                                                    X3 = simdata[simdata$t == 1,]$lagX3,
                                                                                    X4 = simdata[simdata$t == 1,]$lagX4),
                                                               type = 'response'))
      simdata[simdata$t == 2,]$pr_c <- as.vector(1-predict.glm(censor_model_1,
                                                               newdata = data.frame(X1 = simdata[simdata$t == 2,]$lagX1,
                                                                                    X2 = simdata[simdata$t == 2,]$lagX2,
                                                                                    X3 = simdata[simdata$t == 2,]$lagX3,
                                                                                    X4 = simdata[simdata$t == 2,]$lagX4),
                                                               type = 'response'))
      
      simdata$treat_weight <- ave(simdata$ps, simdata$ID, FUN = function(X) 1/cumprod(X))
      simdata$censor_weight <- ave(simdata$pr_c, simdata$ID, FUN = function(X) 1/cumprod(X))
      simdata$weight <- simdata$treat_weight*simdata$censor_weight
      
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
      
      
      calibrate_always_treated <- calibration_by_time_from_baseline(simdata, 
                                                                    var = c("X1", "X2", "X3", "X4"), 
                                                                    censor = TRUE, 
                                                                    c_var = c("X1", "X2", "X3", "X4"),
                                                                    weights_var = 'weight')
      calibrate_always_treated_aggr <- aggregated_calibration_from_baseline(simdata, 
                                                                            var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), 
                                                                            censor = TRUE, 
                                                                            c_var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"),
                                                                            weights_var = 'weight')
      
      simdata$Cweights <- calibrate_always_treated$data$Cweights
      simdata$Cweights_aggr <- calibrate_always_treated_aggr$data$Cweights
      
      simdata$RA <- 1
      simdata[simdata$t == 0 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 1 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 2 & !(simdata$CA == 0),]$RA <- 0
      
      calibrate_never_treated <- calibration_by_time_from_baseline(simdata, 
                                                                   var = c("X1", "X2", "X3", "X4"), 
                                                                   censor = TRUE, 
                                                                   c_var = c("X1", "X2", "X3", "X4"),
                                                                   weights_var = 'Cweights')
      calibrate_never_treated_aggr <- aggregated_calibration_from_baseline(simdata, 
                                                                           var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), 
                                                                           censor = TRUE, 
                                                                           c_var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"),
                                                                           weights_var = 'Cweights_aggr')
      
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
      wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, value.var = c("A", "X1", "X2","X3", "X4", "CA", "Y","C", "weights", "Cweights", "Cweights_aggr"))
      
      wideSimdata$g_treat_1_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_1, wideSimdata$X2_1,wideSimdata$X3_1, wideSimdata$X4_1)) %*% treatment_model_pooled$coefficients)
      wideSimdata$g_treat_2_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_2, wideSimdata$X2_2,wideSimdata$X3_2, wideSimdata$X4_2)) %*% treatment_model_pooled$coefficients)
      
      wideSimdata$g_control_1_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_1, wideSimdata$X2_1,wideSimdata$X3_1, wideSimdata$X4_1)) %*% treatment_model_pooled$coefficients)
      wideSimdata$g_control_2_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_2, wideSimdata$X2_2,wideSimdata$X3_2, wideSimdata$X4_2)) %*% treatment_model_pooled$coefficients)
      
      a = min(simdata$Y,na.rm = TRUE); b = max(simdata$Y, na.rm = TRUE)
      
      wideSimdata$Y_2_scaled = (wideSimdata$Y_2-a)/(b-a)
      wideSimdata$Y_1_scaled = (wideSimdata$Y_1-a)/(b-a)
      wideSimdata$Y_0_scaled = (wideSimdata$Y_0-a)/(b-a)
      
      wideSimdata$RC_0 <- 1.0
      wideSimdata$RC_1 <- 1-wideSimdata$C_0
      wideSimdata$RC_2 <- ifelse(wideSimdata$RC_1 == 0, 0, 1-wideSimdata$C_1 )
      
      ########### t = 2 ##########
      
      Q2_2_fit <- glm(data = wideSimdata[wideSimdata$RC_2 ==1,], formula = Y_2_scaled ~ A_2 + A_1 + X1_2 + X2_2 + X3_2 +X4_2 + X1_1 + X2_1 + X3_1 +X4_1, family = 'quasibinomial')
      
      logitQ2_2 <- predict.glm(Q2_2_fit, newdata = data.frame(A_2 = c(rep(1,sample_size), rep(0,sample_size)),
                                                              A_1 = c(rep(1, sample_size), rep(0, sample_size)),
                                                              X1_2 = rep(wideSimdata$X1_2, 2), 
                                                              X2_2 = rep(wideSimdata$X2_2, 2),
                                                              X3_2 = rep(wideSimdata$X3_2, 2), 
                                                              X4_2 = rep(wideSimdata$X4_2, 2),
                                                              X1_1 = rep(wideSimdata$X1_1, 2), 
                                                              X2_1 = rep(wideSimdata$X2_1, 2),
                                                              X3_1 = rep(wideSimdata$X3_1, 2), 
                                                              X4_1 = rep(wideSimdata$X4_1, 2)),type = 'link')
      
      #------------- update Q2_2-----------------------
      regimen_ind <- c(as.numeric(wideSimdata$CA_2 ==3), as.numeric(wideSimdata$CA_2 ==0))
      regimen_ind[is.na(regimen_ind)] <- 0
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
      
      Q2_2star <- predict.glm(Q2_2star_fit, newdata = update_data, type = 'response')
      
      Q2_2star_fit_cali <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q2_2star_cali <- predict.glm(Q2_2star_fit_cali, newdata = update_data, type = 'response')
      
      Q2_2star_fit_cali_aggr <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q2_2star_cali_aggr <- predict.glm(Q2_2star_fit_cali_aggr, newdata = update_data, type = 'response')
      
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
      
      Q2_1_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
      Q2_1_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
      
      logitQ2_1 <- as.matrix(c(predict.glm(Q2_1_fit_d1, 
                                           newdata = data.frame(CA_1 = rep(2,sample_size),
                                                                X1_1 = wideSimdata$X1_1, 
                                                                X2_1 = wideSimdata$X2_1,
                                                                X3_1 = wideSimdata$X3_1, 
                                                                X4_1 = wideSimdata$X4_1), 
                                           type = 'link'),
                               predict.glm(Q2_1_fit_d0, 
                                           newdata = data.frame(CA_1 = rep(0,sample_size),
                                                                X1_1 = wideSimdata$X1_1, 
                                                                X2_1 = wideSimdata$X2_1,
                                                                X3_1 = wideSimdata$X3_1, 
                                                                X4_1 = wideSimdata$X4_1), 
                                           type = 'link')))
      
      Q2_1_fit_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star_cali ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
      Q2_1_fit_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star_cali ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
      
      logitQ2_1_cali <- as.matrix(c(predict.glm(Q2_1_fit_d1_cali, 
                                                newdata = data.frame(CA_1 = rep(2,sample_size),
                                                                     X1_1 = wideSimdata$X1_1, 
                                                                     X2_1 = wideSimdata$X2_1,
                                                                     X3_1 = wideSimdata$X3_1, 
                                                                     X4_1 = wideSimdata$X4_1), 
                                                type = 'link'),
                                    predict.glm(Q2_1_fit_d0_cali, 
                                                newdata = data.frame(CA_1 = rep(0,sample_size),
                                                                     X1_1 = wideSimdata$X1_1, 
                                                                     X2_1 = wideSimdata$X2_1,
                                                                     X3_1 = wideSimdata$X3_1, 
                                                                     X4_1 = wideSimdata$X4_1), 
                                                type = 'link')))
      
      Q2_1_fit_d1_cali_aggr <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star_cali_aggr ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
      Q2_1_fit_d0_cali_aggr <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star_cali_aggr ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
      
      logitQ2_1_cali_aggr <- as.matrix(c(predict.glm(Q2_1_fit_d1_cali_aggr, 
                                                     newdata = data.frame(CA_1 = rep(2,sample_size),
                                                                          X1_1 = wideSimdata$X1_1, 
                                                                          X2_1 = wideSimdata$X2_1,
                                                                          X3_1 = wideSimdata$X3_1, 
                                                                          X4_1 = wideSimdata$X4_1), 
                                                     type = 'link'),
                                         predict.glm(Q2_1_fit_d0_cali_aggr, 
                                                     newdata = data.frame(CA_1 = rep(0,sample_size),
                                                                          X1_1 = wideSimdata$X1_1, 
                                                                          X2_1 = wideSimdata$X2_1,
                                                                          X3_1 = wideSimdata$X3_1, 
                                                                          X4_1 = wideSimdata$X4_1), 
                                                     type = 'link')))
      
      Q1_1_fit <- glm(data = wideSimdata, formula = Y_1_scaled ~ A_1 + A_0 + X1_1 + X2_1 + X3_1 + X4_1 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
      
      logitQ1_1 <- predict.glm(Q1_1_fit, newdata = data.frame(A_1 = c(rep(1,sample_size), rep(0,sample_size)),
                                                              A_0 = c(rep(1,sample_size), rep(0,sample_size)),
                                                              X1_1 = rep(wideSimdata$X1_1, 2), 
                                                              X2_1 = rep(wideSimdata$X2_1, 2),
                                                              X3_1 = rep(wideSimdata$X3_1, 2), 
                                                              X4_1 = rep(wideSimdata$X4_1, 2),
                                                              X1_0 = rep(wideSimdata$X1_0, 2), 
                                                              X2_0 = rep(wideSimdata$X2_0, 2),
                                                              X3_0 = rep(wideSimdata$X3_0, 2), 
                                                              X4_0 = rep(wideSimdata$X4_0, 2)), type = 'link')
      
      
      #------------- update Q1_1-----------------------
      regimen_ind = c(as.numeric(wideSimdata$CA_1 ==2), as.numeric(wideSimdata$CA_1 ==0))
      regimen_ind[is.na(regimen_ind)] <- 0
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
      
      Q1_1star <- predict.glm(Q1_1star_fit, newdata = update_data, type = 'response')
      
      Q1_1star_fit_cali <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q1_1star_cali <- predict.glm(Q1_1star_fit_cali, newdata = update_data, type = 'response')
      
      Q1_1star_fit_cali_aggr <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q1_1star_cali_aggr <- predict.glm(Q1_1star_fit_cali_aggr, newdata = update_data, type = 'response')
      
      #------------- update Q2_1-----------------------
      
      update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                Q2_2star = Q2_2star,
                                Q2_2star_cali = Q2_2star_cali,
                                Q2_2star_cali_aggr = Q2_2star_cali_aggr,
                                off = logitQ2_1,
                                off_cali = logitQ2_1_cali,
                                off_cali_aggr = logitQ2_1_cali_aggr,
                                intercept = rep(1,sample_size*2), 
                                cumA = c(rep(3,sample_size), rep(0,sample_size)))
      
      Q2_1star_fit <- glm(Q2_2star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                          data = update_data[regimen_ind == 1,],
                          weights = as.vector(weight[regimen_ind == 1]))
      
      Q2_1star <- predict.glm(Q2_1star_fit, newdata = update_data, type = 'response')
      
      Q2_1star_fit_cali <- glm(Q2_2star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q2_1star_cali <- predict.glm(Q2_1star_fit_cali, newdata = update_data, type = 'response')
      
      Q2_1star_fit_cali_aggr <- glm(Q2_2star_cali_aggr ~ intercept + cumA + offset(off_cali_aggr) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q2_1star_cali_aggr <- predict.glm(Q2_1star_fit_cali_aggr, newdata = update_data, type = 'response')
      
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
      
      
      Q2_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star ~ A_0, family = 'quasibinomial')
      Q2_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star ~ A_0, family = 'quasibinomial')
      
      logitQ2_0 <- as.matrix(c(predict.glm(Q2_0_fit_d1, 
                                           newdata = data.frame(A_0 = rep(1,sample_size)),
                                           type = 'link'),
                               predict.glm(Q2_0_fit_d0, 
                                           newdata = data.frame(A_0 = rep(0,sample_size)),
                                           type = 'link')))
      
      Q2_0_fit_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star_cali ~ A_0, family = 'quasibinomial')
      Q2_0_fit_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star_cali ~ A_0, family = 'quasibinomial')
      
      logitQ2_0_cali <- as.matrix(c(predict.glm(Q2_0_fit_d1_cali, 
                                                newdata = data.frame(A_0 = rep(1,sample_size)),
                                                type = 'link'),
                                    predict.glm(Q2_0_fit_d0_cali,
                                                newdata = data.frame(A_0 = rep(0,sample_size)),
                                                type = 'link')))
      
      Q2_0_fit_d1_cali_aggr <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star_cali_aggr ~ A_0, family = 'quasibinomial')
      Q2_0_fit_d0_cali_aggr <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star_cali_aggr ~ A_0, family = 'quasibinomial')
      
      logitQ2_0_cali_aggr <- as.matrix(c(predict.glm(Q2_0_fit_d1_cali_aggr, 
                                                     newdata = data.frame(A_0 = rep(1,sample_size)),
                                                     type = 'link'),
                                         predict.glm(Q2_0_fit_d0_cali_aggr,
                                                     newdata = data.frame(A_0 = rep(0,sample_size)),
                                                     type = 'link')))
      
      Q1_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
      Q1_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
      
      logitQ1_0 <- as.matrix(c(predict.glm(Q1_0_fit_d1, 
                                           newdata = data.frame(A_0 = rep(1,sample_size),
                                                                X1_0 = wideSimdata$X1_0, 
                                                                X2_0 = wideSimdata$X2_0,
                                                                X3_0 = wideSimdata$X3_0, 
                                                                X4_0 = wideSimdata$X4_0),
                                           type = 'link'),
                               predict.glm(Q1_0_fit_d0, 
                                           newdata = data.frame(A_0 = rep(0,sample_size),
                                                                X1_0 = wideSimdata$X1_0, 
                                                                X2_0 = wideSimdata$X2_0,
                                                                X3_0 = wideSimdata$X3_0, 
                                                                X4_0 = wideSimdata$X4_0),
                                           type = 'link')))
      
      Q1_0_fit_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star_cali ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0 , family = 'quasibinomial')
      Q1_0_fit_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star_cali ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
      
      logitQ1_0_cali <- as.matrix(c(predict.glm(Q1_0_fit_d1_cali, 
                                                newdata = data.frame(A_0 = rep(1,sample_size),
                                                                     X1_0 = wideSimdata$X1_0, 
                                                                     X2_0 = wideSimdata$X2_0,
                                                                     X3_0 = wideSimdata$X3_0, 
                                                                     X4_0 = wideSimdata$X4_0),
                                                type = 'link'),
                                    predict.glm(Q1_0_fit_d0_cali, 
                                                newdata = data.frame(A_0 = rep(0,sample_size),
                                                                     X1_0 = wideSimdata$X1_0, 
                                                                     X2_0 = wideSimdata$X2_0,
                                                                     X3_0 = wideSimdata$X3_0, 
                                                                     X4_0 = wideSimdata$X4_0),
                                                type = 'link')))
      
      Q1_0_fit_d1_cali_aggr <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star_cali_aggr ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0 , family = 'quasibinomial')
      Q1_0_fit_d0_cali_aggr <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star_cali_aggr ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
      
      logitQ1_0_cali_aggr <- as.matrix(c(predict.glm(Q1_0_fit_d1_cali_aggr, 
                                                     newdata = data.frame(A_0 = rep(1,sample_size),
                                                                          X1_0 = wideSimdata$X1_0, 
                                                                          X2_0 = wideSimdata$X2_0,
                                                                          X3_0 = wideSimdata$X3_0, 
                                                                          X4_0 = wideSimdata$X4_0),
                                                     type = 'link'),
                                         predict.glm(Q1_0_fit_d0_cali_aggr, 
                                                     newdata = data.frame(A_0 = rep(0,sample_size),
                                                                          X1_0 = wideSimdata$X1_0, 
                                                                          X2_0 = wideSimdata$X2_0,
                                                                          X3_0 = wideSimdata$X3_0, 
                                                                          X4_0 = wideSimdata$X4_0),
                                                     type = 'link')))
      
      Q0_0_fit <- glm(data = wideSimdata, formula = Y_0_scaled ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
      
      logitQ0_0 <- predict.glm(Q0_0_fit, newdata = data.frame(A_0 = c(rep(1,sample_size), rep(0,sample_size)), 
                                                              X1_0 = rep(wideSimdata$X1_0, 2), 
                                                              X2_0 = rep(wideSimdata$X2_0, 2),
                                                              X3_0 = rep(wideSimdata$X3_0, 2), 
                                                              X4_0 = rep(wideSimdata$X4_0, 2)), type = 'link')
      
      
      #------------- update Q0_0-----------------------
      regimen_ind = c(as.numeric(wideSimdata$CA_0 ==1), as.numeric(wideSimdata$CA_0 ==0))
      regimen_ind[is.na(regimen_ind)] <- 0
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
      
      Q0_0star <- predict.glm(Q0_0star_fit, newdata = update_data, type = 'response')
      
      Q0_0star_fit_cali <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q0_0star_cali <- predict.glm(Q0_0star_fit_cali, newdata = update_data, type = 'response')
      
      Q0_0star_fit_cali_aggr <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q0_0star_cali_aggr <- predict.glm(Q0_0star_fit_cali_aggr, newdata = update_data, type = 'response')
      
      #------------- update Q2_0-----------------------
      
      update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                Q2_1star = Q2_1star,
                                Q2_1star_cali = Q2_1star_cali,
                                Q2_1star_cali_aggr = Q2_1star_cali_aggr,
                                off = logitQ2_0,
                                off_cali = logitQ2_0_cali,
                                off_cali_aggr = logitQ2_0_cali_aggr,
                                intercept = rep(1,sample_size*2), 
                                cumA = c(rep(3,sample_size), rep(0,sample_size)))
      
      Q2_0star_fit <- glm(Q2_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                          data = update_data[regimen_ind == 1,],
                          weights = as.vector(weight[regimen_ind == 1]))
      
      Q2_0star <- predict.glm(Q2_0star_fit, newdata = update_data, type = 'response')
      
      Q2_0star_fit_cali <- glm(Q2_1star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q2_0star_cali <- predict.glm(Q2_0star_fit_cali, newdata = update_data, type = 'response')
      
      Q2_0star_fit_cali_aggr <- glm(Q2_1star_cali_aggr ~ intercept + cumA + offset(off_cali_aggr) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q2_0star_cali_aggr <- predict.glm(Q2_0star_fit_cali_aggr, newdata = update_data, type = 'response')
      
      #------------- update Q1_0-----------------------
      
      update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                Q1_1star = Q1_1star,
                                Q1_1star_cali = Q1_1star_cali,
                                Q1_1star_cali_aggr = Q1_1star_cali_aggr,
                                off = logitQ1_0,
                                off_cali = logitQ1_0_cali,
                                off_cali_aggr = logitQ1_0_cali_aggr,
                                intercept = rep(1,sample_size*2), 
                                cumA = c(rep(2,sample_size), rep(0,sample_size)))
      
      Q1_0star_fit <- glm(Q1_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                          data = update_data[regimen_ind == 1,],
                          weights = as.vector(weight[regimen_ind == 1]))
      
      Q1_0star <- predict.glm(Q1_0star_fit, newdata = update_data, type = 'response')
      
      Q1_0star_fit_cali <- glm(Q1_1star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                               data = update_data[regimen_ind == 1,],
                               weights = as.vector(Cweight[regimen_ind == 1]))
      
      Q1_0star_cali <- predict.glm(Q1_0star_fit_cali, newdata = update_data, type = 'response')
      
      Q1_0star_fit_cali_aggr <- glm(Q1_1star_cali_aggr ~ intercept + cumA + offset(off_cali_aggr) - 1, family = 'quasibinomial', 
                                    data = update_data[regimen_ind == 1,],
                                    weights = as.vector(Cweight_aggr[regimen_ind == 1]))
      
      Q1_0star_cali_aggr <- predict.glm(Q1_0star_fit_cali_aggr, newdata = update_data, type = 'response')
      
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
      
      msm <- glm(Y ~ cumA, data = msm_fitting_data, family = 'quasibinomial')
      msm_cali <- glm(Y_cali ~ cumA, data = msm_fitting_data, family = 'quasibinomial')
      msm_cali_aggr <- glm(Y_cali_aggr ~ cumA, data = msm_fitting_data, family = 'quasibinomial')
      
      msm_fitting_data_transformed <- data.frame(id = rep(1:sample_size,6), 
                                                 t = c(rep(0,sample_size*2), 
                                                       rep(1,sample_size*2), 
                                                       rep(2,sample_size*2)),
                                                 Y = c(Q0_0star, Q1_0star, Q2_0star)*(b-a) + a, 
                                                 Y_cali = c(Q0_0star_cali, Q1_0star_cali, Q2_0star_cali)*(b-a) +a,
                                                 Y_cali_aggr = c(Q0_0star_cali_aggr, Q1_0star_cali_aggr, Q2_0star_cali_aggr)*(b-a)+a,
                                                 cumA = c(rep(1,sample_size), rep(0,sample_size), 
                                                          rep(2,sample_size), rep(0,sample_size), 
                                                          rep(3,sample_size), rep(0,sample_size)))
      
      msm_transformed <- glm(Y ~ cumA, data = msm_fitting_data_transformed)
      msm_transformed_cali <- glm(Y_cali ~ cumA, data = msm_fitting_data_transformed)
      msm_transformed_cali_aggr <- glm(Y_cali_aggr ~ cumA, data = msm_fitting_data_transformed)
      
      
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
      msm_transformed$coefficients[1],
      msm_transformed_cali$coefficients[1],
      msm_transformed_cali_aggr$coefficients[1],
      PP_strat$coefficients[2],
      PP_cali$coefficients[2],
      PP_cali_aggr$coefficients[2],
      msm_transformed$coefficients[2],
      msm_transformed_cali$coefficients[2],
      msm_transformed_cali_aggr$coefficients[2]
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

sink("3_visits_with_censoring_logit_models.txt")

cat(paste('Table 2 results \n'))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 300, conf = 0.2, seeds = seeds)$MSM,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 300, conf = 0.2, seeds = seeds)$MSM),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 500, conf = 0.2, seeds = seeds)$MSM,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 500, conf = 0.2, seeds = seeds)$MSM),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 1000, conf = 0.2, seeds = seeds)$MSM,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 1000, conf = 0.2, seeds = seeds)$MSM),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 2500, conf = 0.2, seeds = seeds)$MSM,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 2500, conf = 0.2, seeds = seeds)$MSM),
             type = "latex"))

cat(paste('E(Y_T) results \n'))

print(xtable(cbind(simulation_code(iters = iters, sample_size = 300, conf = 0.2, seeds = seeds)$EYT,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 300, conf = 0.2, seeds = seeds)$EYT),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 500, conf = 0.2, seeds = seeds)$EYT,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 500, conf = 0.2, seeds = seeds)$EYT),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 1000, conf = 0.2, seeds = seeds)$EYT,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 1000, conf = 0.2, seeds = seeds)$EYT),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 2500, conf = 0.2, seeds = seeds)$EYT,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 2500, conf = 0.2, seeds = seeds)$EYT),
             type = "latex"))
cat(paste('%-------------------------------------------- \n'))
cat(paste('Table 3 results \n'))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 300,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$MSM,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 300,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$MSM),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$MSM,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$MSM),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 1000,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$MSM,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 1000,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$MSM),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 2500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$MSM,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 2500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$MSM),
             type = "latex"))

cat(paste('E(Y_T) results \n'))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 300,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$EYT,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 300,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$EYT),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$EYT,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$EYT),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 1000,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$EYT,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 1000,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$EYT),
             type = "latex"))
print(xtable(cbind(simulation_code(iters = iters, sample_size = 2500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$EYT,
                   simulation_code(iters = iters, transformed = TRUE, sample_size = 2500,  seeds = seeds,conf = 5, treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1)$EYT),
             type = "latex"))




