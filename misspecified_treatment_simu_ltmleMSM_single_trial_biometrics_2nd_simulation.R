#!/usr/bin R
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

set.seed(14052025)
seeds <- floor(runif(1000)*10^8)

iters = 1000
registerDoParallel(cores = 10)
sample_size <- 200
nv <- 3
treat_prev <- 1
conf <- 0.5

time <- proc.time()
simulation <- foreach(i = 1:iters, .combine=cbind) %dopar% {
  set.seed(seeds[i])
  suppressMessages(suppressWarnings({
    simdata<-DATA_GEN(ns = sample_size, nv = nv, treat_prev = treat_prev, conf = conf)
    
    # simdata$X1 <- simdata$TX1
    # simdata$X2 <- simdata$TX2
    # simdata$X3 <- simdata$TX3
    # simdata$X4 <- simdata$TX4
    
    simdata$RA <- 1
    simdata[simdata$t == 1 & !(simdata$CA == 2 | simdata$CA == 0),]$RA <- 0
    simdata[simdata$t == 2 & !(simdata$CA == 3 | simdata$CA == 0),]$RA <- 0
    
    
    treatment_model_pooled <- glm(A~Ap+ X1, data = simdata[simdata$t !=0,], family = 'binomial')
    
    treat_model_A_0 <- glm(A~X1, data = simdata[simdata$t == 0,], family = 'binomial')
    treat_model_A_1 <- glm(A~ Ap+ X1, data = simdata[simdata$t == 1,], family = 'binomial')
    treat_model_A_2 <- glm(A~ Ap+ X1 , data = simdata[simdata$t == 2,], family = 'binomial')
    
    
    simdata$ps <- 1.0
    simdata[simdata$t == 0& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 1,], type = 'response'))
    simdata[simdata$t == 1& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 1,], type = 'response'))
    simdata[simdata$t == 2& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 1,], type = 'response'))
    
    simdata[simdata$t == 0& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 0,], type = 'response'))
    simdata[simdata$t == 1& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 0,], type = 'response'))
    simdata[simdata$t == 2& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 0,], type = 'response'))
    
    simdata$weight <- ave(simdata$ps, simdata$ID, FUN = function(X) 1/cumprod(X))
    simdata$weights <- simdata$weight
    simdata$tall <- simdata$t
    
    simdata$RA <- 1
    simdata[simdata$t == 0 & !(simdata$CA == 1),]$RA <- 0
    simdata[simdata$t == 1 & !(simdata$CA == 2),]$RA <- 0
    simdata[simdata$t == 2 & !(simdata$CA == 3),]$RA <- 0
    
    calibrate_always_treated <- calibration_by_time_from_baseline(simdata, var = c("X1"))
    
    simdata <- calibrate_always_treated$data
    simdata$weights <- simdata$Cweights
    
    simdata$RA <- 1
    simdata[simdata$t == 0 & !(simdata$CA == 0),]$RA <- 0
    simdata[simdata$t == 1 & !(simdata$CA == 0),]$RA <- 0
    simdata[simdata$t == 2 & !(simdata$CA == 0),]$RA <- 0
    
    calibrate_never_treated <- calibration_by_time_from_baseline(simdata, var = c("X1"))
    
    simdata <- calibrate_never_treated$data
    
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
    wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, value.var = c("A", "X1", "X2","X3", "X4", "CA", "Y", "weights", "Cweights"))
    
    wideSimdata$g_treat_1_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_1)) %*% treatment_model_pooled$coefficients)
    wideSimdata$g_treat_2_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_2)) %*% treatment_model_pooled$coefficients)
    
    wideSimdata$g_control_1_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_1)) %*% treatment_model_pooled$coefficients)
    wideSimdata$g_control_2_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_2)) %*% treatment_model_pooled$coefficients)
    
    a = min(simdata$Y,na.rm = TRUE); b = max(simdata$Y, na.rm = TRUE)
    
    wideSimdata$Y_2_scaled = (wideSimdata$Y_2-a)/(b-a)
    wideSimdata$Y_1_scaled = (wideSimdata$Y_1-a)/(b-a)
    wideSimdata$Y_0_scaled = (wideSimdata$Y_0-a)/(b-a)
    
    
    
    ########### t = 2 ##########
    
    Q2_2_fit <- glm(data = wideSimdata, formula = Y_2_scaled ~ A_2 + X1_2 + X2_2 + X3_2 +X4_2, family = 'quasibinomial')
    
    logitQ2_2 <- predict.glm(Q2_2_fit, newdata = data.frame(A_2 = c(rep(1,sample_size), rep(0,sample_size)), 
                                                            X1_2 = rep(wideSimdata$X1_2, 2), 
                                                            X2_2 = rep(wideSimdata$X2_2, 2),
                                                            X3_2 = rep(wideSimdata$X3_2, 2), 
                                                            X4_2 = rep(wideSimdata$X4_2, 2)), type = 'link')
    
    #------------- update Q2_2-----------------------
    regimen_ind <- c(as.numeric(wideSimdata$CA_2 ==3), as.numeric(wideSimdata$CA_2 ==0))
    weight <- rep(wideSimdata$weights_2,2)
    Cweight <- rep(wideSimdata$Cweights_2,2)
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Y_2 = rep(wideSimdata$Y_2_scaled, 2), 
                              off = logitQ2_2,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(3,sample_size), rep(0,sample_size)))
    
    Q2_2star_fit <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q2_2star <- predict.glm(Q2_2star_fit, newdata = update_data, type = 'response')
    
    Q2_2star_fit_cali <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(scale(Cweight[regimen_ind == 1], center=FALSE)))
    
    Q2_2star_cali <- predict.glm(Q2_2star_fit_cali, newdata = update_data, type = 'response')
    
    ########### t = 1 ##########
    #------------- Get Q2_1, Q1_1 ---------------------
    fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                               Q2_2star = Q2_2star, 
                               Q2_2star_cali = Q2_2star_cali,
                               CA_1 = rep(wideSimdata$CA_1, 2))
    
    Q2_1_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star ~ CA_1, family = 'quasibinomial')
    Q2_1_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star ~ CA_1, family = 'quasibinomial')
    
    logitQ2_1 <- as.matrix(c(predict.glm(Q2_1_fit_d1, 
                                         newdata = data.frame(CA_1 = rep(2,sample_size)), 
                                         type = 'link'),
                             predict.glm(Q2_1_fit_d0, 
                                         newdata = data.frame(CA_1 = rep(0,sample_size)), 
                                         type = 'link')))
    
    Q2_1_fit_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star_cali ~ CA_1, family = 'quasibinomial')
    Q2_1_fit_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star_cali ~ CA_1, family = 'quasibinomial')
    
    logitQ2_1_cali <- as.matrix(c(predict.glm(Q2_1_fit_d1_cali, 
                                              newdata = data.frame(CA_1 = rep(2,sample_size)), 
                                              type = 'link'),
                                  predict.glm(Q2_1_fit_d0_cali, 
                                              newdata = data.frame(CA_1 = rep(0,sample_size)), 
                                              type = 'link')))
    
    Q1_1_fit <- glm(data = wideSimdata, formula = Y_1_scaled ~ A_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
    
    logitQ1_1 <- predict.glm(Q1_1_fit, newdata = data.frame(A_1 = c(rep(1,sample_size), rep(0,sample_size)), 
                                                            X1_1 = rep(wideSimdata$X1_1, 2), 
                                                            X2_1 = rep(wideSimdata$X2_1, 2),
                                                            X3_1 = rep(wideSimdata$X3_1, 2), 
                                                            X4_1 = rep(wideSimdata$X4_1, 2)), type = 'link')
    
    
    #------------- update Q1_1-----------------------
    regimen_ind = c(as.numeric(wideSimdata$CA_1 ==2), as.numeric(wideSimdata$CA_1 ==0))
    weight <- rep(wideSimdata$weights_1,2)
    Cweight <- rep(wideSimdata$Cweights_1,2)
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Y_1 = rep(wideSimdata$Y_1_scaled, 2), 
                              off = logitQ1_1,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(2,sample_size), rep(0,sample_size)))
    
    Q1_1star_fit <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q1_1star <- predict.glm(Q1_1star_fit, newdata = update_data, type = 'response')
    
    Q1_1star_fit_cali <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(scale(Cweight[regimen_ind == 1], center=FALSE)))
    
    Q1_1star_cali <- predict.glm(Q1_1star_fit_cali, newdata = update_data, type = 'response')
    
    #------------- update Q2_1-----------------------
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Q2_2star = Q2_2star,
                              Q2_2star_cali = Q2_2star_cali,
                              off = logitQ2_1,
                              off_cali = logitQ2_1_cali,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(3,sample_size), rep(0,sample_size)))
    
    Q2_1star_fit <- glm(Q2_2star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q2_1star <- predict.glm(Q2_1star_fit, newdata = update_data, type = 'response')
    
    Q2_1star_fit_cali <- glm(Q2_2star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(scale(Cweight[regimen_ind == 1], center=FALSE)))
    
    Q2_1star_cali <- predict.glm(Q2_1star_fit_cali, newdata = update_data, type = 'response')
    
    ########### t = 0 ##########
    #------------- Get Q4_0, Q3_0, Q2_0, Q1_0, Q0_0 ---------------------
    fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                               Q2_1star = Q2_1star, 
                               Q1_1star = Q1_1star, 
                               Q2_1star_cali = Q2_1star_cali, 
                               Q1_1star_cali = Q1_1star_cali, 
                               CA_0 = rep(wideSimdata$CA_0, 2))
    
    
    Q2_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star ~ 1, family = 'quasibinomial')
    Q2_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star ~ 1, family = 'quasibinomial')
    
    logitQ2_0 <- as.matrix(c(predict.glm(Q2_0_fit_d1, 
                                         type = 'link'),
                             predict.glm(Q2_0_fit_d0, 
                                         type = 'link')))
    
    Q2_0_fit_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star_cali ~ 1, family = 'quasibinomial')
    Q2_0_fit_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star_cali ~ 1, family = 'quasibinomial')
    
    logitQ2_0_cali <- as.matrix(c(predict.glm(Q2_0_fit_d1_cali, 
                                              type = 'link'),
                                  predict.glm(Q2_0_fit_d0_cali, 
                                              type = 'link')))
    
    Q1_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star ~ CA_0, family = 'quasibinomial')
    Q1_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star ~ CA_0, family = 'quasibinomial')
    
    logitQ1_0 <- as.matrix(c(predict.glm(Q1_0_fit_d1, 
                                         newdata = data.frame(CA_0 = rep(1,sample_size)),
                                         type = 'link'),
                             predict.glm(Q1_0_fit_d0, 
                                         newdata = data.frame(CA_0 = rep(0,sample_size)),
                                         type = 'link')))
    
    Q1_0_fit_d1_cali <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star_cali ~ CA_0, family = 'quasibinomial')
    Q1_0_fit_d0_cali <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star_cali ~ CA_0, family = 'quasibinomial')
    
    logitQ1_0_cali <- as.matrix(c(predict.glm(Q1_0_fit_d1_cali, 
                                              newdata = data.frame(CA_0 = rep(1,sample_size)),
                                              type = 'link'),
                                  predict.glm(Q1_0_fit_d0_cali, 
                                              newdata = data.frame(CA_0 = rep(0,sample_size)),
                                              type = 'link')))
    
    Q0_0_fit <- glm(data = wideSimdata, formula = Y_0_scaled ~ CA_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
    
    logitQ0_0 <- predict.glm(Q0_0_fit, newdata = data.frame(CA_0 = c(rep(1,sample_size), rep(0,sample_size)), 
                                                            X1_0 = rep(wideSimdata$X1_0, 2), 
                                                            X2_0 = rep(wideSimdata$X2_0, 2),
                                                            X3_0 = rep(wideSimdata$X3_0, 2), 
                                                            X4_0 = rep(wideSimdata$X4_0, 2)), type = 'link')
    
    
    #------------- update Q0_0-----------------------
    regimen_ind = c(as.numeric(wideSimdata$CA_0 ==1), as.numeric(wideSimdata$CA_0 ==0))
    weight <- rep(wideSimdata$weights_0,2)
    Cweight <- rep(wideSimdata$Cweights_0,2)
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Y_0 = rep(wideSimdata$Y_0_scaled, 2), 
                              off = logitQ0_0,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(1,sample_size), rep(0,sample_size)))
    
    Q0_0star_fit <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q0_0star <- predict.glm(Q0_0star_fit, newdata = update_data, type = 'response')
    
    Q0_0star_fit_cali <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(scale(Cweight[regimen_ind == 1], center=FALSE)))
    
    Q0_0star_cali <- predict.glm(Q0_0star_fit_cali, newdata = update_data, type = 'response')
    
    #------------- update Q2_0-----------------------
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Q2_1star = Q2_1star,
                              Q2_1star_cali = Q2_1star_cali,
                              off = logitQ2_0,
                              off_cali = logitQ2_0_cali,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(3,sample_size), rep(0,sample_size)))
    
    Q2_0star_fit <- glm(Q2_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q2_0star <- predict.glm(Q2_0star_fit, newdata = update_data, type = 'response')
    
    Q2_0star_fit_cali <- glm(Q2_1star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(scale(Cweight[regimen_ind == 1], center=FALSE)))
    
    Q2_0star_cali <- predict.glm(Q2_0star_fit_cali, newdata = update_data, type = 'response')
    
    #------------- update Q1_0-----------------------
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Q1_1star = Q1_1star,
                              Q1_1star_cali = Q1_1star_cali,
                              off = logitQ1_0,
                              off_cali = logitQ1_0_cali,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(2,sample_size), rep(0,sample_size)))
    
    Q1_0star_fit <- glm(Q1_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q1_0star <- predict.glm(Q1_0star_fit, newdata = update_data, type = 'response')
    
    Q1_0star_fit_cali <- glm(Q1_1star_cali ~ intercept + cumA + offset(off_cali) - 1, family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(scale(Cweight[regimen_ind == 1], center=FALSE)))
    
    Q1_0star_cali <- predict.glm(Q1_0star_fit_cali, newdata = update_data, type = 'response')
    
    ################# Fit MSM ##########################
    
    msm_fitting_data <- data.frame(id = rep(1:sample_size,6), 
                                   t = c(rep(0,sample_size*2), 
                                         rep(1,sample_size*2), 
                                         rep(2,sample_size*2)),
                                   Y = c(Q0_0star, Q1_0star, Q2_0star), 
                                   Y_cali = c(Q0_0star_cali, Q1_0star_cali, Q2_0star_cali),
                                   cumA = c(rep(1,sample_size), rep(0,sample_size), 
                                            rep(2,sample_size), rep(0,sample_size), 
                                            rep(3,sample_size), rep(0,sample_size)))
    
    msm <- glm(Y ~ cumA, data = msm_fitting_data, family = 'quasibinomial')
    msm_cali <- glm(Y_cali ~ cumA, data = msm_fitting_data, family = 'quasibinomial')
    
    msm_fitting_data_transformed <- data.frame(id = rep(1:sample_size,6), 
                                               t = c(rep(0,sample_size*2), 
                                                     rep(1,sample_size*2), 
                                                     rep(2,sample_size*2)),
                                               Y = c(Q0_0star, Q1_0star, Q2_0star)*(b-a) + a, 
                                               Y_cali = c(Q0_0star_cali, Q1_0star_cali, Q2_0star_cali)*(b-a) +a,
                                               cumA = c(rep(1,sample_size), rep(0,sample_size), 
                                                        rep(2,sample_size), rep(0,sample_size), 
                                                        rep(3,sample_size), rep(0,sample_size)))
    
    msm_transformed <- glm(Y ~ cumA, data = msm_fitting_data_transformed)
    msm_transformed_cali <- glm(Y_cali ~ cumA, data = msm_fitting_data_transformed)
    
    ################## PACKAGE LTMLE MSM #########################################
    regime_data <- array(1, dim = c(sample_size,3,2))
    regime_data[,,2] <- 0.0*regime_data[,,2]
    summary_measures <- array(, dim = c(2,1,3))
    summary_measures[1,1,] <- c(1,2,3)
    summary_measures[2,1,] <- c(0,0,0)
    colnames(summary_measures) <- 'cumA'
    
    tmle <- ltmleMSM(data = wideSimdata[,.(X1_0, X2_0, X3_0, X4_0,  A_0, Y_0,
                                           X1_1, X2_1, X3_1, X4_1,  A_1, Y_1, 
                                           X1_2, X2_2, X3_2, X4_2,  A_2, Y_2)], 
                     Anodes = c('A_0', 'A_1', 'A_2'), 
                     Lnodes = c('X1_0','X2_0','X3_0','X4_0',
                                'X1_1','X2_1','X3_1','X4_1',
                                'X1_2','X2_2','X3_2','X4_2'), 
                     Ynodes = c('Y_0', 'Y_1', 'Y_2'),
                     Qform = c(X1_0='Q.kplus1 ~ 1', X2_0='Q.kplus1 ~ 1', X3_0='Q.kplus1 ~ 1', X4_0='Q.kplus1 ~ 1',  
                               Y_0='Q.kplus1 ~ X1_0 + X2_0 + X3_0 + X4_0 + A_0',
                               X1_1='Q.kplus1 ~ 1', X2_1='Q.kplus1 ~ 1', X3_1='Q.kplus1 ~ 1', X4_1='Q.kplus1 ~ 1',  
                               Y_1='Q.kplus1 ~ X1_1 + X2_1 + X3_1 + X4_1 + A_0 + A_1',
                               X1_2='Q.kplus1 ~ 1', X2_2='Q.kplus1 ~ 1', X3_2='Q.kplus1 ~ 1', X4_2='Q.kplus1 ~ 1', 
                               Y_2='Q.kplus1 ~ X1_2 + X2_2 + X3_2 + X4_2 + A_0 + A_1 + A_2'),
                     gform = c(A_0 = 'A_0 ~ X1_0 + X2_0 + X3_0 + X4_0', 
                               A_1='A_1 ~ A_0 + X1_1 + X2_1 + X3_1 + X4_1',
                               A_2='A_2 ~ A_1 + X1_2 + X2_2 + X3_2 + X4_2'),
                     gbounds = c(0,1),
                     Yrange = c(a,b),
                     survivalOutcome = FALSE, 
                     regimes = regime_data,
                     working.msm = 'Y ~ cumA',
                     msm.weights = NULL,
                     observation.weights = NULL,
                     summary.measures = summary_measures,
                     final.Ynodes = c('Y_0', 'Y_1', 'Y_2')
    )
    
    ##################### IPW - MSM ##################################
    simdata$RA <- 1
    simdata[simdata$t == 1 & !(simdata$CA == 2 |simdata$CA == 0),]$RA <- 0
    simdata[simdata$t == 2 & !(simdata$CA == 3 |simdata$CA == 0),]$RA <- 0
    
    simdatafinal <- calibration_by_time_from_baseline(as.data.frame(simdata), var = c("X1"))
    
    switch_data <- simdatafinal$data[simdatafinal$data$RA == 1,]
    switch_data$X1 <- ave(switch_data$X1, switch_data$ID, FUN = first)
    switch_data$X2 <- ave(switch_data$X2, switch_data$ID, FUN = first)
    switch_data$X3 <- ave(switch_data$X3, switch_data$ID, FUN = first)
    switch_data$X4 <- ave(switch_data$X4, switch_data$ID, FUN = first)
    
    switch_data$weight <- 1.0
    switch_data[switch_data$t == 1& switch_data$A == 1,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 1,],]$g_treat_1_pooled)
    switch_data[switch_data$t == 2& switch_data$A == 1,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 1,],]$g_treat_1_pooled*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 1,],]$g_treat_2_pooled)
    
    switch_data[switch_data$t == 1& switch_data$A == 0,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 0,],]$g_control_1_pooled)
    switch_data[switch_data$t == 2& switch_data$A == 0,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 0,],]$g_control_1_pooled*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 0,],]$g_control_2_pooled)
    
    PP_pooled <- glm(Y ~ CA + X1, data = switch_data, weights = weight, family = 'gaussian')
    PP_strat <- glm(Y ~ CA, data = switch_data, weights = weights, family = 'gaussian')
    PP_cali <- glm(Y ~ CA, data = switch_data, weights = Cweights, family = 'gaussian')
    
  }))
  
  fitting_data_CA3 <- switch_data[switch_data$t == 0,]
  fitting_data_CA3$CA <- 3
  
  fitting_data_CA0 <- switch_data[switch_data$t == 0,]
  fitting_data_CA0$CA <- 0
  
  c((b-a)*plogis(c(1,3)%*%msm$coefficients) + a - ((b-a)*plogis(c(1,0)%*%msm$coefficients) + a), 
    (b-a)*plogis(c(1,3)%*%msm_cali$coefficients) + a - ((b-a)*plogis(c(1,0)%*%msm_cali$coefficients) + a),
    c(1,3)%*%msm_transformed$coefficients - c(1,0)%*%msm_transformed$coefficients,
    c(1,3)%*%msm_transformed_cali$coefficients - c(1,0)%*%msm_transformed_cali$coefficients,
    (b-a)*plogis(c(1,3)%*%tmle$beta) + a - ((b-a)*plogis(c(1,0)%*%tmle$beta) + a),
    mean(predict.glm(PP_pooled, newdata = fitting_data_CA3, type = 'response')) - mean(predict.glm(PP_pooled, newdata = fitting_data_CA0, type = 'response')),
    predict.glm(PP_strat, newdata = data.frame(CA = 3), type = 'response') - predict.glm(PP_strat, newdata = data.frame(CA = 0), type = 'response'),
    predict.glm(PP_cali, newdata = data.frame(CA = 3), type = 'response') - predict.glm(PP_cali, newdata = data.frame(CA = 0), type = 'response'))
}

cat('Misspecified treat, correctly specified outcome')
cat(paste('Sample size = ', sample_size))
cat(paste('Treat prev = ', treat_prev))
cat(paste('Confounding = ', conf))

cat(paste('Estimation of ATE_2 = E(Y_2(always treated )) - E(Y_2(never treated))\n'))

cat(paste('Time taken for iters = ', iters))
print(proc.time() - time)

cat('\n')
rownames(simulation) <- c('Manual LTMLE-MSM with MSM fitted on Q*_0', 
                          'Manual calibrated LTMLE-MSM with MSM fitted on Q*_0',
                          'Manual LTMLE-MSM with MSM fitted on transformed Q*_0',
                          'Manual calibrated LTMLE-MSM with MSM fitted on transformed Q*_0',
                          'Package ltmleMSM', 
                          'IPW-MSM, pooled treatment model', 
                          'IPW-MSM, stratified treatment model',
                          'Calibrated IPW-MSM, stratified treatment model')
cat('Bias: \n')
print(rowMeans(simulation, na.rm = TRUE)-15)

cat('SD: \n')
print(rowSds(simulation, na.rm = TRUE))

cat('rootMSE: \n')
print(sqrt((rowMeans(simulation, na.rm = TRUE)-15)^2+ rowSds(simulation, na.rm = TRUE)^2))
