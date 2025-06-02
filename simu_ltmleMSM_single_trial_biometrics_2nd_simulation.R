#!/usr/bin R
library(dplyr)
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/dgm_2nd_simulation_biometrics.R")
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

iters = 500
registerDoParallel(cores = 10)
sample_size <- 100

time <- proc.time()
simulation <- foreach(i = 1:iters, .combine=cbind) %dopar% {
  set.seed(seeds[i])
  suppressMessages(suppressWarnings({
    simdata<-DATA_GEN(ns = sample_size, treat_prev = 0, conf = 0.2)

    treatment_model_pooled <- glm(A~Ap+ X1+ X2+ X3 + X4, data = simdata[simdata$t !=0,], family = 'quasibinomial')
    
    treat_model_A_0 <- glm(A~X1 + X2 + X3 + X4, data = simdata[simdata$t == 0,], family = 'quasibinomial')
    treat_model_A_1 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 1,], family = 'quasibinomial')
    treat_model_A_2 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 2,], family = 'quasibinomial')
    treat_model_A_3 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 3,], family = 'quasibinomial')
    treat_model_A_4 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 4,], family = 'quasibinomial')
    
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
    wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, value.var = c("A", "X1", "X2","X3", "X4", "CA", "Y"))
   
    regime_data <- array(1, dim = c(sample_size,5,2))
    regime_data[,,2] <- 0.0*regime_data[,,2]
    summary_measures <- array(, dim = c(2,1,5))
    summary_measures[1,1,] <- c(1,2,3,4,5)
    summary_measures[2,1,] <- c(0,0,0,0,0)
    colnames(summary_measures) <- 'cumA'
    
    
    wideSimdata$g_treat_1_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_1, wideSimdata$X2_1,wideSimdata$X3_1, wideSimdata$X4_1)) %*% treatment_model_pooled$coefficients)
    wideSimdata$g_treat_2_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_2, wideSimdata$X2_2,wideSimdata$X3_2, wideSimdata$X4_2)) %*% treatment_model_pooled$coefficients)
    wideSimdata$g_treat_3_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_3, wideSimdata$X2_3,wideSimdata$X3_3, wideSimdata$X4_3)) %*% treatment_model_pooled$coefficients)
    wideSimdata$g_treat_4_pooled <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_4, wideSimdata$X2_4,wideSimdata$X3_4, wideSimdata$X4_4)) %*% treatment_model_pooled$coefficients)
    
    wideSimdata$g_control_1_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_1, wideSimdata$X2_1,wideSimdata$X3_1, wideSimdata$X4_1)) %*% treatment_model_pooled$coefficients)
    wideSimdata$g_control_2_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_2, wideSimdata$X2_2,wideSimdata$X3_2, wideSimdata$X4_2)) %*% treatment_model_pooled$coefficients)
    wideSimdata$g_control_3_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_3, wideSimdata$X2_3,wideSimdata$X3_3, wideSimdata$X4_3)) %*% treatment_model_pooled$coefficients)
    wideSimdata$g_control_4_pooled <- 1-plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_4, wideSimdata$X2_4,wideSimdata$X3_4, wideSimdata$X4_4)) %*% treatment_model_pooled$coefficients)
    
    wideSimdata$g_treat_0 <- plogis(as.matrix(cbind(rep(1,sample_size), wideSimdata$X1_0, wideSimdata$X2_0,wideSimdata$X3_0, wideSimdata$X4_0)) %*% treat_model_A_0$coefficients)
    wideSimdata$g_treat_1 <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_1, wideSimdata$X2_1,wideSimdata$X3_1, wideSimdata$X4_1)) %*% treat_model_A_1$coefficients)
    wideSimdata$g_treat_2 <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_2, wideSimdata$X2_2,wideSimdata$X3_2, wideSimdata$X4_2)) %*% treat_model_A_2$coefficients)
    wideSimdata$g_treat_3 <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_3, wideSimdata$X2_3,wideSimdata$X3_3, wideSimdata$X4_3)) %*% treat_model_A_3$coefficients)
    wideSimdata$g_treat_4 <- plogis(as.matrix(cbind(rep(1,sample_size), rep(1,sample_size), wideSimdata$X1_4, wideSimdata$X2_4,wideSimdata$X3_4, wideSimdata$X4_4)) %*% treat_model_A_4$coefficients)
    
    wideSimdata$g_control_0 <- 1 - plogis(as.matrix(cbind(rep(1,sample_size), wideSimdata$X1_0, wideSimdata$X2_0,wideSimdata$X3_0, wideSimdata$X4_0)) %*% treat_model_A_0$coefficients)
    wideSimdata$g_control_1 <- 1 - plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_1, wideSimdata$X2_1,wideSimdata$X3_1, wideSimdata$X4_1)) %*% treat_model_A_1$coefficients)
    wideSimdata$g_control_2 <- 1 - plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_2, wideSimdata$X2_2,wideSimdata$X3_2, wideSimdata$X4_2)) %*% treat_model_A_2$coefficients)
    wideSimdata$g_control_3 <- 1 - plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_3, wideSimdata$X2_3,wideSimdata$X3_3, wideSimdata$X4_3)) %*% treat_model_A_3$coefficients)
    wideSimdata$g_control_4 <- 1 - plogis(as.matrix(cbind(rep(1,sample_size), rep(0,sample_size), wideSimdata$X1_4, wideSimdata$X2_4,wideSimdata$X3_4, wideSimdata$X4_4)) %*% treat_model_A_4$coefficients)
    
    a = min(simdata$Y,na.rm = TRUE); b = max(simdata$Y, na.rm = TRUE)
    
    wideSimdata$Y_4_scaled = (wideSimdata$Y_4-a)/(b-a)
    wideSimdata$Y_3_scaled = (wideSimdata$Y_3-a)/(b-a)
    wideSimdata$Y_2_scaled = (wideSimdata$Y_2-a)/(b-a)
    wideSimdata$Y_1_scaled = (wideSimdata$Y_1-a)/(b-a)
    wideSimdata$Y_0_scaled = (wideSimdata$Y_0-a)/(b-a)
    
    
    ########### t = 4 ##########
    Q4_4_fit <- glm(data = wideSimdata, formula = Y_4_scaled ~ A_4 + X1_4 + X2_4 + X3_4 + X4_4, family = 'quasibinomial')
    
    logitQ4_4 <- predict.glm(Q4_4_fit, newdata = data.frame(A_4 = c(rep(1,sample_size), rep(0,sample_size)), X1_4 = rep(wideSimdata$X1_4, 2), X2_4 = rep(wideSimdata$X2_4, 2), X3_4 = rep(wideSimdata$X3_4, 2), X4_4 = rep(wideSimdata$X4_4, 2)), type = 'link')
    
    #------------- update Q2_2-----------------------
    regimen_ind = c(as.numeric(wideSimdata$CA_4 ==5), as.numeric(wideSimdata$CA_4 ==0))
    weight = 1/c(wideSimdata$g_treat_4*wideSimdata$g_treat_3*wideSimdata$g_treat_2*wideSimdata$g_treat_1*wideSimdata$g_treat_0,  wideSimdata$g_control_4*wideSimdata$g_control_3*wideSimdata$g_control_2*wideSimdata$g_control_1*wideSimdata$g_control_0)
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Y_4 = rep(wideSimdata$Y_4_scaled, 2), 
                              off = logitQ4_4,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(5,sample_size), rep(0,sample_size)))
    
    Q4_4star_fit <- glm(Y_4 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q4_4star <- predict.glm(Q4_4star_fit, newdata = update_data, type = 'response')
    
    ########### t = 3 ##########
    Q3_3_fit <- glm(data = wideSimdata, formula = Y_3_scaled ~ A_3 + X1_3 + X2_3 + X3_3 + X4_3, family = 'quasibinomial')
    
    logitQ3_3 <- predict.glm(Q3_3_fit, newdata = data.frame(A_3 = c(rep(1,sample_size), rep(0,sample_size)), X1_3 = rep(wideSimdata$X1_3, 2), X2_3 = rep(wideSimdata$X2_3, 2), X3_3 = rep(wideSimdata$X3_3, 2), X4_3 = rep(wideSimdata$X4_3, 2)), type = 'link')
    
    fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), Q4_4star = Q4_4star, CA_3 = rep(wideSimdata$CA_3,2))
    
    Q4_3_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q4_4star ~ CA_3, family = 'quasibinomial')
    Q4_3_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q4_4star ~ CA_3, family = 'quasibinomial')
    
    logitQ4_3 <- as.matrix(c(predict.glm(Q4_3_fit_d1, 
                                         newdata = data.frame(CA_3 = rep(4,sample_size)), 
                                         type = 'link'),
                             predict.glm(Q4_3_fit_d0, 
                                         newdata = data.frame(CA_3 = rep(0,sample_size)), 
                                         type = 'link')))
    
    #------------- update Q3_3-----------------------
    regimen_ind = c(as.numeric(wideSimdata$CA_3 == 4), as.numeric(wideSimdata$CA_3 ==0))
    weight = 1/c(wideSimdata$g_treat_3*wideSimdata$g_treat_2*wideSimdata$g_treat_1*wideSimdata$g_treat_0, wideSimdata$g_control_3*wideSimdata$g_control_2*wideSimdata$g_control_1*wideSimdata$g_control_0)
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Y_3 = rep(wideSimdata$Y_3_scaled, 2), 
                              off = logitQ3_3,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(4,sample_size), rep(0,sample_size)))
    
    Q3_3star_fit <- glm(Y_3 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q3_3star <- predict.glm(Q3_3star_fit, newdata = update_data, type = 'response')
    
    #------------- update Q4_3-----------------------
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Q4_4star = Q4_4star,
                              off = logitQ4_3,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(5,sample_size), rep(0,sample_size)))
    
    Q4_3star_fit <- glm(Q4_4star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q4_3star <- predict.glm(Q4_3star_fit, newdata = update_data, type = 'response')
    
    ########### t = 2 ##########
    
    Q2_2_fit <- glm(data = wideSimdata, formula = Y_2_scaled ~ A_2 + cumX1_2 + cumX2_2, family = 'quasibinomial')
    
    logitQ2_2 <- predict.glm(Q2_2_fit, newdata = data.frame(A_2 = c(rep(1,sample_size), rep(0,sample_size)), cumX1_2 = rep(wideSimdata$cumX1_2, 2), cumX2_2 = rep(wideSimdata$cumX2_2, 2)), type = 'link')
    
    fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), Q4_3star = Q4_3star, Q3_3star = Q3_3star, CA_2 = rep(wideSimdata$CA_2,2))
    
    Q4_2_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q4_3star ~ 1, family = 'quasibinomial')
    Q4_2_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q4_3star ~ 1, family = 'quasibinomial')
    
    logitQ4_2 <- as.matrix(c(predict.glm(Q4_2_fit_d1, 
                                         #newdata = data.frame(CA_3 = rep(4,sample_size)), 
                                         type = 'link'),
                             predict.glm(Q4_2_fit_d0, 
                                         newdata = data.frame(CA_3 = rep(0,sample_size)), 
                                         type = 'link')))
    
    Q3_2_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q3_3star ~ CA_2, family = 'quasibinomial')
    Q3_2_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q3_3star ~ CA_2, family = 'quasibinomial')
    
    logitQ3_2 <- as.matrix(c(predict.glm(Q3_2_fit_d1, 
                                         newdata = data.frame(CA_2 = rep(3,sample_size)), 
                                         type = 'link'),
                             predict.glm(Q3_2_fit_d0, 
                                         newdata = data.frame(CA_2 = rep(0,sample_size)), 
                                         type = 'link')))
    
    
    #------------- update Q2_2-----------------------
    regimen_ind = c(as.numeric(wideSimdata$CA_2 ==3), as.numeric(wideSimdata$CA_2 ==0))
    weight = 1/c(wideSimdata$g_treat_2*wideSimdata$g_treat_1*wideSimdata$g_treat_0, wideSimdata$g_control_2*wideSimdata$g_control_1*wideSimdata$g_control_0)
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Y_2 = rep(wideSimdata$Y_2_scaled, 2), 
                              off = logitQ2_2,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(3,sample_size), rep(0,sample_size)))
    
    Q2_2star_fit <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q2_2star <- predict.glm(Q2_2star_fit, newdata = update_data, type = 'response')
    
    ########### t = 1 ##########
    #------------- Get Q2_1, Q1_1 ---------------------
    fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), Q2_2star = Q2_2star, A_1 = rep(wideSimdata$A_1, 2), cumX1_1 = rep(wideSimdata$cumX1_1,2), cumX2_1 = rep(wideSimdata$cumX2_1,2))
    
    Q2_1_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star ~ A_1 + cumX1_1 + cumX2_1, family = 'quasibinomial')
    Q2_1_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star ~ A_1 + cumX1_1 + cumX2_1, family = 'quasibinomial')
    
    logitQ2_1 <- as.matrix(c(predict.glm(Q2_1_fit_d1, 
                                         newdata = data.frame(A_1 = rep(1,sample_size),cumX1_1 = fitting_data[1:sample_size,]$cumX1_1, cumX2_1 = fitting_data[1:sample_size,]$cumX2_1), 
                                         type = 'link'),
                             predict.glm(Q2_1_fit_d0, 
                                         newdata = data.frame(A_1 = rep(0,sample_size),cumX1_1 = fitting_data[(sample_size+1):(sample_size*2),]$cumX1_1, cumX2_1 = fitting_data[(sample_size+1):(sample_size*2),]$cumX2_1), 
                                         type = 'link')))
    
    Q1_1_fit <- glm(data = wideSimdata, formula = Y_1_scaled ~ A_1 + cumX1_1 + cumX2_1, family = 'quasibinomial')
    
    logitQ1_1 <- predict.glm(Q1_1_fit, newdata = data.frame(A_1 = c(rep(1,sample_size), rep(0,sample_size)), cumX1_1 = rep(wideSimdata$cumX1_1, 2), cumX2_1 = rep(wideSimdata$cumX2_1, 2)), type = 'link')
    
    
    #------------- update Q1_1-----------------------
    regimen_ind = c(as.numeric(wideSimdata$CA_1 ==2), as.numeric(wideSimdata$CA_1 ==0))
    weight = 1/c(wideSimdata$g_treat_1*wideSimdata$g_treat_0, wideSimdata$g_control_1*wideSimdata$g_control_0)
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Y_1 = rep(wideSimdata$Y_1_scaled, 2), 
                              off = logitQ1_1,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(2,sample_size), rep(0,sample_size)))
    
    Q1_1star_fit <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q1_1star <- predict.glm(Q1_1star_fit, newdata = update_data, type = 'response')
    
    #------------- update Q2_1-----------------------
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Q2_2star = Q2_2star,
                              off = logitQ2_1,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(3,sample_size), rep(0,sample_size)))
    
    Q2_1star_fit <- glm(Q2_2star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q2_1star <- predict.glm(Q2_1star_fit, newdata = update_data, type = 'response')
    
    ########### t = 0 ##########
    #------------- Get Q2_0, Q1_0, Q0_0 ---------------------
    fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), Q2_1star = Q2_1star, Q1_1star = Q1_1star, A_0 = rep(wideSimdata$A_0, 2), X1_0 = rep(wideSimdata$X1_0,2), X2_0 = rep(wideSimdata$X2_0,2))
    
    Q2_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star ~ A_0 + X1_0 + X2_0, family = 'quasibinomial')
    Q2_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star ~ A_0 + X1_0 + X2_0, family = 'quasibinomial')
    
    logitQ2_0 <- as.matrix(c(predict.glm(Q2_0_fit_d1, 
                                         newdata = data.frame(A_0 = rep(1,sample_size),X1_0 = fitting_data[1:sample_size,]$X1_0, X2_0 = fitting_data[1:sample_size,]$X2_0), 
                                         type = 'link'),
                             predict.glm(Q2_0_fit_d0, 
                                         newdata = data.frame(A_0 = rep(0,sample_size),X1_0 = fitting_data[(sample_size+1):(sample_size*2),]$X1_0, X2_0 = fitting_data[(sample_size+1):(sample_size*2),]$X2_0), 
                                         type = 'link')))
    
    Q1_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star ~ A_0 + X1_0 + X2_0, family = 'quasibinomial')
    Q1_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star ~ A_0 + X1_0 + X2_0, family = 'quasibinomial')
    
    logitQ1_0 <- as.matrix(c(predict.glm(Q1_0_fit_d1, 
                                         newdata = data.frame(A_0 = rep(1,sample_size),X1_0 = fitting_data[1:sample_size,]$X1_0, X2_0 = fitting_data[1:sample_size,]$X2_0), 
                                         type = 'link'),
                             predict.glm(Q1_0_fit_d0, 
                                         newdata = data.frame(A_0 = rep(0,sample_size),X1_0 = fitting_data[(sample_size+1):(sample_size*2),]$X1_0, X2_0 = fitting_data[(sample_size+1):(sample_size*2),]$X2_0), 
                                         type = 'link')))
    
    Q0_0_fit <- glm(data = wideSimdata, formula = Y_0_scaled ~ A_0 + X1_0 + X2_0, family = 'quasibinomial')
    
    logitQ0_0 <- predict.glm(Q0_0_fit, newdata = data.frame(A_0 = c(rep(1,sample_size), rep(0,sample_size)), X1_0 = rep(wideSimdata$X1_0, 2), X2_0 = rep(wideSimdata$X2_0, 2)), type = 'link')
    
    
    
    #------------- update Q0_0-----------------------
    regimen_ind = c(as.numeric(wideSimdata$A_0 ==1), as.numeric(wideSimdata$A_0 ==0))
    weight = 1/c(wideSimdata$g_treat_0, wideSimdata$g_control_0)
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Y_0 = rep(wideSimdata$Y_0_scaled, 2), 
                              off = logitQ0_0,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(1,sample_size), rep(0,sample_size)))
    
    Q0_0star_fit <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q0_0star <- predict.glm(Q0_0star_fit, newdata = update_data, type = 'response')
    #------------- update Q2_0-----------------------
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Q2_1star = Q2_1star,
                              off = logitQ2_0,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(3,sample_size), rep(0,sample_size)))
    
    Q2_0star_fit <- glm(Q2_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q2_0star <- predict.glm(Q2_0star_fit, newdata = update_data, type = 'response')
    
    #------------- update Q1_0-----------------------
    
    update_data <- data.frame(id = rep(wideSimdata$ID,2),
                              Q1_1star = Q1_1star,
                              off = logitQ1_0,
                              intercept = rep(1,sample_size*2), 
                              cumA = c(rep(2,sample_size), rep(0,sample_size)))
    
    Q1_0star_fit <- glm(Q1_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                        data = update_data[regimen_ind == 1,],
                        weights = as.vector(scale(weight[regimen_ind == 1], center=FALSE)))
    
    Q1_0star <- predict.glm(Q1_0star_fit, newdata = update_data, type = 'response')
    
    ################# Fit MSM ##########################
    
    msm_fitting_data <- data.frame(id = rep(1:sample_size,6), t = c(rep(0,sample_size*2), rep(1,sample_size*2), rep(2,sample_size*2)),Y = c(Q0_0star, Q1_0star, Q2_0star), cumA = c(rep(1,sample_size), rep(0,sample_size), rep(2,sample_size), rep(0,sample_size), rep(3,sample_size), rep(0,sample_size)))
    msm <- glm(Y ~ cumA, data = msm_fitting_data, family = 'quasibinomial')
    #print(msm$coefficients)
    #print(paste('E(Y_2(1)) =', (b-a)*plogis(c(1,3)%*%msm$coefficients) + a))
    #print(paste('E(Y_2(0)) =', (b-a)*plogis(c(1,0)%*%msm$coefficients) + a))
    #print(paste('ATE(t = 2) =', (b-a)*plogis(c(1,3)%*%msm$coefficients) + a - ((b-a)*plogis(c(1,0)%*%msm$coefficients) + a)))
    
    msm_fitting_data_transformed <- data.frame(id = rep(1:sample_size,6), t = c(rep(0,sample_size*2), rep(1,sample_size*2), rep(2,sample_size*2)),Y = c(Q0_0star, Q1_0star, Q2_0star)*(b-a)+a, cumA = c(rep(1,sample_size), rep(0,sample_size), rep(2,sample_size), rep(0,sample_size), rep(3,sample_size), rep(0,sample_size)))
    msm_transformed <- glm(Y ~ cumA, data = msm_fitting_data_transformed)
    #print(msm_transformed$coefficients)
    
    #print(paste('SCALED MSM E(Y_2(1)) =', c(1,3)%*%msm_transformed$coefficients))
    #print(paste('SCALED MSM E(Y_2(0)) =', c(1,0)%*%msm_transformed$coefficients))
    #print(paste('SCALED MSM ATE(t = 2) =', c(1,3)%*%msm_transformed$coefficients - c(1,0)%*%msm_transformed$coefficients))
    
    ################## PACKAGE LTMLE MSM #########################################
    regime_data <- array(1, dim = c(sample_size,3,2))
    regime_data[,,2] <- 0.0*regime_data[,,2]
    summary_measures <- array(, dim = c(2,1,3))
    summary_measures[1,1,] <- c(1,2,3)
    summary_measures[2,1,] <- c(0,0,0)
    colnames(summary_measures) <- 'cumA'
    
    tmle <- ltmleMSM(data = wideSimdata[,.(X1_0,X2_0,A_0,Y_0,X1_1,X2_1,cumX1_1, cumX2_1, A_1, Y_1, X1_2, X2_2, cumX1_2, cumX2_2,A_2, Y_2)], 
                     Anodes = c('A_0', 'A_1', 'A_2'), 
                     Lnodes = c('X1_0','X2_0', 'X1_1', 'X2_1', 'cumX1_1', 'cumX2_1', 'X1_2', 'X2_2','cumX1_2', 'cumX2_2'), 
                     Ynodes = c('Y_0', 'Y_1', 'Y_2'),
                     Qform = c(X1_0 = 'Q.kplus1 ~ 1',
                               X2_0 = 'Q.kplus1 ~ 1',
                               Y_0 = 'Q.kplus1 ~ X1_0 + X2_0 + A_0',
                               X1_1 = 'Q.kplus1 ~ 1',
                               X2_1 = 'Q.kplus1 ~ 1',
                               cumX1_1 = 'Q.kplus1 ~ 1',
                               cumX2_1 = 'Q.kplus1 ~ 1',
                               Y_1 = 'Q.kplus1 ~ cumX1_1 + cumX2_1 + A_1',
                               X1_2 = 'Q.kplus1 ~ 1',
                               X2_2 = 'Q.kplus1 ~ 1',
                               cumX1_2 = 'Q.kplus1 ~ 1',
                               cumX2_2 = 'Q.kplus1 ~ 1',
                               Y_2 = 'Q.kplus1 ~ cumX1_2 + cumX2_2 + A_2'),
                     gform = c(A_0 = 'A_0 ~ 1', A_1 = 'A_1 ~ A_0 + X1_1 + X2_1', A_2 = 'A_2 ~ A_1 + X1_2 + X2_2'),
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
    #print(paste('PACKAGE E(Y_2(1)) =', (b-a)*plogis(c(1,3)%*%tmle$beta) + a))
    #print(paste('PACKAGE E(Y_2(0)) =', (b-a)*plogis(c(1,0)%*%tmle$beta) + a))
    #print(paste('PACKAGE ATE(t = 2) =', (b-a)*plogis(c(1,3)%*%tmle$beta) + a - ((b-a)*plogis(c(1,0)%*%tmle$beta) + a)))
    
    ##################### IPW - MSM ##################################
    switch_data <- rbind(simdata[simdata$t == 0,],
                         simdata[simdata$t == 1 & (simdata$CA ==2 | simdata$CA ==0),],
                         simdata[simdata$t == 2 & (simdata$CA ==3 | simdata$CA ==0),])
    switch_data$weight <- 1.0
    switch_data[switch_data$t == 1& switch_data$A == 1,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 1,],]$g_treat_1_pooled)
    switch_data[switch_data$t == 2& switch_data$A == 1,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 1,],]$g_treat_1_pooled*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 1,],]$g_treat_2_pooled)
    
    switch_data[switch_data$t == 1& switch_data$A == 0,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 0,],]$g_control_1_pooled)
    switch_data[switch_data$t == 2& switch_data$A == 0,]$weight <- 1/(wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 0,],]$g_control_1_pooled*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 0,],]$g_control_2_pooled)
    
    switch_data$weight_strat <- 1.0
    switch_data[switch_data$t == 0& switch_data$A == 1,]$weight_strat <- 1/(wideSimdata[switch_data[switch_data$t == 0& switch_data$A == 1,],]$g_treat_0)
    switch_data[switch_data$t == 1& switch_data$A == 1,]$weight_strat <- 1/(wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 1,],]$g_treat_0*wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 1,],]$g_treat_1)
    switch_data[switch_data$t == 2& switch_data$A == 1,]$weight_strat <- 1/(wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 1,],]$g_treat_0*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 1,],]$g_treat_1*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 1,],]$g_treat_2)
    
    switch_data[switch_data$t == 0& switch_data$A == 0,]$weight_strat <- 1/(wideSimdata[switch_data[switch_data$t == 0& switch_data$A == 0,],]$g_control_0)
    switch_data[switch_data$t == 1& switch_data$A == 0,]$weight_strat <- 1/(wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 0,],]$g_control_0*wideSimdata[switch_data[switch_data$t == 1& switch_data$A == 0,],]$g_control_1)
    switch_data[switch_data$t == 2& switch_data$A == 0,]$weight_strat <- 1/(wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 0,],]$g_control_0*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 0,],]$g_control_1*wideSimdata[switch_data[switch_data$t == 2& switch_data$A == 0,],]$g_control_2)
    
    PP_pooled <- glm(Y ~ CA, data = switch_data, weights = weight, family = 'gaussian')
    PP_strat <- glm(Y ~ CA, data = switch_data, weights = weight_strat, family = 'gaussian')
  }))
  
  c((b-a)*plogis(c(1,3)%*%msm$coefficients) + a - ((b-a)*plogis(c(1,0)%*%msm$coefficients) + a), 
    c(1,3)%*%msm_transformed$coefficients - c(1,0)%*%msm_transformed$coefficients,
    (b-a)*plogis(c(1,3)%*%tmle$beta) + a - ((b-a)*plogis(c(1,0)%*%tmle$beta) + a),
    predict.glm(PP_pooled, newdata = data.frame(CA = 3), type = 'response') - predict.glm(PP_pooled, newdata = data.frame(CA = 0), type = 'response'),
    predict.glm(PP_strat, newdata = data.frame(CA = 3), type = 'response') - predict.glm(PP_strat, newdata = data.frame(CA = 0), type = 'response'))
}

cat(paste('Estimation of ATE_2 = E(Y_2(always treated )) - E(Y_2(never treated))\n'))

cat(paste('Time taken for iters = ', iters))
print(proc.time() - time)

cat('\n')
rownames(simulation) <- c('Manual LTMLE-MSM with MSM fitted on Q*_0', 'Manual LTMLE-MSM with MSM fitted on transformed Q*_0', 'Package ltmleMSM', 'IPW-MSM, pooled treatment model', 'IPW-MSM, stratified treatment model')
cat('Bias: \n')
print(rowMeans(simulation, na.rm = TRUE)+9)

cat('SD: \n')
print(rowSds(simulation, na.rm = TRUE))

cat('rootMSE: \n')
print(sqrt((rowMeans(simulation, na.rm = TRUE)+9)^2+ rowSds(simulation, na.rm = TRUE)^2))
