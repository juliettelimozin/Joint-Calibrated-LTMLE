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
library(xtable)

set.seed(160625)
seeds <- floor(runif(1000)*10^8)

check_weights <- function(transformed = FALSE, sample_size, conf,seeds){
  time <- proc.time()
  if(transformed){
    print("Transformed covariates")
  } else {print("Correct covariates")}
    set.seed(seeds[1])
    suppressMessages(suppressWarnings({
      simdata<-DATA_GEN(ns = sample_size, nv = 3, conf = conf)
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
      
      calibrate_always_treated <- calibration_by_time_from_baseline(simdata, var = c("X1", "X2", "X3", "X4"), weights_var = 'weight')
      calibrate_always_treated_aggr <- aggregated_calibration_from_baseline(simdata, var = c("X1", "X2", "X3", "X4"), weights_var = 'weight')
      
      simdata$Cweights <- calibrate_always_treated$data$Cweights
      simdata$Cweights_aggr <- calibrate_always_treated_aggr$data$Cweights
      
      simdata$RA <- 1
      simdata[simdata$t == 0 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 1 & !(simdata$CA == 0),]$RA <- 0
      simdata[simdata$t == 2 & !(simdata$CA == 0),]$RA <- 0
      
      calibrate_never_treated <- calibration_by_time_from_baseline(simdata, var = c("X1", "X2", "X3", "X4"), weights_var = 'Cweights')
      calibrate_never_treated_aggr <- aggregated_calibration_from_baseline(simdata, var = c("X1", "X2", "X3", "X4"), weights_var = 'Cweights_aggr')
      
      simdata$Cweights <- calibrate_never_treated$data$Cweights
      simdata$Cweights_aggr <- calibrate_never_treated_aggr$data$Cweights
      
      simdata$weights <- simdata$weight
      
      print(summary(simdata$weights), digits = 2)
      print(summary(simdata$Cweights), digits = 2)
      print(summary(simdata$Cweights_aggr), digits = 2)
    }))
}

check_weights(transformed = FALSE, sample_size = 200, conf = -0.2, seeds = seeds)
check_weights(transformed = FALSE, sample_size = 1000, conf = -0.2, seeds = seeds)
check_weights(transformed = FALSE, sample_size = 2500, conf = -0.2, seeds = seeds)      

check_weights(transformed = T, sample_size = 200, conf = -0.2, seeds = seeds)
check_weights(transformed = T, sample_size = 1000, conf = -0.2, seeds = seeds)
check_weights(transformed = T, sample_size = 2500, conf = -0.2, seeds = seeds)      

check_weights(transformed = FALSE, sample_size = 200, conf = -5, seeds = seeds)
check_weights(transformed = FALSE, sample_size = 1000, conf = -5, seeds = seeds)
check_weights(transformed = FALSE, sample_size = 2500, conf = -5, seeds = seeds)      

check_weights(transformed = T, sample_size = 200, conf = -5, seeds = seeds)
check_weights(transformed = T, sample_size = 1000, conf = -5, seeds = seeds)
check_weights(transformed = T, sample_size = 2500, conf = -5, seeds = seeds)      
