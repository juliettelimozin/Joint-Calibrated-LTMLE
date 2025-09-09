#!/usr/bin R
###############################
# This script runs simulations with Q functions on logit form
###############################
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/dgm_biometrics_censoring.R")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/calibration_func_trials.R")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/ltmle_utils.R")
library(dplyr)
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
conf = 0.2
treat_prev_0 = 0
treat_prev_d1_1 = 1
treat_prev_d0_1 = -1.25
treat_prev_d1_2 =0.8
treat_prev_d0_2 = -1.25
sample_size = 300
transformed = FALSE

set.seed(160625)
seeds <- floor(runif(1000)*10^8)

iters = 1000
registerDoParallel(cores = 10)

combine_matrices <- function(x, y) {
  list(
    mat1 = cbind(x$modified_bootstrap_CIs, y$modified_bootstrap_CIs),
    mat2 = cbind(x$normal_bootstrap_CIs, y$normal_bootstrap_CIs)
  )
}

simulation <- foreach(i = 1:iters, .combine=combine_matrices) %dopar% {
  set.seed(seeds[i])
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
  
  point_estimate_TMLE <-  MyLtmleMSM(simdata)
  point_estimate_modifiedTMLE <- MyLtmleMSM_modified(simdata, initial_Q_t = NA, refit_weights = TRUE)
  
  simdata_with_weights <- point_estimate_modifiedTMLE$simdata
  ########## Bootstrap ############
  bootstrap_iter = 500
  boot_samples <- vector("list", bootstrap_iter)
  for (b in 1:bootstrap_iter) {
    boot_ids <- sort(sample(unique(simdata_with_weights$ID), length(unique(simdata_with_weights$ID)), replace = TRUE))
    boot_samples[[b]] <- do.call(rbind, lapply(seq_along(boot_ids), function(i) {
      tmp <- simdata_with_weights[simdata_with_weights$ID == boot_ids[i], ]
      tmp$bootstrap_rep <- b
      tmp$ID <- paste0(tmp$ID, "_boot", i)  # assign unique ID
      tmp
    }))
  }
  
  
  ############################# MODIFIED BOOTSTRAP TMLE #############################
  
  
  
  modified_bootstrap_LTMLE <-as.data.frame(matrix(,length(point_estimate_modifiedTMLE$estimates),bootstrap_iter))
  
  modified_bootstrap_LTMLE <- foreach(k = 1:bootstrap_iter, .combine=cbind) %do% {
    
    boot_design_data <- boot_samples[[k]]
    
    bootstrap_tmle <- MyLtmleMSM_modified(simdata = boot_design_data,initial_Q_t = point_estimate_modifiedTMLE$Q_t, refit_weights = FALSE)
    
    bootstrap_tmle$estimates
  }
  
  modified_bootstrap_CIs <- cbind(rowQuantiles(modified_bootstrap_LTMLE, probs = 0.025), rowQuantiles(modified_bootstrap_LTMLE, probs = 0.975))
  
  ############################# NORMAL BOOTSTRAP TMLE #############################
  
  normal_bootstrap_LTMLE <-as.data.frame(matrix(,length(point_estimate_TMLE$estimates),bootstrap_iter))
  normal_bootstrap_LTMLE <- foreach(k = 1:bootstrap_iter, .combine=cbind) %do% {
    
    boot_design_data <- as.data.frame(boot_samples[[k]])[, colnames(simdata)]
    boot_design_data$ID <- as.integer(factor(boot_design_data$ID))
    bootstrap_tmle <- MyLtmleMSM(simdata = boot_design_data)
    
    bootstrap_tmle$estimates
  }
  
  normal_bootstrap_CIs <- cbind(rowQuantiles(normal_bootstrap_LTMLE, probs = 0.025), rowQuantiles(normal_bootstrap_LTMLE, probs = 0.975))
  
  list(modified_bootstrap_CIs = modified_bootstrap_CIs, normal_bootstrap_CIs = normal_bootstrap_CIs)
}

saveRDS(simulation, file = 'ci_coverage_simu_result.rds')