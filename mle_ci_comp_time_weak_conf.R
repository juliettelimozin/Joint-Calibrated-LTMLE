#!/usr/bin R
###############################
# This script runs simulations with Q functions on logit form
###############################
setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation/dgm_2nd_simulation_biometrics.R")
source("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation/calibration_func_trials.R")
source("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation/comp_time_ltmle_utils.R")
library(dplyr)
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
conf = 0.2
treat_prev_0 = 0
treat_prev_d1_1 = 1
treat_prev_d0_1 = -1.25
treat_prev_d1_2 =0.8
treat_prev_d0_2 = -1.25

l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
sample_size = l
transformed = FALSE

set.seed(160625)
seeds <- floor(runif(1000)*10^8)


iters = 1000
bootstrap_iter = 500
registerDoParallel(cores = 67)

combine_matrices <- function(x, y) {
  list(
    computation_time = cbind(x$computation_time, y$computation_time))
}

#sink("job_output.txt", append = TRUE, split = TRUE)

simulation <- foreach(i = 1:iters, .combine=combine_matrices, .init =  list(
  computation_time = NULL)) %dopar% {
    set.seed(seeds[i])
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
    
    point_estimate_TMLE <-  MyLtmleMSM(simdata)
    point_estimate_modifiedTMLE <- MyLtmleMSM_modified(simdata, initial_Q_t = NA, refit_weights = TRUE)

    simdata_with_weights <- point_estimate_modifiedTMLE$simdata
    ########## Bootstrap ############
    
    
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
    time <- proc.time()
    
    modified_bootstrap_LTMLE <- foreach(k = 1:bootstrap_iter, .combine = cbind) %do% {
      
      boot_design_data <- boot_samples[[k]]
      
      bootstrap_tmle <- MyLtmleMSM_modified(simdata = boot_design_data,initial_Q_t = point_estimate_modifiedTMLE$Q_t, refit_weights = FALSE)
      
      bootstrap_tmle$estimates
    }
    
    computation_time <- (proc.time()-time)[[3]]
    
    ############################# NORMAL BOOTSTRAP TMLE #############################
    time <- proc.time()
    normal_bootstrap_LTMLE <- foreach(k = 1:bootstrap_iter, .combine=cbind) %do% {
      
      boot_design_data <- as.data.frame(boot_samples[[k]])[, colnames(simdata)]
      boot_design_data$ID <- as.integer(factor(boot_design_data$ID))
      bootstrap_tmle <- MyLtmleMSM(simdata = boot_design_data)
      
      bootstrap_tmle$estimates
    }
    
    computation_time <- rbind(computation_time,(proc.time()-time)[[3]])
    message("Finished iteration ", i)
    flush.console()
    list(computation_time = computation_time)
  }
#sink()
saveRDS(simulation, file = paste0('mle_ci_comp_time_weak_conf_result_', as.character(l),'.rds'))
