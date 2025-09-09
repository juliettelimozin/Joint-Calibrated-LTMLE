#!/usr/bin R
###############################
# This script runs simulations with Q functions on logit form
###############################
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
load("/home/juliette/Calibrated-weights-sequential-trial-emulation/hers.Rdata")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/calibration_func_trials.R")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/hers_ltmle_utils.R")
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

# Data preparation
simdata <- HERS

simdata$Y <- as.factor(simdata$Y)
simdata$t <- simdata$visit - 8
simdata$SITE2 <- as.integer(simdata$SITE2)
simdata$SITE3 <- as.integer(simdata$SITE3)
simdata$WHITE <- as.integer(simdata$WHITE)
simdata$OTHER <- as.integer(simdata$OTHER)

simdata$CD4 <- (sqrt(as.numeric(simdata$CD4)) - mean(sqrt(simdata$CD4)))/sd(sqrt(simdata$CD4))
simdata$CD4_1 <- (sqrt(as.numeric(simdata$CD4_1)) - mean(sqrt(simdata$CD4_1)))/sd(sqrt(simdata$CD4_1))
simdata$CD4_2 <- (sqrt(as.numeric(simdata$CD4_2)) - mean(sqrt(simdata$CD4_2)))/sd(sqrt(simdata$CD4_2))

simdata$viral <- (log10(simdata$viral) - mean(log10(simdata$viral)))/sd(log10(simdata$viral))
simdata$viral_1 <- (log10(simdata$viral_1) - mean(log10(simdata$viral_1)))/sd(log10(simdata$viral_1))
simdata$viral_2 <- (log10(simdata$viral_2) - mean(log10(simdata$viral_2)))/sd(log10(simdata$viral_2))

simdata <- simdata %>% 
  dplyr::arrange(id,t) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(CAp = cumsum(as.numeric(haart_1 == 1 | haart_2 == 1))) %>% 
  dplyr::ungroup()

simdata$A <- simdata$haart
simdata$Ap <- simdata$haart_1
simdata$App <- simdata$haart_2
simdata[,'ID'] <- simdata$id
simdata <- simdata %>% 
  dplyr::select(ID, t, A, Ap, App,CAp, CD4, CD4_1,CD4_2,
                viral,viral_1,viral_2,HIVsym,HIVsym_1,HIVsym_2,
                SITE2, SITE3, WHITE, OTHER, Y, C)
simdata$eligible <- as.numeric(simdata$CAp == 0)
simdata$CA <-ave(simdata$A,simdata$ID,FUN=cumsum)

simdata <- as.data.frame(simdata)

# Point estimate 
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

registerDoParallel(cores = 10)
modified_bootstrap_LTMLE <- foreach(k = 1:bootstrap_iter, .combine=cbind) %dopar% {
  
  boot_design_data <- boot_samples[[k]]
  
  bootstrap_tmle <- MyLtmleMSM_modified(simdata = boot_design_data,initial_Q_t = point_estimate_modifiedTMLE$Q_t, refit_weights = FALSE)
  
  bootstrap_tmle$estimates
}

modified_bootstrap_CIs <- cbind(rowQuantiles(modified_bootstrap_LTMLE, probs = 0.025), rowQuantiles(modified_bootstrap_LTMLE, probs = 0.975))

############################# NORMAL BOOTSTRAP TMLE #############################

normal_bootstrap_LTMLE <-as.data.frame(matrix(,length(point_estimate_TMLE$estimates),bootstrap_iter))
normal_bootstrap_LTMLE <- foreach(k = 1:bootstrap_iter, .combine=cbind) %dopar% {
  
  boot_design_data <- as.data.frame(boot_samples[[k]])[, colnames(simdata)]
  boot_design_data$ID <- as.integer(factor(boot_design_data$ID))
  bootstrap_tmle <- MyLtmleMSM(simdata = boot_design_data)
  
  bootstrap_tmle$estimates
}

normal_bootstrap_CIs <- cbind(rowQuantiles(normal_bootstrap_LTMLE, probs = 0.025), rowQuantiles(normal_bootstrap_LTMLE, probs = 0.975))

