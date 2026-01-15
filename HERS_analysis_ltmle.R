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
library(data.table)
library(matrixStats)
library(xtable)
library(ggplot2)
library(ggpubr)
library(tidyr)

set.seed(160625)

# Data preparation
simdata <- HERS

simdata$C <- simdata$C + simdata$Y
simdata$t <- simdata$visit - 8
simdata$SITE2 <- as.integer(simdata$SITE2)
simdata$SITE3 <- as.integer(simdata$SITE3)
simdata$WHITE <- as.integer(simdata$WHITE)
simdata$OTHER <- as.integer(simdata$OTHER)

simdata$sqrtCD4 <- (sqrt(as.numeric(simdata$CD4)) - mean(sqrt(simdata$CD4)))/sd(sqrt(simdata$CD4))
simdata$sqrtCD4_1 <- (sqrt(as.numeric(simdata$CD4_1)) - mean(sqrt(simdata$CD4_1)))/sd(sqrt(simdata$CD4_1))
simdata$sqrtCD4_2 <- (sqrt(as.numeric(simdata$CD4_2)) - mean(sqrt(simdata$CD4_2)))/sd(sqrt(simdata$CD4_2))

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
  dplyr::select(ID, t, A, Ap, App,CAp, CD4, CD4_1,CD4_2,sqrtCD4, sqrtCD4_1,sqrtCD4_2,
                viral,viral_1,viral_2,HIVsym,HIVsym_1,HIVsym_2,
                SITE2, SITE3, WHITE, OTHER, Y, C)
simdata$eligible <- as.numeric(simdata$CAp == 0 & simdata$t == 0)
simdata$CA <-ave(simdata$A,simdata$ID,FUN=cumsum)
simdata <- as.data.frame(simdata)
simdata$eligible <- ave(simdata$eligible,simdata$ID,FUN=cumsum)

simdata <- simdata[simdata$eligible ==1,]

simdata$RA <- 1
simdata[simdata$t == 0 & !(simdata$CA == 1),]$RA <- 0
simdata[simdata$t == 1 & !(simdata$CA == 2),]$RA <- 0
simdata[simdata$t == 2 & !(simdata$CA == 3),]$RA <- 0
simdata[simdata$t == 3 & !(simdata$CA == 4),]$RA <- 0
simdata[simdata$t == 4 & !(simdata$CA == 5),]$RA <- 0

tab <- xtabs(RA ~ t, data = simdata)

print(tab)

simdata$RA <- 1
simdata[simdata$t == 0 & !(simdata$CA == 0),]$RA <- 0
simdata[simdata$t == 1 & !(simdata$CA == 0),]$RA <- 0
simdata[simdata$t == 2 & !(simdata$CA == 0),]$RA <- 0
simdata[simdata$t == 3 & !(simdata$CA == 0),]$RA <- 0
simdata[simdata$t == 4 & !(simdata$CA == 0),]$RA <- 0

tab <- xtabs(RA ~ t, data = simdata)

print(tab)

simdata <- simdata %>% 
  group_by(ID) %>% 
  mutate(D = factor(case_when(
    CD4_1[t == 0][1] < 200              ~ 0,
    CD4_1[t == 0][1] >= 200 & CD4_1[t == 0][1] <= 500 ~ 1,
    CD4_1[t == 0][1] > 500             ~ 2
  ),
  levels = c(0,1,2))
  ) %>%
  ungroup()

simdata <- as.data.frame(simdata)

################# MODIFIED TMLE ######################################
# Point estimate 
point_estimate_modifiedTMLE <- MyLtmleMSM_modified_hers(simdata, initial_Q_t = NA, refit_weights = TRUE)

simdata_with_weights <- point_estimate_modifiedTMLE$data_with_weights
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

registerDoParallel(cores = 10)

time <- proc.time()

modified_bootstrap_LTMLE <- foreach(k = 1:bootstrap_iter, .combine = function(x, y) {
  # combine list of results
  list(
    coeffs = rbind(x$coeffs, y$coeffs)
  )
}, .init = list(coeffs = NULL)) %dopar% {

  boot_design_data <- boot_samples[[k]]

  bootstrap_tmle <- MyLtmleMSM_modified_hers(simdata = boot_design_data,
                                             initial_Q_t = point_estimate_modifiedTMLE$Q_fits,
                                             refit_weights = FALSE,
                                             treat_models_n = point_estimate_modifiedTMLE$treat_models_n,
                                             censor_models_n = point_estimate_modifiedTMLE$censor_models_n)

  list(coeffs = bootstrap_tmle$MSM_estimates)
}


coeffs_boot <- modified_bootstrap_LTMLE$coeffs

modified_bootstrap_CIs_coeffs_LB <- rbind(colQuantiles(coeffs_boot[rownames(coeffs_boot) == "MLE LTMLE",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))

rownames(modified_bootstrap_CIs_coeffs_LB) <- c('MLE LTMLE',
                                                'CMLE LTMLE treat. only',
                                                'Aggr. CMLE LTMLE treat. only',
                                                'CMLE LTMLE',
                                                'Aggr. CMLE LTMLE')

modified_bootstrap_CIs_coeffs_UB <- rbind(colQuantiles(coeffs_boot[rownames(coeffs_boot) == "MLE LTMLE",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))
rownames(modified_bootstrap_CIs_coeffs_UB) <- c('MLE LTMLE',
                                                'CMLE LTMLE treat. only',
                                                'Aggr. CMLE LTMLE treat. only',
                                                'CMLE LTMLE',
                                                'Aggr. CMLE LTMLE')


est <- point_estimate_modifiedTMLE$MSM_estimates[,(ncol(point_estimate_modifiedTMLE$MSM_estimates)-2):ncol(point_estimate_modifiedTMLE$MSM_estimates)]
lower <- modified_bootstrap_CIs_coeffs_LB[,(ncol(modified_bootstrap_CIs_coeffs_LB)-2):ncol(modified_bootstrap_CIs_coeffs_LB)]
upper <- modified_bootstrap_CIs_coeffs_UB[,(ncol(modified_bootstrap_CIs_coeffs_UB)-2):ncol(modified_bootstrap_CIs_coeffs_UB)]

format_matrix <- function(est, lower, upper) {
  m <- matrix("", nrow = nrow(est), ncol = ncol(est),
              dimnames = list(rownames(est), colnames(est)))
  for (i in 1:nrow(est)) {
    for (j in 1:ncol(est)) {
      m[i, j] <- sprintf("%.2f (%.2f, %.2f)", est[i, j], lower[i, j], upper[i, j])
    }
  }
  m
}

summary_table <- format_matrix(est, lower, upper)
print(xtable(summary_table, caption = "Estimates and 95% CIs by modified LTMLE"))

print('Modified bootstrap CI time \n')
print((proc.time() - time)[[3]])
################# NORMAL TMLE ######################################
# Point estimate 
point_estimate_TMLE <- MyLtmleMSM_hers(simdata)

############################# Normal BOOTSTRAP TMLE #############################
registerDoParallel(cores = 10)

time <- proc.time()
normal_bootstrap_LTMLE <- foreach(k = 1:bootstrap_iter, .combine = function(x, y) {
  # combine list of results
  list(
    coeffs = rbind(x$coeffs, y$coeffs)
  )
}, .init = list(coeffs = NULL)) %dopar% {

  boot_design_data <- as.data.frame(boot_samples[[k]])[,names(simdata)]
  boot_design_data$ID <- as.integer(factor(boot_design_data$ID))
  bootstrap_tmle <- MyLtmleMSM_hers(simdata = boot_design_data)

  list(coeffs = bootstrap_tmle$MSM_estimates)
}

coeffs_boot <- normal_bootstrap_LTMLE$coeffs

normal_bootstrap_CIs_coeffs_LB <- rbind(colQuantiles(coeffs_boot[rownames(coeffs_boot) == "MLE LTMLE",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))

rownames(normal_bootstrap_CIs_coeffs_LB) <- c('MLE LTMLE',
                                                'CMLE LTMLE treat. only',
                                                'Aggr. CMLE LTMLE treat. only',
                                                'CMLE LTMLE',
                                                'Aggr. CMLE LTMLE')

normal_bootstrap_CIs_coeffs_UB <- rbind(colQuantiles(coeffs_boot[rownames(coeffs_boot) == "MLE LTMLE",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))
rownames(normal_bootstrap_CIs_coeffs_UB) <- c('MLE LTMLE',
                                                'CMLE LTMLE treat. only',
                                                'Aggr. CMLE LTMLE treat. only',
                                                'CMLE LTMLE',
                                                'Aggr. CMLE LTMLE')

est <- point_estimate_TMLE$MSM_estimates[,(ncol(point_estimate_TMLE$MSM_estimates)-2):ncol(point_estimate_TMLE$MSM_estimates)]
lower <- normal_bootstrap_CIs_coeffs_LB[,(ncol(normal_bootstrap_CIs_coeffs_LB)-2):ncol(normal_bootstrap_CIs_coeffs_LB)]
upper <- normal_bootstrap_CIs_coeffs_UB[,(ncol(normal_bootstrap_CIs_coeffs_UB)-2):ncol(normal_bootstrap_CIs_coeffs_UB)]

format_matrix <- function(est, lower, upper) {
  m <- matrix("", nrow = nrow(est), ncol = ncol(est),
              dimnames = list(rownames(est), colnames(est)))
  for (i in 1:nrow(est)) {
    for (j in 1:ncol(est)) {
      m[i, j] <- sprintf("%.2f (%.2f, %.2f)", est[i, j], lower[i, j], upper[i, j])
    }
  }
  m
}

summary_table <- format_matrix(est, lower, upper)
print(xtable(summary_table, caption = "Estimates and 95% CIs by normal LTMLE"))
print('Normal bootstrap CI time \n')
print((proc.time() - time)[[3]])

################# Cali Boot TMLE ######################################
# Point estimate 
point_estimate_cali_boot_TMLE <- MyLtmleMSM_cali_boot_hers(simdata)

############################# Cali Boot BOOTSTRAP TMLE #############################
registerDoParallel(cores = 10)
time <- proc.time()
normal_bootstrap_LTMLE <- foreach(k = 1:bootstrap_iter, .combine = function(x, y) {
  # combine list of results
  list(
    coeffs = rbind(x$coeffs, y$coeffs),
    ATE_strata0 = rbind(x$ATE_strata0, y$ATE_strata0),
    ATE_strata1 = rbind(x$ATE_strata1, y$ATE_strata1),
    ATE_strata2 = rbind(x$ATE_strata2, y$ATE_strata2)
  )
}, .init = list(coeffs = NULL, ATE_strata0 = NULL,ATE_strata1 = NULL,ATE_strata2 = NULL)) %dopar% {
  
  boot_design_data <- as.data.frame(boot_samples[[k]])[,c(names(simdata), "weight")]
  boot_design_data$ID <- as.integer(factor(boot_design_data$ID))
  bootstrap_tmle <- MyLtmleMSM_cali_boot_hers(simdata = boot_design_data,
                                              initial_Q_t = point_estimate_cali_boot_TMLE$Q_fits, 
                                              refit_weights = FALSE,
                                              treat_models_n = point_estimate_cali_boot_TMLE$treat_models_n,
                                              censor_models_n = point_estimate_cali_boot_TMLE$censor_models_n)
  
  list(coeffs = bootstrap_tmle$MSM_estimates,
       ATE_strata0 = bootstrap_tmle$counterfactual_means$EY_always_treated_strata0 - bootstrap_tmle$counterfactual_means$EY_never_treated_strata0,
       ATE_strata1 = bootstrap_tmle$counterfactual_means$EY_always_treated_strata1 - bootstrap_tmle$counterfactual_means$EY_never_treated_strata1,
       ATE_strata2 = bootstrap_tmle$counterfactual_means$EY_always_treated_strata2 - bootstrap_tmle$counterfactual_means$EY_never_treated_strata2
       )
}

coeffs_boot <- normal_bootstrap_LTMLE$coeffs
ATE_strata0_boot <- normal_bootstrap_LTMLE$ATE_strata0
ATE_strata1_boot <- normal_bootstrap_LTMLE$ATE_strata1
ATE_strata2_boot <- normal_bootstrap_LTMLE$ATE_strata2


normal_bootstrap_CIs_coeffs_LB <- rbind(colQuantiles(coeffs_boot[rownames(coeffs_boot) == "MLE LTMLE",], probs = 0.025),
                                        colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                        colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                        colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE",], probs = 0.025),
                                        colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))


normal_bootstrap_CIs_coeffs_UB <- rbind(colQuantiles(coeffs_boot[rownames(coeffs_boot) == "MLE LTMLE",], probs = 0.975),
                                        colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                        colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                        colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE",], probs = 0.975),
                                        colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))

normal_bootstrap_CIs_ATE_strata0_LB <- rbind(colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "MLE LTMLE",], probs = 0.025),
                                             colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                             colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                             colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "CMLE LTMLE",], probs = 0.025),
                                             colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))


normal_bootstrap_CIs_ATE_strata0_UB <- rbind(colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "MLE LTMLE",], probs = 0.975),
                                             colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                             colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                             colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "CMLE LTMLE",], probs = 0.975),
                                             colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))

normal_bootstrap_CIs_ATE_strata1_LB <- rbind(colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "MLE LTMLE",], probs = 0.025),
                                             colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                             colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                             colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "CMLE LTMLE",], probs = 0.025),
                                             colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))


normal_bootstrap_CIs_ATE_strata1_UB <- rbind(colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "MLE LTMLE",], probs = 0.975),
                                             colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                             colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                             colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "CMLE LTMLE",], probs = 0.975),
                                             colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))

normal_bootstrap_CIs_ATE_strata2_LB <- rbind(colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "MLE LTMLE",], probs = 0.025),
                                             colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                             colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                             colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "CMLE LTMLE",], probs = 0.025),
                                             colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))


normal_bootstrap_CIs_ATE_strata2_UB <- rbind(colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "MLE LTMLE",], probs = 0.975),
                                             colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                             colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                             colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "CMLE LTMLE",], probs = 0.975),
                                             colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))

rownames(normal_bootstrap_CIs_coeffs_LB) <- rownames(normal_bootstrap_CIs_coeffs_UB) <- 
  rownames(normal_bootstrap_CIs_ATE_strata0_LB) <- rownames(normal_bootstrap_CIs_ATE_strata0_UB) <-
  rownames(normal_bootstrap_CIs_ATE_strata1_LB) <- rownames(normal_bootstrap_CIs_ATE_strata1_UB) <-
  rownames(normal_bootstrap_CIs_ATE_strata2_LB) <- rownames(normal_bootstrap_CIs_ATE_strata2_UB) <- c('MLE LTMLE',
                                                                                                      'CMLE LTMLE treat. only',
                                                                                                      'Aggr. CMLE LTMLE treat. only',
                                                                                                      'CMLE LTMLE',
                                                                                                      'Aggr. CMLE LTMLE')


est <- point_estimate_TMLE$MSM_estimates[,(ncol(point_estimate_TMLE$MSM_estimates)-2):ncol(point_estimate_TMLE$MSM_estimates)]
lower <- normal_bootstrap_CIs_coeffs_LB[,(ncol(normal_bootstrap_CIs_coeffs_LB)-2):ncol(normal_bootstrap_CIs_coeffs_LB)]
upper <- normal_bootstrap_CIs_coeffs_UB[,(ncol(normal_bootstrap_CIs_coeffs_UB)-2):ncol(normal_bootstrap_CIs_coeffs_UB)]

format_matrix <- function(est, lower, upper) {
  m <- matrix("", nrow = nrow(est), ncol = ncol(est), 
              dimnames = list(rownames(est), colnames(est)))
  for (i in 1:nrow(est)) {
    for (j in 1:ncol(est)) {
      m[i, j] <- sprintf("%.2f (%.2f, %.2f)", est[i, j], lower[i, j], upper[i, j])
    }
  }
  m
}

summary_table <- format_matrix(est, lower, upper)
print(xtable(summary_table, caption = "Estimates and 95% CIs by normal LTMLE with cali boot"))
print('Cali bootstrap CI time \n')
print((proc.time() - time)[[3]])

point_estimate_TMLE$counterfactual_means
outcome_plots <- list()

outcome_plots[[1]] <- ggplot() +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata0[1,], colour = "Strata 0", shape = 'Always-treated')) +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata1[1,], colour = "Strata 1", shape = 'Always-treated')) +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata2[1,], colour = "Strata 2", shape = 'Always-treated')) + 
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_never_treated_strata0[1,], colour = "Strata 0", shape = 'Never-treated')) +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_never_treated_strata1[1,], colour = "Strata 1", shape = 'Never-treated')) +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_never_treated_strata2[1,], colour = "Strata 2", shape = 'Never-treated')) + 
  scale_color_manual(name = "Stratum", values = c("Strata 0"= "red", "Strata 1" = "blue",
                                                  "Strata 2" = "orange"), labels = c("0","1", "2")) +
  scale_shape_manual(name = 'Treatment protocol', values = c('Never-treated' = 1, 'Always-treated' = 2)) +
  theme(legend.text = element_text(size=14),
        title = element_text(size=12)) +
  ylim(100,700) + 
  labs(title = 'MLE', x = 'Visit', y = 'Mean counterfactual CD4 cell count by strata')

outcome_plots[[2]] <- ggplot() +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata0[4,], colour = "Strata 0", shape = 'Always-treated')) +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata1[4,], colour = "Strata 1", shape = 'Always-treated')) +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata2[4,], colour = "Strata 2", shape = 'Always-treated')) + 
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_never_treated_strata0[4,], colour = "Strata 0", shape = 'Never-treated')) +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_never_treated_strata1[4,], colour = "Strata 1", shape = 'Never-treated')) +
  geom_point(aes(x = 8:12, y = point_estimate_TMLE$counterfactual_means$EY_never_treated_strata2[4,], colour = "Strata 2", shape = 'Never-treated')) + 
  scale_color_manual(name = "Stratum", values = c("Strata 0"= "red", "Strata 1" = "blue",
                                                 "Strata 2" = "orange"), labels = c("0","1", "2")) +
  scale_shape_manual(name = 'Treatment protocol', values = c('Never-treated' = 1, 'Always-treated' = 2)) +
  theme(legend.text = element_text(size=14),
        title = element_text(size=12)) +
  ylim(100,700) + 
  labs(title = 'CMLE', x = 'Visit', y = 'Mean counterfactual CD4 cell count by strata')
  
png('HERS_outcome_plots.png', width = 10, height = 5, units = 'in', res = 300)
ggarrange(plotlist = outcome_plots, nrow = 1, ncol = 2, common.legend = T,
          legend = 'bottom')
dev.off()


CI_MLE_ATE <- data.frame(Visit = 8:12,ATE_strata0 = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata0[1,] - point_estimate_TMLE$counterfactual_means$EY_never_treated_strata0[1,],
                         ATE_strata1 = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata1[1,] - point_estimate_TMLE$counterfactual_means$EY_never_treated_strata1[1,],
                         ATE_strata2 = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata2[1,] - point_estimate_TMLE$counterfactual_means$EY_never_treated_strata2[1,],
                         LB_strata0 = normal_bootstrap_CIs_ATE_strata0_LB[1,],
                         UB_strata0 = normal_bootstrap_CIs_ATE_strata0_UB[1,],
                         LB_strata1 = normal_bootstrap_CIs_ATE_strata1_LB[1,],
                         UB_strata1 = normal_bootstrap_CIs_ATE_strata1_UB[1,],
                         LB_strata2 = normal_bootstrap_CIs_ATE_strata2_LB[1,],
                         UB_strata2 = normal_bootstrap_CIs_ATE_strata2_UB[1,]
                                  )

CI_MLE_ATE_long <- pivot_longer(data = CI_MLE_ATE,cols = -Visit,
                                names_to = c("Stat", "Strata"),
                                names_pattern = "(.*)_(strata\\d+)",
                                values_to = "Value"
                                )
CI_MLE_ATE <- pivot_wider(CI_MLE_ATE_long,names_from = Stat, values_from = Value)

CI_CMLE_ATE <- data.frame(Visit = 8:12,ATE_strata0 = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata0[4,] - point_estimate_TMLE$counterfactual_means$EY_never_treated_strata0[4,],
                         ATE_strata1 = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata1[4,] - point_estimate_TMLE$counterfactual_means$EY_never_treated_strata1[4,],
                         ATE_strata2 = point_estimate_TMLE$counterfactual_means$EY_always_treated_strata2[4,] - point_estimate_TMLE$counterfactual_means$EY_never_treated_strata2[4,],
                         LB_strata0 = normal_bootstrap_CIs_ATE_strata0_LB[4,],
                         UB_strata0 = normal_bootstrap_CIs_ATE_strata0_UB[4,],
                         LB_strata1 = normal_bootstrap_CIs_ATE_strata1_LB[4,],
                         UB_strata1 = normal_bootstrap_CIs_ATE_strata1_UB[4,],
                         LB_strata2 = normal_bootstrap_CIs_ATE_strata2_LB[4,],
                         UB_strata2 = normal_bootstrap_CIs_ATE_strata2_UB[4,]
)

CI_CMLE_ATE_long <- pivot_longer(data = CI_CMLE_ATE,cols = -Visit,
                                names_to = c("Stat", "Strata"),
                                names_pattern = "(.*)_(strata\\d+)",
                                values_to = "Value"
)
CI_CMLE_ATE <- pivot_wider(CI_CMLE_ATE_long,names_from = Stat, values_from = Value)

ATE_plots <- list()
ATE_plots[[1]] <- ggplot(CI_MLE_ATE, aes(x = Visit,colour = Strata)) +
  geom_point(aes(y = ATE),position = position_dodge(width = 0.2)) + 
  geom_errorbar(aes(ymin = LB, ymax = UB), width = 0.2,position = position_dodge(width = 0.2)) +
  scale_color_manual(name = "Stratum", values = c("strata0"= "red", "strata1" = "blue",
                                                 "strata2" = "orange"), labels = c("0","1", "2")) +
  theme(legend.text = element_text(size=14),
        title = element_text(size=12)) +
  ylim(-150,350) + 
  labs(title = 'MLE', x = 'Visit', y = 'ATE by strata')

ATE_plots[[2]] <- ggplot(CI_CMLE_ATE, aes(x = Visit,colour = Strata)) +
  geom_point(aes(y = ATE),position = position_dodge(width = 0.2)) + 
  geom_errorbar(aes(ymin = LB, ymax = UB), width = 0.2,position = position_dodge(width = 0.2)) +
  scale_color_manual(name = "Stratum", values = c("strata0"= "red", "strata1" = "blue",
                                                 "strata2" = "orange"), labels = c("0","1", "2")) +
  theme(legend.text = element_text(size=14),
        title = element_text(size=12)) +
  ylim(-150,350) + 
  labs(title = 'CMLE', x = 'Visit', y = 'ATE by strata')

png('HERS_ATE_plots.png', width = 10, height = 5, units = 'in', res = 300)
ggarrange(plotlist = ATE_plots, nrow = 1, ncol = 2, common.legend = T,
          legend = 'bottom')
dev.off()
