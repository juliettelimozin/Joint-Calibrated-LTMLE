#!/usr/bin R
###############################
# This script runs simulations with Q functions on logit form
###############################
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
load("/home/juliette/hers.Rdata")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/calibration_func_trials.R")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/hers_ltmle_utils_SL.R")
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
library(SuperLearner)
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

point_estimate_TMLE <- MyLtmleMSM_hers(simdata)

point_estimate_cali_boot_TMLE <- MyLtmleMSM_cali_boot_hers(simdata)

simdata_with_weights <- point_estimate_cali_boot_TMLE$data_with_weights

SL <- list("SL.glm",
     "SL.stepAIC", "SL.bayesglm", c("SL.glm","screen.corP"), c("SL.step", "screen.corP"),
     c("SL.step.forward","screen.corP"), c("SL.stepAIC", "screen.corP"), c("SL.step.interaction","screen.corP"),c("SL.bayesglm", "screen.corP"))
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

############################# Cali Boot BOOTSTRAP TMLE #############################
registerDoParallel(cores = 12)
time <- proc.time()
cali_bootstrap_modLTMLE <- foreach(k = 1:bootstrap_iter, .combine = function(x, y) {
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
                                              treat_model_n = point_estimate_cali_boot_TMLE$treat_model_n,
                                              cense_model_n = point_estimate_cali_boot_TMLE$cense_model_n)
  
  list(coeffs = bootstrap_tmle$MSM_estimates,
       ATE_strata0 = bootstrap_tmle$counterfactual_means$EY_always_treated_strata0 - bootstrap_tmle$counterfactual_means$EY_never_treated_strata0,
       ATE_strata1 = bootstrap_tmle$counterfactual_means$EY_always_treated_strata1 - bootstrap_tmle$counterfactual_means$EY_never_treated_strata1,
       ATE_strata2 = bootstrap_tmle$counterfactual_means$EY_always_treated_strata2 - bootstrap_tmle$counterfactual_means$EY_never_treated_strata2
  )
}


coeffs_boot <- cali_bootstrap_modLTMLE$coeffs
ATE_strata0_boot <- cali_bootstrap_modLTMLE$ATE_strata0
ATE_strata1_boot <- cali_bootstrap_modLTMLE$ATE_strata1
ATE_strata2_boot <- cali_bootstrap_modLTMLE$ATE_strata2


cali_bootstrap_mod_CIs_coeffs_LB <- rbind(colQuantiles(coeffs_boot[rownames(coeffs_boot) == "MLE LTMLE",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE",], probs = 0.025),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))


cali_bootstrap_mod_CIs_coeffs_UB <- rbind(colQuantiles(coeffs_boot[rownames(coeffs_boot) == "MLE LTMLE",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "CMLE LTMLE",], probs = 0.975),
                                          colQuantiles(coeffs_boot[rownames(coeffs_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))

cali_bootstrap_mod_CIs_ATE_strata0_LB <- rbind(colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "MLE LTMLE",], probs = 0.025),
                                               colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                               colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                               colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "CMLE LTMLE",], probs = 0.025),
                                               colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))


cali_bootstrap_mod_CIs_ATE_strata0_UB <- rbind(colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "MLE LTMLE",], probs = 0.975),
                                               colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                               colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                               colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "CMLE LTMLE",], probs = 0.975),
                                               colQuantiles(ATE_strata0_boot[rownames(ATE_strata0_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))

cali_bootstrap_mod_CIs_ATE_strata1_LB <- rbind(colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "MLE LTMLE",], probs = 0.025),
                                               colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                               colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                               colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "CMLE LTMLE",], probs = 0.025),
                                               colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))


cali_bootstrap_mod_CIs_ATE_strata1_UB <- rbind(colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "MLE LTMLE",], probs = 0.975),
                                               colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                               colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                               colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "CMLE LTMLE",], probs = 0.975),
                                               colQuantiles(ATE_strata1_boot[rownames(ATE_strata1_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))

cali_bootstrap_mod_CIs_ATE_strata2_LB <- rbind(colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "MLE LTMLE",], probs = 0.025),
                                               colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "CMLE LTMLE treat. only",], probs = 0.025),
                                               colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.025),
                                               colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "CMLE LTMLE",], probs = 0.025),
                                               colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "Aggr. CMLE LTMLE",], probs = 0.025))


cali_bootstrap_mod_CIs_ATE_strata2_UB <- rbind(colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "MLE LTMLE",], probs = 0.975),
                                               colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "CMLE LTMLE treat. only",], probs = 0.975),
                                               colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "Aggr. CMLE LTMLE treat. only",], probs = 0.975),
                                               colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "CMLE LTMLE",], probs = 0.975),
                                               colQuantiles(ATE_strata2_boot[rownames(ATE_strata2_boot) == "Aggr. CMLE LTMLE",], probs = 0.975))

rownames(cali_bootstrap_mod_CIs_coeffs_LB) <- rownames(cali_bootstrap_mod_CIs_coeffs_UB) <- 
  rownames(cali_bootstrap_mod_CIs_ATE_strata0_LB) <- rownames(cali_bootstrap_mod_CIs_ATE_strata0_UB) <-
  rownames(cali_bootstrap_mod_CIs_ATE_strata1_LB) <- rownames(cali_bootstrap_mod_CIs_ATE_strata1_UB) <-
  rownames(cali_bootstrap_mod_CIs_ATE_strata2_LB) <- rownames(cali_bootstrap_mod_CIs_ATE_strata2_UB) <- c('MLE LTMLE',
                                                                                                          'CMLE LTMLE treat. only',
                                                                                                          'Aggr. CMLE LTMLE treat. only',
                                                                                                          'CMLE LTMLE',
                                                                                                          'Aggr. CMLE LTMLE')


est <- point_estimate_TMLE$MSM_estimates[,(ncol(point_estimate_TMLE$MSM_estimates)-2):ncol(point_estimate_TMLE$MSM_estimates)]
lower <- cali_bootstrap_mod_CIs_coeffs_LB[,(ncol(cali_bootstrap_mod_CIs_coeffs_LB)-2):ncol(cali_bootstrap_mod_CIs_coeffs_LB)]
upper <- cali_bootstrap_mod_CIs_coeffs_UB[,(ncol(cali_bootstrap_mod_CIs_coeffs_UB)-2):ncol(cali_bootstrap_mod_CIs_coeffs_UB)]

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

summary_table <- format_matrix(est[c(1,4),], lower[c(1,4),], upper[c(1,4),])
print(xtable(summary_table, caption = "Estimates and 95% CIs by normal LTMLE with cali boot"))
print('Cali bootstrap CI time \n')

save.image(file = 'HERS.RData')
