source("dgm_2nd_simulation_biometrics.R")
source("calibration_func_trials.R")
source("ltmle_utils.R")
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


for(i in c(300,500,1000,2500)){
  coverage <- matrix(,9,3)
  rownames(coverage) <- c('MLE LTMLE beta_0',
                          'CMLE LTMLE beta_0',
                          'Aggr. CMLE LTMLE beta_0',
                          'MLE LTMLE beta_1',
                          'CMLE LTMLE beta_1',
                          'Aggr. CMLE LTMLE beta_1',
                          'MLE LTMLE E(Y^1_2)',
                          'CMLE LTMLE E(Y^1_2)',
                          'Aggr. CMLE LTMLE E(Y^1_2)')
  colnames(coverage) <- c('Modified boot.', 'Modified boot. with cali', 'Normal boot.')
  
  ci_coverage_simu_result <- readRDS(paste0("ci_coverage_simu_result_",i,".rds"))
  mle_comp_result <- readRDS(paste0("mle_ci_comp_time_weak_conf_result_",i,".rds"))
  
  thr <- c(200,200,200,10,10,10,230,230,230)
  indLB <- sweep(ci_coverage_simu_result$modified_bootstrap_CIs_LB, 1, thr, FUN = "<=")
  indUB <- sweep(ci_coverage_simu_result$modified_bootstrap_CIs_UB, 1, thr, FUN = ">=")
  
  ind <- indLB & indUB
  
  coverage[,1] <- rowSums(ind)/1000
  
  indLB <- sweep(ci_coverage_simu_result$modified_cali_bootstrap_CIs_LB, 1, thr, FUN = "<=")
  indUB <- sweep(ci_coverage_simu_result$modified_cali_bootstrap_CIs_UB, 1, thr, FUN = ">=")
  
  ind <- indLB & indUB
  
  coverage[,2] <- rowSums(ind)/1000
  
  indLB <- sweep(ci_coverage_simu_result$normal_bootstrap_CIs_LB, 1, thr, FUN = "<=")
  indUB <- sweep(ci_coverage_simu_result$normal_bootstrap_CIs_UB, 1, thr, FUN = ">=")
  
  ind <- indLB & indUB
  
  coverage[,3] <- rowSums(ind)/1000
  
  
  width <- coverage
  
  width[,1] <- rowMeans(ci_coverage_simu_result$modified_bootstrap_CIs_UB - ci_coverage_simu_result$modified_bootstrap_CIs_LB)
  width[,2] <- rowMeans(ci_coverage_simu_result$modified_cali_bootstrap_CIs_UB - ci_coverage_simu_result$modified_cali_bootstrap_CIs_LB)
  width[,3] <- rowMeans(ci_coverage_simu_result$normal_bootstrap_CIs_UB - ci_coverage_simu_result$normal_bootstrap_CIs_LB)
  
  computation_time <- rowMeans(ci_coverage_simu_result$computation_time)
  mle_computation_time <- rowMeans(mle_comp_result$computation_time)
  
  names(computation_time) <- c('Modified boot.', 'Modified boot. with cali', 'Normal boot.')
  
  coverage <- t(coverage)
  
  coverage <- rbind(coverage[c(3,1),c('MLE LTMLE beta_0',
                                'MLE LTMLE beta_1')],
                    coverage[c(3,1,2),c('CMLE LTMLE beta_0',
                          'CMLE LTMLE beta_1')])
  
  width <- t(width)
  
  width <- rbind(width[c(3,1),c('MLE LTMLE beta_0',
                             'MLE LTMLE beta_1')],
                 width[c(3,1,2),c('CMLE LTMLE beta_0',
                             'CMLE LTMLE beta_1')])
  computation_time <- c(mle_computation_time[c(2,1)],computation_time[c(3,1,2)])
  
  results_table <- cbind(coverage,width, computation_time)
  
  ########################
  
  coverage <- matrix(,9,3)
  rownames(coverage) <- c('MLE LTMLE beta_0',
                          'CMLE LTMLE beta_0',
                          'Aggr. CMLE LTMLE beta_0',
                          'MLE LTMLE beta_1',
                          'CMLE LTMLE beta_1',
                          'Aggr. CMLE LTMLE beta_1',
                          'MLE LTMLE E(Y^1_2)',
                          'CMLE LTMLE E(Y^1_2)',
                          'Aggr. CMLE LTMLE E(Y^1_2)')
  colnames(coverage) <- c('Modified boot.', 'Modified boot. with cali', 'Normal boot.')
  
  ci_coverage_simu_result <- readRDS(paste0("ci_coverage_simu_strong_conf_result_",i,".rds"))
  mle_comp_result <- readRDS(paste0("mle_ci_comp_time_strong_conf_result_",i,".rds"))
  
  thr <- c(200,200,200,10,10,10,230,230,230)
  indLB <- sweep(ci_coverage_simu_result$modified_bootstrap_CIs_LB, 1, thr, FUN = "<=")
  indUB <- sweep(ci_coverage_simu_result$modified_bootstrap_CIs_UB, 1, thr, FUN = ">=")
  
  ind <- indLB & indUB
  
  coverage[,1] <- rowSums(ind)/1000
  
  indLB <- sweep(ci_coverage_simu_result$modified_cali_bootstrap_CIs_LB, 1, thr, FUN = "<=")
  indUB <- sweep(ci_coverage_simu_result$modified_cali_bootstrap_CIs_UB, 1, thr, FUN = ">=")
  
  ind <- indLB & indUB
  
  coverage[,2] <- rowSums(ind)/1000
  
  indLB <- sweep(ci_coverage_simu_result$normal_bootstrap_CIs_LB, 1, thr, FUN = "<=")
  indUB <- sweep(ci_coverage_simu_result$normal_bootstrap_CIs_UB, 1, thr, FUN = ">=")
  
  ind <- indLB & indUB
  
  coverage[,3] <- rowSums(ind)/1000
  
  
  width <- coverage
  
  width[,1] <- rowMeans(ci_coverage_simu_result$modified_bootstrap_CIs_UB - ci_coverage_simu_result$modified_bootstrap_CIs_LB)
  width[,2] <- rowMeans(ci_coverage_simu_result$modified_cali_bootstrap_CIs_UB - ci_coverage_simu_result$modified_cali_bootstrap_CIs_LB)
  width[,3] <- rowMeans(ci_coverage_simu_result$normal_bootstrap_CIs_UB - ci_coverage_simu_result$normal_bootstrap_CIs_LB)
  
  computation_time <- rowMeans(ci_coverage_simu_result$computation_time)
  mle_computation_time <- rowMeans(mle_comp_result$computation_time)
  
  names(computation_time) <- c('Modified boot.', 'Modified boot. with cali', 'Normal boot.')
  
  coverage <- t(coverage)
  
  coverage <- rbind(coverage[c(3,1),c('MLE LTMLE beta_0',
                                      'MLE LTMLE beta_1')],
                    coverage[c(3,1,2),c('CMLE LTMLE beta_0',
                                'CMLE LTMLE beta_1')])
  
  width <- t(width)
  
  width <- rbind(width[c(3,1),c('MLE LTMLE beta_0',
                                'MLE LTMLE beta_1')],
                 width[c(3,1,2),c('CMLE LTMLE beta_0',
                          'CMLE LTMLE beta_1')])
  computation_time <- c(mle_computation_time[c(2,1)],computation_time[c(3,1,2)])
  
  results_table <- cbind(results_table, cbind(coverage,width, computation_time))
  
  rownames(results_table)<- c('MLE Normal.boot.','MLE Modified.boot.', 'CMLE Normal boot.', 'CMLE Modified boot.', 'CMLE Modified boot. with cali')
  
  print(xtable(results_table, type = 'latex'))
}




