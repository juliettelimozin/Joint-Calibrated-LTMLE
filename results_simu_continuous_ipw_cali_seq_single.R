library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(modelr)
library(tidyverse)
library(tidyr)
library(MASS)
library(sandwich)
library(foreach)
library(doParallel)
library(parallel)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)
library(doRNG)
library(matrixStats)
library(rlist)
library(xtable)
library(janitor)
library(berryFunctions)

size <- c(200,500,1000)
treat <- c(-1,0,1)
conf <- c(1,1.5,2)

scenarios <- tidyr::crossing(size, treat,conf)
iters <- 500

bias_mr_all <- array(,dim = c(3*10,5,27))
bias_hr_all <- array(, dim = c(10,3,27))
sd_mr_all <- array(,dim = c(3*10,5,27))
sd_hr_all <- array(, dim = c(10,3,27))
mae_hr_all <-  array(, dim = c(10,3,27))
rootmse_mr_all <-  array(,dim = c(3*10,5,27))
rootmse_hr_all<- array(, dim = c(10,3,27))
for (i in 1:27){
  load(paste0("Simulation results/result_simu_continuous_ipw_cali_seq_single_", as.character(i), ".rda"))
  simu.t <- as.data.frame(tidyr::crossing(1:iters, 0:4))
  mr_all <- as.data.frame(list.rbind(oper[14,]))
  mr_all$simu <- simu.t[,1]
  mr_all$visit <- simu.t[,2]
  simu.scenario <- as.data.frame(tidyr::crossing(1:iters, 1:10))
  hr_all <- as.data.frame(list.rbind(oper[13,])) %>% 
    dplyr::mutate(simu = simu.scenario[,1], scenario =  simu.scenario[,2])
  for (t in 1:5){
    if (t == 1){
     bias_mr_all[,t,i] <- colMeans(mr_all[mr_all$visit == t-1,1:30],na.rm = T) - c(95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5)
    } else {
      bias_mr_all[,t,i] <- colMeans(mr_all[mr_all$visit == t-1,1:30],na.rm = T) - c(92,100,-8,
                                                                                    92,100,-8,
                                                                                    92,100,-8,
                                                                                    92,100,-8,
                                                                                    92,100,-8,
                                                                                    92,100,-8,
                                                                                    92,100,-8,
                                                                                    92,100,-8,
                                                                                    92,100,-8,
                                                                                    92,100,-8)
    }
    sd_mr_all[,t,i] <-colSds(as.matrix(mr_all[mr_all$visit == t-1,1:30]),na.rm = T)
    rootmse_mr_all[,t,i]<- sqrt(bias_mr_all[,t,i]^2 + sd_mr_all[,t,i]^2)
  }
  for (k in 1:10){
    bias_hr_all[k,,i]<- colMeans(hr_all[hr_all$scenario == k, 1:3], na.rm = TRUE) - c(100,-5,-8)
    sd_hr_all[k,,i]<- colSds(as.matrix(hr_all[hr_all$scenario == k, 1:3]), na.rm = TRUE)
    rootmse_hr_all[k,,i]<- sqrt(bias_hr_all[k,,i]^2 + sd_hr_all[k,,i]^2)
    mae_hr_all[k,,i] <- colMeans(abs(hr_all[hr_all$scenario == k, 1:3] %>% mutate(V1 = V1 - 100, V2 = V2+5, V3 = V3 + 8)), na.rm = TRUE)
  }
  
}

bias_hr_table <- bias_hr_all[,,1]
sd_hr_table <- sd_hr_all[,,1]
rootmse_hr_table <- rootmse_hr_all[,,1]
mae_hr_table <- mae_hr_all[,,1]
for (i in 2:27){
  bias_hr_table <- rbind(bias_hr_table,bias_hr_all[,,i])
  sd_hr_table <- rbind(sd_hr_table,sd_hr_all[,,i])
  rootmse_hr_table <- rbind(rootmse_hr_table,rootmse_hr_all[,,i])
  mae_hr_table <- rbind(mae_hr_table, mae_hr_all[,,i])
}

colnames(bias_hr_table) <- c('Bias 1', 'Bias 2', 'Bias 3')
colnames(sd_hr_table) <- c('SD 1', 'SD 2', 'SD 3')
colnames(rootmse_hr_table) <- c('rootMSE 1', 'rootMSE 2', 'rootMSE 3')
colnames(mae_hr_table) <- c('MAE 1', 'MAE 2', 'MAE 3')
size <- c(200,500,1000)
treat <- c(-1,0,1)
conf <- c(1,1.5,2)
var <- c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
         'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
scenarios <- tidyr::crossing(size, treat,conf,var) %>% 
  dplyr::arrange(size, treat, conf,
                 match(var, c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
                              'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
                 ))

hr_table <- cbind(scenarios,bias_hr_table, sd_hr_table,rootmse_hr_table, mae_hr_table) %>% 
  dplyr::select(size, treat, conf,var, 'Bias 1', 'SD 1', 'rootMSE 1', 'MAE 1','Bias 2', 'SD 2', 'rootMSE 2','MAE 2','Bias 3', 'SD 3', 'rootMSE 3', 'MAE 3')


print(xtable(hr_table[hr_table$size == 200,],digits=c(0,0,0,1,0,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')

print(xtable(hr_table[hr_table$size == 500,]),include.rownames=FALSE, type = 'latex')

print(xtable(hr_table[hr_table$size == 1000,]),include.rownames=FALSE, type = 'latex')

write_csv(hr_table, file = 'simulation_results_hr_table.csv')

bias_mr_table <- bias_mr_all[c(3,6,9,12,15,18,21,24,27,30),1:2,1]
sd_mr_table <- sd_mr_all[c(3,6,9,12,15,18,21,24,27,30),1:2,1]
rootmse_mr_table <- rootmse_mr_all[c(3,6,9,12,15,18,21,24,27,30),1:2,1]
for (i in 2:27){
  bias_mr_table <- rbind(bias_mr_table,bias_mr_all[c(3,6,9,12,15,18,21,24,27,30),1:2,i])
  sd_mr_table <- rbind(sd_mr_table,sd_mr_all[c(3,6,9,12,15,18,21,24,27,30),1:2,i])
  rootmse_mr_table <- rbind(rootmse_mr_table,rootmse_mr_all[c(3,6,9,12,15,18,21,24,27,30),1:2,i])
}

colnames(bias_mr_table) <- c('Bias 0', 'Bias k')
colnames(sd_mr_table) <- c('SD 0', 'SD k')
colnames(rootmse_mr_table) <- c('rootMSE 0', 'rootMSE k')
size <- c(200,500,1000)
treat <- c(-1,0,1)
conf <- c(1,1.5,2)
var <- c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
         'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
scenarios <- tidyr::crossing(size, treat,conf,var)%>% 
  dplyr::arrange(size, treat, conf,
                 match(var, c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
                              'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')))


mr_table <- cbind(scenarios,bias_mr_table, sd_mr_table,rootmse_mr_table) %>% 
  dplyr::select(size, treat, conf,var, 'Bias 0', 'SD 0', 'rootMSE 0','Bias k', 'SD k', 'rootMSE k')

write_csv(mr_table, file = 'simulation_results_mr_table.csv')
print(xtable(mr_table[mr_table$size == 200,]),include.rownames=FALSE, type = 'latex')

na_list <- c(6)
for (r in 2:54){
  na_list <- cbind(na_list, r*5 + r)
}
big_results_table <- hr_table %>% 
  insertRows(r = na_list, new = NA)

print(xtable(big_results_table[big_results_table$size == 200,],digits=c(0,0,0,1,0,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')

print(xtable(big_results_table[big_results_table$size == 500,],digits=c(0,0,0,1,0,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')

print(xtable(big_results_table[big_results_table$size == 1000,],digits=c(0,0,0,1,0,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')


wide_mr <-reshape(mr_table,
                  idvar = c('size', 'treat','conf'),
                  timevar = "var",
                  direction = "wide") %>% 
  clean_names() %>% 
  dplyr::mutate(diff_rootmseipw_cali_by_time = root_mse_k_ipw - root_mse_k_calibration_by_time,
                diff_rootmseipw_cali_by_time_miss = root_mse_k_ipw_miss - root_mse_k_calibration_by_time_miss) %>% 
  dplyr::filter(diff_rootmseipw_cali_by_time_miss >= 0)




