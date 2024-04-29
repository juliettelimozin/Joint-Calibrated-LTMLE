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
treat <- c(0.25,0.5,0.75)
conf <- c(1,1.5,2)

scenarios <- tidyr::crossing(size, treat,conf)
iters <- 500

bias_ate_all <- array(,dim = c(10,5,27))
bias_mean1_all <- array(,dim = c(10,5,27))
bias_hr_all <- array(, dim = c(10,3,27))
sd_ate_all <- array(,dim = c(10,5,27))
sd_hr_all <- array(, dim = c(10,3,27))
sd_mean1_all <- array(,dim = c(10,5,27))

mae_ate_all <- array(,dim = c(10,5,27))
mae_hr_all <-  array(, dim = c(10,3,27))
mae_mean1_all <- array(,dim = c(10,5,27))

medae_ate_all <- array(,dim = c(10,5,27))
medae_mean1_all <- array(,dim = c(10,5,27))
medae_hr_all <-  array(, dim = c(10,3,27))

rootmse_ate_all <-  array(,dim = c(10,5,27))
rootmse_mean1_all <-  array(,dim = c(10,5,27))
rootmse_hr_all<- array(, dim = c(10,3,27))
for (i in 1:27){
  load(paste0("Simulation results/result_simu_continuous_timevarying_ipw_cali_seq_single_", as.character(i), ".rda"))
  simu.t <- as.data.frame(tidyr::crossing(1:iters, 0:4))
  ate_all <- as.data.frame(list.rbind(oper[14,]))
  mean1_all <- ate_all + 100
  ate_all$simu <- simu.t[,1]
  ate_all$visit <- simu.t[,2]
  mean1_all$simu <- simu.t[,1]
  mean1_all$visit <- simu.t[,2]
  simu.scenario <- as.data.frame(tidyr::crossing(1:iters, 1:10))
  hr_all <- as.data.frame(list.rbind(oper[13,])) %>% 
    dplyr::mutate(simu = simu.scenario[,1], scenario =  simu.scenario[,2])
  for (t in 1:5){
    bias_ate_all[,t,i] <- colMeans(cbind(ate_all[ate_all$visit == t-1,2:6],ate_all[ate_all$visit == t-1,8:12]),na.rm = T) - c(-5-3*(t-1),
                                                                                                                           -5-3*(t-1),
                                                                                                                           -5-3*(t-1),
                                                                                                                           -5-3*(t-1),
                                                                                                                           -5-3*(t-1),
                                                                                                                           -5-3*(t-1),
                                                                                                                           -5-3*(t-1),
                                                                                                                           -5-3*(t-1),
                                                                                                                           -5-3*(t-1),
                                                                                                                           -5-3*(t-1))
    sd_ate_all[,t,i] <-colSds(as.matrix(cbind(ate_all[ate_all$visit == t-1,2:6],ate_all[ate_all$visit == t-1,8:12])),na.rm = T)
    rootmse_ate_all[,t,i]<- sqrt(bias_ate_all[,t,i]^2 + sd_ate_all[,t,i]^2)
    mae_ate_all[,t,i] <- colMeans(abs(cbind(ate_all[ate_all$visit == t-1,2:6],ate_all[ate_all$visit == t-1,8:12])-(-5-3*(t-1))), na.rm = TRUE)
    medae_ate_all[,t,i] <- colMedians(as.matrix(abs(cbind(ate_all[ate_all$visit == t-1,2:6],ate_all[ate_all$visit == t-1,8:12])-(-5-3*(t-1)))), na.rm = TRUE)
    
    bias_mean1_all[,t,i] <- colMeans(cbind(mean1_all[mean1_all$visit == t-1,2:6],mean1_all[mean1_all$visit == t-1,8:12]),na.rm = T) - c(100-5-3*(t-1),
                                                                                                                                        100-5-3*(t-1),
                                                                                                                                        100-5-3*(t-1),
                                                                                                                                        100-5-3*(t-1),
                                                                                                                                        100-5-3*(t-1),
                                                                                                                                        100-5-3*(t-1),
                                                                                                                                        100-5-3*(t-1),
                                                                                                                                        100-5-3*(t-1),
                                                                                                                                        100-5-3*(t-1),
                                                                                                                                        100-5-3*(t-1))
    sd_mean1_all[,t,i] <-colSds(as.matrix(cbind(mean1_all[mean1_all$visit == t-1,2:6],mean1_all[mean1_all$visit == t-1,8:12])),na.rm = T)
    rootmse_mean1_all[,t,i]<- sqrt(bias_mean1_all[,t,i]^2 + sd_mean1_all[,t,i]^2)
    mae_mean1_all[,t,i] <- colMeans(abs(cbind(mean1_all[mean1_all$visit == t-1,2:6],mean1_all[mean1_all$visit == t-1,8:12])-(100-5-3*(t-1))), na.rm = TRUE)
    medae_mean1_all[,t,i] <- colMedians(as.matrix(abs(cbind(mean1_all[mean1_all$visit == t-1,2:6],mean1_all[mean1_all$visit == t-1,8:12])-(100-5-3*(t-1)))), na.rm = TRUE)
  }
  for (k in 1:10){
    bias_hr_all[k,,i]<- colMeans(hr_all[hr_all$scenario == k,1:3], na.rm = TRUE) - c(100,-5,-3)
    sd_hr_all[k,,i]<- colSds(as.matrix(hr_all[hr_all$scenario == k,1:3]), na.rm = TRUE)
    rootmse_hr_all[k,,i]<- sqrt(bias_hr_all[k,,i]^2 + sd_hr_all[k,,i]^2)
    mae_hr_all[k,,i] <- colMeans(abs(hr_all[hr_all$scenario == k,1:3] %>% mutate(V1 = V1 - 100, V2 = V2+5, V3 = V3 + 3)), na.rm = TRUE)
    medae_hr_all[k,,i] <- colMedians(as.matrix(abs(hr_all[hr_all$scenario == k,1:3] %>% mutate(V1 = V1 - 100, V2 = V2+5, V3 = V3 + 3))), na.rm = TRUE)
  }
  
}

bias_hr_table <- bias_hr_all[,,1]
sd_hr_table <- sd_hr_all[,,1]
rootmse_hr_table <- rootmse_hr_all[,,1]
mae_hr_table <- mae_hr_all[,,1]
medae_hr_table <- medae_hr_all[,,1]
for (i in 2:27){
  bias_hr_table <- rbind(bias_hr_table,bias_hr_all[,,i])
  sd_hr_table <- rbind(sd_hr_table,sd_hr_all[,,i])
  rootmse_hr_table <- rbind(rootmse_hr_table,rootmse_hr_all[,,i])
  mae_hr_table <- rbind(mae_hr_table, mae_hr_all[,,i])
  medae_hr_table <- rbind(medae_hr_table, medae_hr_all[,,i])
}

colnames(bias_hr_table) <- c('Bias 1', 'Bias 2', 'Bias 3')
colnames(sd_hr_table) <- c('SD 1', 'SD 2', 'SD 3')
colnames(rootmse_hr_table) <- c('rootMSE 1', 'rootMSE 2', 'rootMSE 3')
colnames(mae_hr_table) <- c('MAE 1', 'MAE 2', 'MAE 3')
colnames(medae_hr_table) <- c('MedAE 1', 'MedAE 2', 'MedAE 3')
size <- c(200,500,1000)
treat <- c(0.25,0.5,0.75)
conf <- c(1,1.5,2)
var <- c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
         'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
scenarios <- tidyr::crossing(size, treat,conf,var) %>% 
  dplyr::arrange(size, treat, conf,
                 match(var, c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
                              'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
                 ))

na_list <- c(6)
for (r in 2:54){
  na_list <- cbind(na_list, r*5 + r)
}
hr_table <- cbind(scenarios,bias_hr_table, sd_hr_table,rootmse_hr_table, mae_hr_table,medae_hr_table) %>% 
  dplyr::select(size, treat, conf,var, 'Bias 1', 'SD 1', 'rootMSE 1', 'MAE 1','MedAE 1','Bias 2', 'SD 2', 'rootMSE 2','MAE 2','MedAE 2','Bias 3', 'SD 3', 'rootMSE 3', 'MAE 3','MedAE 3')%>% 
  insertRows(r = na_list, new = NA)


print(xtable(hr_table[hr_table$size == 200,],digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')

print(xtable(hr_table[hr_table$size == 500,],digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')

print(xtable(hr_table[hr_table$size == 1000,],digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),include.rownames=FALSE, type = 'latex')

write_csv(hr_table, file = 'simulation_results_hr_table.csv')

bias_ate_table <- bias_ate_all[,,1]
sd_ate_table <- sd_ate_all[,,1]
rootmse_ate_table <- rootmse_ate_all[,,1]
mae_ate_table <- mae_ate_all[,,1]
medae_ate_table <- medae_ate_all[,,1]
for (i in 2:27){
  bias_ate_table <- rbind(bias_ate_table,bias_ate_all[,,i])
  sd_ate_table <- rbind(sd_ate_table,sd_ate_all[,,i])
  rootmse_ate_table <- rbind(rootmse_ate_table,rootmse_ate_all[,,i])
  mae_ate_table <- rbind(mae_ate_table, mae_ate_all[,,i])
  medae_ate_table <- rbind(medae_ate_table, medae_ate_all[,,i])
}

colnames(bias_ate_table) <- c('Bias 0', 'Bias 1','Bias 2', 'Bias 3','Bias 4')
colnames(sd_ate_table) <- c('SD 0', 'SD 1','SD 2', 'SD 3','SD 4')
colnames(rootmse_ate_table) <- c('rootMSE 0', 'rootMSE 1','rootMSE 2', 'rootMSE 3','rootMSE 4')
colnames(mae_ate_table) <- c('MAE 0','MAE 1', 'MAE 2', 'MAE 3','MAE 4')
colnames(medae_ate_table) <- c('MedAE 0','MedAE 1', 'MedAE 2', 'MedAE 3','MedAE 4')
size <- c(200,500,1000)
treat <- c(0.25,0.50,0.75)
conf <- c(1,1.5,2)
var <- c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
         'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
scenarios <- tidyr::crossing(size, treat,conf,var)%>% 
  dplyr::arrange(size, treat, conf,
                 match(var, c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
                              'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')))


ate_table <- cbind(scenarios,bias_ate_table, sd_ate_table,rootmse_ate_table,mae_ate_table,medae_ate_table) %>% 
  dplyr::select(size, treat, conf,var,
                'Bias 0', 'SD 0', 'rootMSE 0', 'MAE 0','MedAE 0',
                'Bias 1', 'SD 1', 'rootMSE 1', 'MAE 1','MedAE 1',
                'Bias 2', 'SD 2', 'rootMSE 2','MAE 2','MedAE 2',
                'Bias 3', 'SD 3', 'rootMSE 3', 'MAE 3','MedAE 3',
                'Bias 4', 'SD 4', 'rootMSE 4', 'MAE 4','MedAE 4')%>% 
  insertRows(r = na_list, new = NA)

write_csv(ate_table, file = 'simulation_results_ate_table.csv')
print(xtable(ate_table[ate_table$size == 200,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
             include.rownames=FALSE, type = 'latex')
print(xtable(ate_table[ate_table$size == 500,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')
print(xtable(ate_table[ate_table$size == 1000,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')


bias_mean1_table <- bias_mean1_all[,,1]
sd_mean1_table <- sd_mean1_all[,,1]
rootmse_mean1_table <- rootmse_mean1_all[,,1]
mae_mean1_table <- mae_mean1_all[,,1]
medae_mean1_table <- medae_mean1_all[,,1]
for (i in 2:27){
  bias_mean1_table <- rbind(bias_mean1_table,bias_mean1_all[,,i])
  sd_mean1_table <- rbind(sd_mean1_table,sd_mean1_all[,,i])
  rootmse_mean1_table <- rbind(rootmse_mean1_table,rootmse_mean1_all[,,i])
  mae_mean1_table <- rbind(mae_mean1_table, mae_mean1_all[,,i])
  medae_mean1_table <- rbind(medae_mean1_table, medae_mean1_all[,,i])
}

colnames(bias_mean1_table) <- c('Bias 0', 'Bias 1','Bias 2', 'Bias 3','Bias 4')
colnames(sd_mean1_table) <- c('SD 0', 'SD 1','SD 2', 'SD 3','SD 4')
colnames(rootmse_mean1_table) <- c('rootMSE 0', 'rootMSE 1','rootMSE 2', 'rootMSE 3','rootMSE 4')
colnames(mae_mean1_table) <- c('MAE 0','MAE 1', 'MAE 2', 'MAE 3','MAE 4')
colnames(medae_mean1_table) <- c('MedAE 0','MedAE 1', 'MedAE 2', 'MedAE 3','MedAE 4')
size <- c(200,500,1000)
treat <- c(0.25,0.50,0.75)
conf <- c(1,1.5,2)
var <- c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
         'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')
scenarios <- tidyr::crossing(size, treat,conf,var)%>% 
  dplyr::arrange(size, treat, conf,
                 match(var, c('Naive', 'IPW', 'Aggregated cali.', 'Calibration by time','Aggregated cali. with time', 
                              'Naive miss.', 'IPW miss.','Aggregated cali. miss.', 'Calibration by time miss.', 'Aggregated cali. with time miss.')))


mean1_table <- cbind(scenarios,bias_mean1_table, sd_mean1_table,rootmse_mean1_table,mae_mean1_table,medae_mean1_table) %>% 
  dplyr::select(size, treat, conf,var,
                'Bias 0', 'SD 0', 'rootMSE 0', 'MAE 0','MedAE 0',
                'Bias 1', 'SD 1', 'rootMSE 1', 'MAE 1','MedAE 1',
                'Bias 2', 'SD 2', 'rootMSE 2','MAE 2','MedAE 2',
                'Bias 3', 'SD 3', 'rootMSE 3', 'MAE 3','MedAE 3',
                'Bias 4', 'SD 4', 'rootMSE 4', 'MAE 4','MedAE 4')%>% 
  insertRows(r = na_list, new = NA)

write_csv(mean1_table, file = 'simulation_results_mean1_table.csv')
print(xtable(mean1_table[mean1_table$size == 200,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')
print(xtable(mean1_table[mean1_table$size == 500,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')
print(xtable(mean1_table[mean1_table$size == 1000,],
             digits=c(0,0,2,1,0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)),
      include.rownames=FALSE, type = 'latex')
