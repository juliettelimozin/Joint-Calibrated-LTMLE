library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
load('Simulation results/true_HR_singletrial_scenario2.rda')
load('Simulation results/true_MRD_singletrial_scenario2.rda')
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

treat_pos <- c(-1,-0.8,-0.5,-0.2,0,0.2,0.5,0.8,1)
outcomes <- c("low", 'med', 'high')


size <- c(500,1000,5000)
treat <- c(-1,0,1)

scenarios <- as.data.frame(tidyr::crossing(size, treat))

iters <- 200

bias_mrd <- array(,dim = c(4,5,9))
bias_treated <- array(,dim = c(4,5,9))
bias_control <- array(,dim = c(4,5,9))
point_bias_mrd<- array(,dim = c(4,5,9))
bias_hr <- array(,dim = c(4,9))
sd_mrd <- array(, dim = c(4,5,9))
sd_hr <- array(,dim = c(4,9))
mr_all <- array(,dim = c(4,2,5,iters,9))
meandiffs_all <- array(, dim = c(5,5,6,iters,9))
objectives_all <- array(0,dim = c(4,8,iters,9))


for (i in 1:1){
  load(paste0("Simulation results/simulation_ipw_cali_cali_seq_singletrial", i, ".rda"))
  simu.t <- as.data.frame(tidyr::crossing(1:iters, 0:4))
  mr_all <- list.rbind(oper[10,]) %>% 
    dplyr::mutate(simu =simu.t[,1], visit =  simu.t[,2])
  scenario <- i%%3
  if (scenario ==0){scenario <- 3}
  for (t in 1:5)
  bias_mrd[,t,i] <- rowMeans(mr_all[1,1,,,i] - mr_all[1,2,,,i], na.rm = T) - (true_MRD[,1,scenario,2] - true_MRD[,2,scenario,2])
  bias_treated[,t,i] <- rowMeans(mr_all[1,2,,,i], na.rm = T) - (true_MRD[,2,scenario,2])
  bias_control[1,,i] <- rowMeans(mr_all[1,1,,,i], na.rm = T) - (true_MRD[,1,scenario,2])
  bias_control[2,,i] <- rowMeans(mr_all[2,1,,,i] , na.rm = T) - (true_MRD[,1,scenario,2])
  bias_control[3,,i] <- rowMeans(mr_all[3,1,,,i], na.rm = T) - (true_MRD[,1,scenario,2])
  bias_control[4,,i] <- rowMeans(mr_all[4,1,,,i] , na.rm = T) - (true_MRD[,1,scenario,2])
  point_bias_mrd[1,,i] <- rowMeans(mr_all[1,1,,,i] - mr_all[1,2,,,i], na.rm = T) / (true_MRD[,1,scenario,2] - true_MRD[,2,scenario,2])
  point_bias_mrd[2,,i] <- rowMeans(mr_all[2,1,,,i] - mr_all[2,2,,,i], na.rm = T) / (true_MRD[,1,scenario,2] - true_MRD[,2,scenario,2])
  point_bias_mrd[3,,i] <- rowMeans(mr_all[3,1,,,i] - mr_all[3,2,,,i], na.rm = T) / (true_MRD[,1,scenario,2] - true_MRD[,2,scenario,2])
  point_bias_mrd[4,,i] <- rowMeans(mr_all[4,1,,,i] - mr_all[4,2,,,i], na.rm = T) / (true_MRD[,1,scenario,2] - true_MRD[,2,scenario,2])
  sd_mrd[1,,i] <- rowSds(mr_all[1,1,,,i] - mr_all[1,2,,,i], na.rm = T) 
  sd_mrd[2,,i] <- rowSds(mr_all[2,1,,,i] - mr_all[2,2,,,i], na.rm = T) 
  sd_mrd[3,,i] <- rowSds(mr_all[3,1,,,i] - mr_all[3,2,,,i], na.rm = T) 
  sd_mrd[4,,i] <- rowSds(mr_all[4,1,,,i] - mr_all[4,2,,,i], na.rm = T) 
  bias_hr[,i]<- rowMeans(hr_estimates, na.rm = T) - c(true_HR[scenario,1],true_HR[scenario,1])
  sd_hr[,i]<- rowSds(hr_estimates, na.rm = T) 
}

# mr_plot_low <- lapply(1:9, function(i){
#   scenario <- i%%9
#   if (scenario ==0){scenario <- 9}
#   plot <- ggplot() +
#     scale_color_manual(name = "Weight type", 
#                        values = c("True: never treated"= "black", "True: always treated" = "black",
#                                   'Estimated: never treated' = 'red', 'Estimated: always treated' = 'blue')) +
#     labs(x = paste0('N = ', scenarios[i,1], ',\nMisspecification = ', scenarios[i,2], ', \nTreat. prev. = ',scenarios[i,3]),
#          y = "Empirical SD of MRD estimation") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
#     ylim(0.5,1)
#   cmat0 = t(mr_all[1,1,,,i])
#   rownames(cmat0) = paste("simu", seq(1000), sep="")
#   colnames(cmat0) = paste("time", seq(5), sep="")
#   dat0 = as.data.frame(cmat0)
#   dat0$simu = rownames(dat0)
#   mdat0 = melt(dat0, id.vars="simu")
#   mdat0$time = as.numeric(gsub("time", "", mdat0$variable))
#   
#   cmat1 = t(mr_all[1,2,,,i])
#   rownames(cmat1) = paste("simu", seq(1000), sep="")
#   colnames(cmat1) = paste("time", seq(5), sep="")
#   dat1 = as.data.frame(cmat1)
#   dat1$simu = rownames(dat1)
#   mdat1 = melt(dat1, id.vars="simu")
#   mdat1$time = as.numeric(gsub("time", "", mdat1$variable))
#   
#   plot+
#     geom_line(aes(x=mdat0$time-1, y=mdat0$value, group=mdat0$simu, colour = 'Estimated: never treated'),
#               size=0.1, alpha=0.05 ) +
#     geom_line(aes(x=mdat1$time-1, y=mdat1$value, group=mdat1$simu, colour = 'Estimated: always treated'),
#               size=0.1, alpha=0.05) +
#     geom_line(aes(x = 0:4, y = true_MRD[,1,scenario,2], colour = 'True: never treated'), linetype = 1,size=0.5) +
#     geom_line(aes(x = 0:4, y = true_MRD[,2,scenario,2], colour = 'True: always treated'), linetype = 2,size=0.5) 
#   
# })
# annotate_figure(ggarrange(plotlist = mr_plot_low[1:9], nrow = 3, ncol = 9,common.legend = T , legend = 'bottom'),top = 'Low event rate')

################BIAS, SD, MSE PLOTS ###################
bias_plots_med_mrd <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_mrd[1,,i], colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = bias_mrd[1,,i],colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = bias_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = bias_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_line(aes(x = 0:4, y = bias_mrd[3,,i], colour = 'MLE-IPW mis.')) +
    geom_point(aes(x = 0:4, y = bias_mrd[3,,i],colour = 'MLE-IPW mis.')) +
    geom_line(aes(x = 0:4, y = bias_mrd[4,,i], colour = 'Calibrated weights mis.')) +
    geom_point(aes(x = 0:4, y = bias_mrd[4,,i], colour = 'Calibrated weights mis.')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue",
                                                        'MLE-IPW mis.' = 'purple', 'Calibrated weights mis.' = 'green')) +
    labs(x = paste0('N = ', scenarios[i,1], ', \nTreat. prev. = ',scenarios[i,2]),
         y = "Bias") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(-0.03,0.02)
})
annotate_figure(ggarrange(plotlist = bias_plots_med_mrd[1:9], nrow = 3, ncol = 3,common.legend = T , legend = 'bottom'),top = 'MRD bias, Medium event rate')

bias_plots_med_mrd <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = point_bias_mrd[1,,i], colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = point_bias_mrd[1,,i],colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = point_bias_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = point_bias_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_line(aes(x = 0:4, y = point_bias_mrd[3,,i], colour = 'MLE-IPW mis.')) +
    geom_point(aes(x = 0:4, y = point_bias_mrd[3,,i],colour = 'MLE-IPW mis.')) +
    geom_line(aes(x = 0:4, y = point_bias_mrd[4,,i], colour = 'Calibrated weights mis.')) +
    geom_point(aes(x = 0:4, y = point_bias_mrd[4,,i], colour = 'Calibrated weights mis.')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue",
                                                        'MLE-IPW mis.' = 'purple', 'Calibrated weights mis.' = 'green')) +
    labs(x = paste0('N = ', scenarios[i,1], ', \nTreat. prev. = ',scenarios[i,2]),
         y = "Point Bias") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(0.6,1.4)
})
annotate_figure(ggarrange(plotlist = bias_plots_med_mrd[1:9], nrow = 3, ncol = 3,common.legend = T , legend = 'bottom'),top = 'MRD point-bias, Medium event rate')

bias_plots_med_mrd <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_treated[1,,i], colour = 'MLE-IPW', linetype = 'Always treated')) +
    geom_point(aes(x = 0:4, y = bias_treated[1,,i],colour = 'MLE-IPW', linetype = 'Always treated')) +
    geom_line(aes(x = 0:4, y = bias_treated[2,,i], colour = 'Calibrated weights', linetype = 'Always treated')) +
    geom_point(aes(x = 0:4, y = bias_treated[2,,i], colour = 'Calibrated weights', linetype = 'Always treated')) +
    geom_line(aes(x = 0:4, y = bias_treated[3,,i], colour = 'MLE-IPW mis.', linetype = 'Always treated')) +
    geom_point(aes(x = 0:4, y = bias_treated[3,,i],colour = 'MLE-IPW mis.', linetype = 'Always treated')) +
    geom_line(aes(x = 0:4, y = bias_treated[4,,i], colour = 'Calibrated weights mis.', linetype = 'Always treated')) +
    geom_point(aes(x = 0:4, y = bias_treated[4,,i], colour = 'Calibrated weights mis.', linetype = 'Always treated')) +
    geom_line(aes(x = 0:4, y = bias_control[1,,i], colour = 'MLE-IPW', linetype = 'Never treated')) +
    geom_point(aes(x = 0:4, y = bias_control[1,,i],colour = 'MLE-IPW', linetype = 'Never treated')) +
    geom_line(aes(x = 0:4, y = bias_control[2,,i], colour = 'Calibrated weights', linetype = 'Never treated')) +
    geom_point(aes(x = 0:4, y = bias_control[2,,i], colour = 'Calibrated weights', linetype = 'Never treated')) +
    geom_line(aes(x = 0:4, y = bias_control[3,,i], colour = 'MLE-IPW mis.', linetype = 'Never treated')) +
    geom_point(aes(x = 0:4, y = bias_control[3,,i],colour = 'MLE-IPW mis.', linetype = 'Never treated')) +
    geom_line(aes(x = 0:4, y = bias_control[4,,i], colour = 'Calibrated weights mis.', linetype = 'Never treated')) +
    geom_point(aes(x = 0:4, y = bias_control[4,,i], colour = 'Calibrated weights mis.', linetype = 'Never treated')) +
    scale_linetype_manual(name = 'Treatment protocol', values = c("Always treated" = 1, 'Never treated' = 2)) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue",
                                                        'MLE-IPW mis.' = 'purple', 'Calibrated weights mis.' = 'green')) +
    labs(x = paste0('N = ', scenarios[i,1], ', \nTreat. prev. = ',scenarios[i,2]),
         y = "Point Bias") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +
    ylim(-0.04,0.01)
  
})
annotate_figure(ggarrange(plotlist = bias_plots_med_mrd[1:9], nrow = 3, ncol = 3,common.legend = T , legend = 'bottom'),top = 'MRD point-bias, Medium event rate')

sd_plots_med <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = sd_mrd[1,,i], colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = sd_mrd[1,,i],colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = sd_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = sd_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_line(aes(x = 0:4, y = sd_mrd[3,,i], colour = 'MLE-IPW mis.')) +
    geom_point(aes(x = 0:4, y = sd_mrd[3,,i],colour = 'MLE-IPW mis.')) +
    geom_line(aes(x = 0:4, y = sd_mrd[4,,i], colour = 'Calibrated weights mis.')) +
    geom_point(aes(x = 0:4, y = sd_mrd[4,,i], colour = 'Calibrated weights mis.')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue",
                                                        'MLE-IPW mis.' = 'purple', 'Calibrated weights mis.' = 'green')) +
    labs(x = paste0('N = ', scenarios[i,1], ', \nTreat. prev. = ',scenarios[i,2]),
         y = "SD") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(0,0.2)
})
annotate_figure(ggarrange(plotlist = sd_plots_med[1:9], nrow = 3, ncol = 3,common.legend = T , legend = 'bottom'), top = 'MRD SD, medium event rate')

bias.sd_plots_med <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_mrd[1,,i]/sd_mrd[1,,i], colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = bias_mrd[1,,i]/sd_mrd[1,,i],colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = bias_mrd[2,,i]/sd_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = bias_mrd[2,,i]/sd_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_line(aes(x = 0:4, y = bias_mrd[3,,i]/sd_mrd[3,,i], colour = 'MLE-IPW mis.')) +
    geom_point(aes(x = 0:4, y = bias_mrd[3,,i]/sd_mrd[3,,i],colour = 'MLE-IPW mis.')) +
    geom_line(aes(x = 0:4, y = bias_mrd[4,,i]/sd_mrd[4,,i], colour = 'Calibrated weights mis.')) +
    geom_point(aes(x = 0:4, y = bias_mrd[4,,i]/sd_mrd[4,,i], colour = 'Calibrated weights mis.')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue",
                                                        'MLE-IPW mis.' = 'purple', 'Calibrated weights mis.' = 'green')) +
    labs(x = paste0('N = ', scenarios[i,1], ', \nTreat. prev. = ',scenarios[i,2]),
         y = "SD") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(-0.5,0.5)
})
annotate_figure(ggarrange(plotlist = bias.sd_plots_med[1:9], nrow = 3, ncol = 3,common.legend = T , legend = 'bottom'), top = 'MRD SD, medium event rate')


mse_plots_low <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_mrd[1,,i]^2 + sd_mrd[1,,i]^2, colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = bias_mrd[1,,i]^2 + sd_mrd[1,,i]^2,colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = bias_mrd[2,,i]^2 + sd_mrd[2,,i]^2, colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = bias_mrd[2,,i]^2 + sd_mrd[2,,i]^2, colour = 'Calibrated weights')) +
    geom_line(aes(x = 0:4, y = bias_mrd[3,,i]^2 + sd_mrd[3,,i]^2, colour = 'MLE-IPW mis.')) +
    geom_point(aes(x = 0:4, y = bias_mrd[3,,i]^2 + sd_mrd[3,,i]^2,colour = 'MLE-IPW mis.')) +
    geom_line(aes(x = 0:4, y = bias_mrd[4,,i]^2 + sd_mrd[4,,i]^2, colour = 'Calibrated weights mis.')) +
    geom_point(aes(x = 0:4, y = bias_mrd[4,,i]^2 + sd_mrd[4,,i]^2, colour = 'Calibrated weights mis.')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue",
                                                        'MLE-IPW mis.' = 'purple', 'Calibrated weights mis.' = 'green')) +
    labs(x = paste0('N = ', scenarios[i,1], ', \nTreat. prev. = ',scenarios[i,2]),
         y = "Empirical MSE") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(0,0.1)
})
annotate_figure(ggarrange(plotlist = mse_plots_low[1:9], nrow = 3, ncol = 3,common.legend = T , legend = 'bottom'), top = 'MRD MSE, Medium event rate')

randomiter <- sample(1:200,1)

meandiffsX1_treated_low <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = meandiffs_all[1,,1,randomiter,i], colour = 'Unadjusted')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[1,,1,randomiter,i], colour = 'Unadjusted')) +
    geom_line(aes(x = 0:4, y = meandiffs_all[2,,1,randomiter,i], colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[2,,1,randomiter,i], colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = meandiffs_all[3,,1,randomiter,i], colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[3,,1,randomiter,i], colour = 'Calibrated weights')) +
    geom_line(aes(x = 0:4, y = meandiffs_all[4,,1,randomiter,i], colour = 'MLE-IPW mis.')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[4,,1,randomiter,i], colour = 'MLE-IPW mis.')) +
    geom_line(aes(x = 0:4, y = meandiffs_all[5,,1,randomiter,i], colour = 'Calibrated weights mis.')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[5,,1,randomiter,i], colour = 'Calibrated weights mis.')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue", 'Unadjusted' = 'grey',
                                                        'MLE-IPW mis.' = 'purple', 'Calibrated weights mis.' = 'green')) +
    labs(x = paste0('N = ', scenarios[i,1], ', \nTreat. prev. = ',scenarios[i,2]),
         y = "Mean") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    
    geom_hline(yintercept = meandiffs_all[1,1,1,randomiter,i],linetype = 'dashed')
})
annotate_figure(ggarrange(plotlist = meandiffsX1_treated_low[1:9], nrow = 3, ncol = 3,common.legend = T , legend = 'bottom'), top = 'Mean of X1 in treated')

meandiffsX1_untreated_low <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = meandiffs_all[1,,4,randomiter,i], colour = 'Unadjusted')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[1,,4,randomiter,i], colour = 'Unadjusted')) +
    geom_line(aes(x = 0:4, y = meandiffs_all[2,,4,randomiter,i], colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[2,,4,randomiter,i], colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = meandiffs_all[3,,4,randomiter,i], colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[3,,4,randomiter,i], colour = 'Calibrated weights')) +
    geom_line(aes(x = 0:4, y = meandiffs_all[4,,4,randomiter,i], colour = 'MLE-IPW mis.')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[4,,4,randomiter,i], colour = 'MLE-IPW mis.')) +
    geom_line(aes(x = 0:4, y = meandiffs_all[5,,4,randomiter,i], colour = 'Calibrated weights mis.')) +
    geom_point(aes(x = 0:4, y = meandiffs_all[5,,4,randomiter,i], colour = 'Calibrated weights mis.')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue", 'Unadjusted' = 'grey',
                                                        'MLE-IPW mis.' = 'purple', 'Calibrated weights mis.' = 'green')) +
    labs(x = paste0('N = ', scenarios[i,1], ', \nTreat. prev. = ',scenarios[i,2]),
         y = "Mean") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    geom_hline(yintercept = 0,linetype = 'dashed')
  
})
annotate_figure(ggarrange(plotlist = meandiffsX1_untreated_low[1:9], nrow = 3, ncol = 3,common.legend = T , legend = 'bottom'), top = 'Mean of X1 in untreated')

objectives_low <- lapply(1:9, function(i){
  ggplot() +
    geom_point(aes(x = objectives_all[1,,randomiter,i], 
                   y = 1:8, colour = 'MLE-IPW')) +
    geom_point(aes(x = objectives_all[2,,randomiter,i], 
                   y = 1:8, colour = "Calibrated weights")) +
    geom_point(aes(x = objectives_all[3,,randomiter,i], 
                   y = 1:8, colour = 'MLE-IPW')) +
    geom_point(aes(x = objectives_all[4,,randomiter,i], 
                   y = 1:8, colour = "Calibrated weights")) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue", 'Unadjusted' = 'grey',
                                                        'MLE-IPW mis.' = 'purple', 'Calibrated weights mis.' = 'green')) +
    labs(x = paste0('N = ', scenarios[i,1], ',\nMisspecification = ', scenarios[i,2], ', \nTreat. prev. = ',scenarios[i,3]),
         y = "Restriction") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    scale_y_discrete(limits = c('A1', 'A1X2','A1X2sq', 'A1nextX1','A0','A0X2','A0X2sq','A0nextX1')) +
    geom_hline(yintercept = 0,linetype = 'dashed')
  
})
annotate_figure(ggarrange(plotlist = objectives_low[1:9], nrow = 3, ncol = 9,common.legend = T , legend = 'bottom'), top = 'Calibration restriction objectives')



