library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
load("/Users/juliette/Documents/MPhil PHS 21-22/Multiple-trial-emulation-IPTW-MSM-CIs/Code/HPC output/true_value_red_newsimus.rda")
load('true_HR_singletrial.rda')
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

absmeandiffs_all <- array(,dim = c(3,2,4,500,27))
objectives_all <- array(,dim = c(2,6,500,27))

treat_pos <- c(-1,-0.8,-0.5,-0.2,0,0.2,0.5,0.8,1)
outcomes <- c("low", 'med', 'high')


size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(1,3,5)

scenarios <- as.data.frame(tidyr::crossing(size,conf, treat))
bias_hr <- array(,dim = c(2,27))
sd_hr <- array(,dim = c(2,27))
bias_mrd <- array(,dim = c(2,5,27))
sd_mrd <- array(,dim = c(2,5,27))

for (i in 1:19){
    load(paste0("absmeandiffs_singletrial_", i, ".rda"))
    load(paste0("objectives_singletrial_", i, ".rda"))
    load(paste0("hr_estimates_singletrial_", i, ".rda"))
    load(paste0("mrd_estimates_singletrial_", i, ".rda"))
    absmeandiffs_all[,,,,i] <- absmeandiffs
    objectives_all[,,,i] <- objectives
    scenario <- i%%9
    if (scenario ==0){scenario <- 9}
    bias_mrd[1,,i] <- rowMeans(mrd_estimates[1,,], na.rm = TRUE) - true_value_red[,scenario,1]
    bias_mrd[2,,i] <- rowMeans(mrd_estimates[2,,], na.rm = TRUE) - true_value_red[,scenario,1]
    sd_mrd[1,,i] <- rowSds(mrd_estimates[1,,], na.rm = TRUE)
    sd_mrd[2,,i] <- rowSds(mrd_estimates[2,,], na.rm = TRUE)
    bias_hr[,i]<- rowMeans(hr_estimates, na.rm = T) - c(true_HR[scenario,1],true_HR[scenario,1])
    sd_hr[,i]<- rowSds(hr_estimates, na.rm = T) 
}

colnames(mean_time) <- c('Outcome_prevalence', 'Sample_size', 'Confounding', 'Treatment_prevalence', 'Bootstrap', 'LEF_outcome',
                         'LEF_both', 'Sandwich')

################BIAS, SD, MSE PLOTS ###################
bias_plots_low_mrd <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_mrd[1,,i], colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = bias_mrd[1,,i],colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = bias_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = bias_mrd[2,,i], colour = 'Calibrated weights')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue")) +
    labs(x = paste0('N = ', scenarios[i,1], ',\nConfounding = ', scenarios[i,2], ', \nTreat. prev. = ',scenarios[i,3]),
         y = "Empirical bias of MRD estimation") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(-0.1,0.5)
})
annotate_figure(ggarrange(plotlist = bias_plots_low_mrd[1:27], nrow = 3, ncol = 9,common.legend = T , legend = 'bottom'),top = 'Low event rate')

sd_plots_low <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = sd_mrd[1,,i], colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = sd_mrd[1,,i],colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = sd_mrd[2,,i], colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = sd_mrd[2,,i], colour = 'Calibrated weights')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue")) +
    labs(x = paste0('N = ', scenarios[i,1], ',\nConfounding = ', scenarios[i,2], ', \nTreat. prev. = ',scenarios[i,3]),
         y = "Empirical SD of MRD estimation") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(-0.1,0.5)
})
annotate_figure(ggarrange(plotlist = sd_plots_low[1:27], nrow = 3, ncol = 9,common.legend = T , legend = 'bottom'), top = 'Low event rate')


mse_plots_low <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_mrd[1,,i]^2 + sd_mrd[1,,i]^2, colour = 'MLE-IPW')) +
    geom_point(aes(x = 0:4, y = bias_mrd[1,,i]^2 + sd_mrd[1,,i]^2,colour = 'MLE-IPW')) +
    geom_line(aes(x = 0:4, y = bias_mrd[2,,i]^2 + sd_mrd[2,,i]^2, colour = 'Calibrated weights')) +
    geom_point(aes(x = 0:4, y = bias_mrd[2,,i]^2 + sd_mrd[2,,i]^2, colour = 'Calibrated weights')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue")) +
    labs(x = paste0('N = ', scenarios[i,1], ',\nConfounding = ', scenarios[i,2], ', \nTreat. prev. = ',scenarios[i,3]),
         y = "Empirical MSE of MRD estimation") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(-0.1,0.5)
})
annotate_figure(ggarrange(plotlist = mse_plots_low[1:27], nrow = 3, ncol = 9,common.legend = T , legend = 'bottom'), top = 'Low event rate')

absmeandiffsX2_low <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 1:4, y = rowMeans(absmeandiffs_all[1,1,,,i], na.rm = T), colour = 'Unadjusted')) +
    geom_point(aes(x = 1:4, y = rowMeans(absmeandiffs_all[1,1,,,i], na.rm = T), colour = 'Unadjusted')) +
    geom_line(aes(x = 1:4, y = rowMeans(absmeandiffs_all[2,1,,,i], na.rm = T), colour = 'MLE-IPW')) +
    geom_point(aes(x = 1:4, y = rowMeans(absmeandiffs_all[2,1,,,i], na.rm = T), colour = 'MLE-IPW')) +
    geom_line(aes(x = 1:4, y = rowMeans(absmeandiffs_all[3,1,,,i], na.rm = T), colour = 'Calibrated weights')) +
    geom_point(aes(x = 1:4, y = rowMeans(absmeandiffs_all[3,1,,,i], na.rm = T), colour = 'Calibrated weights')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue", 'Unadjusted' = 'grey')) +
    labs(x = paste0('N = ', scenarios[i,1], ',\nConfounding = ', scenarios[i,2], ', \nTreat. prev. = ',scenarios[i,3]),
         y = "ASMD of X2") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(0,3)
})
annotate_figure(ggarrange(plotlist = absmeandiffsX2_low[1:27], nrow = 3, ncol = 9,common.legend = T , legend = 'bottom'), top = 'Low event rate')

absmeandiffsX4_low <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 1:4, y = rowMeans(absmeandiffs_all[1,2,,,i], na.rm = T), colour = 'Unadjusted')) +
    geom_point(aes(x = 1:4, y = rowMeans(absmeandiffs_all[1,2,,,i], na.rm = T), colour = 'Unadjusted')) +
    geom_line(aes(x = 1:4, y = rowMeans(absmeandiffs_all[2,2,,,i], na.rm = T), colour = 'MLE-IPW')) +
    geom_point(aes(x = 1:4, y = rowMeans(absmeandiffs_all[2,2,,,i], na.rm = T), colour = 'MLE-IPW')) +
    geom_line(aes(x = 1:4, y = rowMeans(absmeandiffs_all[3,2,,,i], na.rm = T), colour = 'Calibrated weights')) +
    geom_point(aes(x = 1:4, y = rowMeans(absmeandiffs_all[3,2,,,i], na.rm = T), colour = 'Calibrated weights')) +
    scale_color_manual(name = "Weight type", values = c("MLE-IPW"= "red", "Calibrated weights" = "blue", 'Unadjusted' = 'grey')) +
    labs(x = paste0('N = ', scenarios[i,1], ',\nConfounding = ', scenarios[i,2], ', \nTreat. prev. = ',scenarios[i,3]),
         y = "ASMD of X4") + theme(aspect.ratio = 1, axis.title = element_text(size = 10)) +  
    ylim(0,3)
})
annotate_figure(ggarrange(plotlist = absmeandiffsX4_low[1:27], nrow = 3, ncol = 9,common.legend = T , legend = 'bottom'), top = 'Low event rate')


