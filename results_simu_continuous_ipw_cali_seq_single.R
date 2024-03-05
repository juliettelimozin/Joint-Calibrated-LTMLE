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

size <- c(500,1000,5000)
treat <- c(-1,0,1)

scenarios <- as.data.frame(tidyr::crossing(size, treat))

iters <- 200

bias_mr_all <- array(,dim = c(3*8,5,3))
bias_hr_all <- array(, dim = c(8,3,3))

for (i in 4:6){
  load(paste0("Simulation results/result_simu_continuous_ipw_cali_seq_single_", i, ".rda"))
  simu.t <- as.data.frame(tidyr::crossing(1:iters, 0:4))
  mr_all <- list.rbind(oper[10,]) %>% 
    dplyr::mutate(simu =simu.t[,1], visit =  simu.t[,2])
  scenario <- i%%3
  simu.scenario <- as.data.frame(tidyr::crossing(1:iters, 1:8))
  hr_all <- as.data.frame(list.rbind(oper[9,])) %>% 
    dplyr::mutate(simu = simu.scenario[,1], scenario =  simu.scenario[,2])
  for (t in 1:5)
    if (t == 1){
     bias_mr_all[,t,i-3] <- colMeans(mr_all[mr_all$visit == t-1,1:24],na.rm = T) - c(95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5,
                                                                               95,100,-5)
    } else {
      bias_mr_all[,t,i-3] <- colMeans(mr_all[mr_all$visit == t-1,1:24],na.rm = T) - c(100-5*t-3,100,-5*t-3,
                                                                                100-5*t-3,100,-5*t-3,
                                                                                100-5*t-3,100,-5*t-3,
                                                                                100-5*t-3,100,-5*t-3,
                                                                                100-5*t-3,100,-5*t-3,
                                                                                100-5*t-3,100,-5*t-3,
                                                                                100-5*t-3,100,-5*t-3,
                                                                                100-5*t-3,100,-5*t-3)
    }
    for (k in 1:8){
      bias_hr_all[k,,i-3]<- colMeans(hr_all[hr_all$scenario == k, 1:3], na.rm = TRUE) - c(100,-5,-3)
    }
}

print(xtable(bias_hr_all[,,1]),include.rownames=FALSE, type = 'latex')
print(xtable(bias_hr_all[,,2]),include.rownames=FALSE, type = 'latex')
print(xtable(bias_hr_all[,,3]),include.rownames=FALSE, type = 'latex')

print(xtable(bias_mr_all[c(3,6,9,12,15,18,21,24),,1]),include.rownames=FALSE, type = 'latex')
print(xtable(bias_mr_all[c(3,6,9,12,15,18,21,24),,2]),include.rownames=FALSE, type = 'latex')
print(xtable(bias_mr_all[c(3,6,9,12,15,18,21,24),,3]),include.rownames=FALSE, type = 'latex')

print(xtable(bias_mr_all[c(3,6,9,12,15,18,21,24),,1]),include.rownames=FALSE, type = 'latex')
print(xtable(bias_mr_all[c(3,6,9,12,15,18,21,24),,2]),include.rownames=FALSE, type = 'latex')
print(xtable(bias_mr_all[c(3,6,9,12,15,18,21,24),,3]),include.rownames=FALSE, type = 'latex')

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



