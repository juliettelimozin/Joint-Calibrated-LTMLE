library(modelr)
library(tidyverse)
library(tidyr)
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("simulate_MSM_treatment_switch.R")
set.seed(NULL)
library(MASS)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)


treat_prev <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)
sample_size <- c(200,1000,5000)
alpha_y <- c(-4.7,-3.8,-3)
summary_table <- array(,dim = c(27,5))

scenarios <- tidyr::crossing(sample_size,conf,alpha_y)

for (i in 1:27){
  data <- DATA_GEN_treatment_switch(as.numeric(scenarios[i,1]), 5, conf = as.numeric(scenarios[i,2]), treat_prev = 1 , 
                                      outcome_prev = as.numeric(scenarios[i,3]), all_treat = FALSE, 
                                      all_control = FALSE, censor = F)
  
  initiators <- data %>% 
    filter(t == 0) %>% 
    summarise(CA = sum(A))
  events <- data %>% 
    summarise(CY = sum(Y))
  
  summary_table[i,] <- do.call(rbind, lapply(c(scenarios[i,1], scenarios[i,2],scenarios[i,3], initiators$CA/scenarios[i,1],events$CY/scenarios[i,1]), as.numeric))
}

colnames(summary_table) <- c('Sample size', 'Conf', 'alpha_y', 'Initiators', 'Events')
