#!/usr/bin R
library(dplyr)
library(tidyr)
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("continuous_outcome_datagen.R")
library(TrialEmulation)
library(MASS)
library(sandwich)
library(foreach)
library(doParallel)
library(doRNG)
library(nleqslv)
source('calibration_func_trials.R')
set.seed(NULL)

iters <- 500
#l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
size <- c(200,500,1000)
treat <- c(-0.7,-0.5,-0.3)
conf <- c(0.5,0.9,1.5)

scenarios <- tidyr::crossing(size, treat,conf)
simdata<-DATA_GEN_continous_outcome_treatment_switch(ns = as.numeric(scenarios[l,1]),nv = 5,treat_prev =  0,
                                                           conf =  0.5,
                                                           censor = F)
simdata <- simdata %>% 
  mutate(switch = ifelse(t == 0, 0,ifelse(A!=Ap,1,0)))
con4<-xtabs(~t + switch, data=simdata)
ftable(con4)
PP_prep <- TrialEmulation::data_preparation(simdata, id='ID', period='t', treatment='A', outcome='Y', cense = 'C',
                                            eligible ='eligible',
                                            estimand_type = 'PP',
                                            switch_d_cov = ~ X1 + X2 + X3,
                                            switch_n_cov = ~ -1,
                                            outcome_cov = ~ X1 + X2 + X3, model_var = c('assigned_treatment'),
                                            quiet = T,
                                            save_weight_models = F,
                                            data_dir = getwd())

switch_data <- PP_prep$data %>% 
  dplyr::mutate(t = trial_period + followup_time) %>%
  merge(simdata[,c('ID', 't', 'Y')], by.x = c('id', 't'), by.y = c('ID', 't')) %>% 
  dplyr::filter(trial_period == 0) %>% 
  dplyr::arrange(id, followup_time) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(Ap = ifelse(followup_time == 0, 0,lag(assigned_treatment)))

con4<-xtabs(~followup_time + assigned_treatment, data=switch_data)
ftable(con4)
