#!/usr/bin R
library(dplyr)
library(tidyr)
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
source("continuous_timevarying_outcome_datagen.R")
library(TrialEmulation)
library(MASS)
library(sandwich)
library(foreach)
library(doParallel)
library(doRNG)
library(nleqslv)
source('calibration_func_trials.R')
set.seed(14052025)

simdata<-DATA_GEN_continous_timevarying_outcome_treatment_switch(ns = 200 ,nv = 3,treat_prev =  0.5,
                                                                 conf =  1.5,
                                                                 censor = F)
simdata <- simdata %>% 
  mutate(switch = ifelse(t == 0, 0,ifelse(A!=Ap,1,0)),
         missX1 = log(abs(X1))/4,
         missX2 = sqrt(abs(X2))/3) %>% 
  group_by(ID) %>% 
  dplyr::mutate(missCX1 = cumsum(missX1),
                missCX2 = cumsum(missX2),
                CA = cumsum(A))
#con4<-xtabs(~t + switch, data=simdata)
#ftable(con4)



######### Correctly specified ############### 
PP_prep <- TrialEmulation::data_preparation(simdata, id='ID', period='t', treatment='A', outcome='Y', cense = 'C',
                                            eligible ='eligible',
                                            estimand_type = 'PP',
                                            switch_d_cov = ~ X1 + X2,
                                            switch_n_cov = ~ -1,
                                            outcome_cov = ~ X1 + X2, model_var = c('assigned_treatment'),
                                            quiet = T,
                                            save_weight_models = F,
                                            data_dir = getwd())

switch_data <- PP_prep$data %>% 
  dplyr::mutate(t = trial_period + followup_time) %>%
  merge(simdata[,c('ID', 't', 'Y')], by.x = c('id', 't'), by.y = c('ID', 't')) %>% 
  dplyr::filter(trial_period == 0) %>% 
  dplyr::arrange(id, followup_time) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(Ap = ifelse(followup_time == 0, 0,lag(assigned_treatment)),
                CAp = cumsum(Ap),
                CA = cumsum(assigned_treatment))

#con4<-xtabs(~followup_time + assigned_treatment, data=switch_data)
#ftable(con4)

########## Outcome imputation correct and miss ############
library(data.table)
library(ltmle)
wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, value.var = c("A", "X1", "X2", "missX1", "missX2", "CA", "Y"))
regime_data <- array(1, dim = c(200,3,2))
regime_data[,,2] <- 0.0*regime_data[,,2]
summary_measures <- array(, dim = c(2,1,3))
summary_measures[1,1,] <- c(1,2,3)
summary_measures[2,1,] <- c(0,0,0)
colnames(summary_measures) <- 'cumA'

tmle <- ltmleMSM(data = wideSimdata[,.(X1_0,X2_0,A_0,Y_0,X1_1,X2_1, A_1, Y_1, X1_2, X2_2, A_2, Y_2)], Anodes = c('A_0', 'A_1', 'A_2'), 
                 Lnodes = c('X1_0','X2_0', 'X1_1', 'X2_1', 'X1_2', 'X2_2'), 
                 Ynodes = c('Y_0', 'Y_1', 'Y_2'),
                 Qform = c(X1_0 = 'Q.kplus1 ~ 1',
                           X2_0 = 'Q.kplus1 ~ 1',
                           Y_0 = 'Q.kplus1 ~ X1_0 + X2_0 + A_0',
                           X1_1 = 'Q.kplus1 ~ 1',
                           X2_1 = 'Q.kplus1 ~ 1',
                           Y_1 = 'Q.kplus1 ~ X1_0 + X2_0  + X1_1 + X2_1 + A_1',
                           X1_2 = 'Q.kplus1 ~ 1',
                           X2_2 = 'Q.kplus1 ~ 1',
                           Y_2 = 'Q.kplus1 ~ X1_0 + X2_0 +  X1_1 + X2_1 +  X1_2 + X2_2 + A_2'),
                 gform = c(A_0 = 'A_0 ~ 1', A_1 = 'A_1 ~ A_0 + X1_1 + X2_1', A_2 = 'A_2 ~ A_1 + X1_2 + X2_2'),
                 survivalOutcome = FALSE, 
                 regimes = regime_data,
                 working.msm = 'Y ~ cumA',
                 summary.measures = summary_measures,
                 final.Ynodes = c('Y_0', 'Y_1', 'Y_2')
                 )
summary(tmle)
print(tmle$beta*(max(simdata$Y) - min(simdata$Y)) + min(simdata$Y))
print(tmle$beta.iptw*(max(simdata$Y) - min(simdata$Y)) + min(simdata$Y))
