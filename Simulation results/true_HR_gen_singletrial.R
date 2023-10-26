library(modelr)
library(tidyverse)
library(tidyr)
source("simulate_MSM_simplified.R")
set.seed(20222022)
library(MASS)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)

treat <- c(-1,0,1)
conf <- c(1,3,5)
outcome_prev <- c(-4.7,-3.8,-3)

scenarios <- tidyr::crossing(conf, treat)

true_HR <- array(,dim = c(9,3))

for (l in 1:9){
  for (j in 1:3){
    simdata_censored<-DATA_GEN_censored_reduced(1000000, 5, 
                                                      conf = as.numeric(scenarios[l,1]), 
                                                      treat_prev = as.numeric(scenarios[l,2]),
                                                      outcome_prev = outcome_prev[j],
                                                      censor = F)
    PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                                eligible ='eligible',
                                                switch_d_cov = ~X2 + X4,
                                                outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                                use_weight=T, use_censor=T, quiet = T,
                                                save_weight_models = F,
                                                data_dir = getwd())
    switch_data <- PP_prep$data %>% 
      dplyr::filter(trial_period == 0) %>% 
      dplyr::mutate(t_1 = ifelse(followup_time == 1,1,0),
                    t_2 = ifelse(followup_time == 2,1,0),
                    t_3 = ifelse(followup_time == 3,1,0),
                    t_4 = ifelse(followup_time == 4,1,0),
                    t_1A = t_1*assigned_treatment,
                    t_2A = t_2*assigned_treatment,
                    t_3A = t_3*assigned_treatment,
                    t_4A = t_4*assigned_treatment,
                    t_1X2 = t_1*X2,
                    t_2X2 = t_2*X2,
                    t_3X2 = t_3*X2,
                    t_4X2 = t_4*X2,
                    t_1X4 = t_1*X4,
                    t_2X4 = t_2*X4,
                    t_3X4 = t_3*X4,
                    t_4X4 = t_4*X4)
    PP <- TrialEmulation::trial_msm(data = switch_data,
                                    outcome_cov = ~ X2 + X4+ assigned_treatment+
                                      t_1 + t_2 + t_3 + t_4 +
                                      t_1A + t_2A + t_3A + t_4A +
                                      t_1X2 + t_2X2 + t_3X2 + t_4X2 +
                                      t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                    model_var = c('assigned_treatment'),
                                    glm_function = 'glm',
                                    include_trial_period = ~1, include_followup_time = ~1,
                                    use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
    true_HR[l,j] <- PP$model$coefficients[2]
  }
}
save(true_HR, file = "true_HR_singletrial.rda")

