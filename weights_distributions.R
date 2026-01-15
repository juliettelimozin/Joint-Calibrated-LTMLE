source("/home/juliette/Calibrated-weights-sequential-trial-emulation/dgm_biometrics_censoring.R")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/calibration_func_trials.R")
library(MASS)
library(foreach)
library(doParallel)
library(doRNG)
library(nleqslv)
library(ltmle)
library(data.table)
library(matrixStats)
library(xtable)

set.seed(160625)
seeds <- floor(runif(1000)*10^8)

weight_dist_generator_censored <- function(sample_size = 1000,  seeds = seeds,conf = 5,treat_prev_0 = 0,  treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1){
  set.seed(seeds[1])
  simdata<-DATA_GEN_censored(ns = sample_size, conf = conf,treat_prev_0 = treat_prev_0, 
                             treat_prev_d1_1 = treat_prev_d1_1, 
                             treat_prev_d0_1 = treat_prev_d0_1, 
                             treat_prev_d1_2 =treat_prev_d1_2, 
                             treat_prev_d0_2 = treat_prev_d0_2)
  
  
  treatment_model_pooled <- glm(A~Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t !=0,], family = 'binomial')
  
  treat_model_A_0 <- glm(A~X1 + X2 + X3 + X4-1, data = simdata[simdata$t == 0,], family = 'binomial')
  treat_model_A_1 <- glm(A~ Ap+X1 + X2 + X3 + X4, data = simdata[simdata$t == 1,], family = 'binomial')
  treat_model_A_2 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 2,], family = 'binomial')
  
  
  simdata$ps <- 1.0
  simdata[simdata$t == 0& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 1,], type = 'response'))
  simdata[simdata$t == 1& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 1,], type = 'response'))
  simdata[simdata$t == 2& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 1,], type = 'response'))
  
  simdata[simdata$t == 0& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 0,], type = 'response'))
  simdata[simdata$t == 1& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 0,], type = 'response'))
  simdata[simdata$t == 2& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 0,], type = 'response'))
  
  censor_model_0 <- glm(C ~ X1 + X2 + X3 + X4, data = simdata[simdata$t == 0,], family = 'binomial')
  censor_model_1 <- glm(C ~ X1 + X2 + X3 + X4, data = simdata[simdata$t == 1,], family = 'binomial')
  
  simdata$lagX1 <- ave(simdata$X1, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
  simdata$lagX2 <- ave(simdata$X2, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
  simdata$lagX3 <- ave(simdata$X3, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
  simdata$lagX4 <- ave(simdata$X4, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
  
  simdata$pr_c <- 1.0
  simdata[simdata$t == 1,]$pr_c <- as.vector(1-predict.glm(censor_model_0, 
                                                           newdata = data.frame(X1 = simdata[simdata$t == 1,]$lagX1,
                                                                                X2 = simdata[simdata$t == 1,]$lagX2,
                                                                                X3 = simdata[simdata$t == 1,]$lagX3,
                                                                                X4 = simdata[simdata$t == 1,]$lagX4),
                                                           type = 'response'))
  simdata[simdata$t == 2,]$pr_c <- as.vector(1-predict.glm(censor_model_1,
                                                           newdata = data.frame(X1 = simdata[simdata$t == 2,]$lagX1,
                                                                                X2 = simdata[simdata$t == 2,]$lagX2,
                                                                                X3 = simdata[simdata$t == 2,]$lagX3,
                                                                                X4 = simdata[simdata$t == 2,]$lagX4),
                                                           type = 'response'))
  
  simdata$treat_weight <- ave(simdata$ps, simdata$ID, FUN = function(X) 1/cumprod(X))
  simdata$censor_weight <- ave(simdata$pr_c, simdata$ID, FUN = function(X) 1/cumprod(X))
  simdata$weight <- simdata$treat_weight*simdata$censor_weight
  
  simdata$tall <- simdata$t
  simdata$sub <- simdata$ID
  simdata$RA <- 1
  simdata[simdata$t == 0 & !(simdata$CA == 1),]$RA <- 0
  simdata[simdata$t == 1 & !(simdata$CA == 2),]$RA <- 0
  simdata[simdata$t == 2 & !(simdata$CA == 3),]$RA <- 0
  simdata$tX1 <- simdata$t*simdata$X1
  simdata$tX2 <- simdata$t*simdata$X2
  simdata$tX3 <- simdata$t*simdata$X3
  simdata$tX4 <- simdata$t*simdata$X4
  
  
  calibrate_always_treated <- calibration_by_time_from_baseline(simdata, 
                                                                var = c("X1", "X2", "X3", "X4"), 
                                                                censor = TRUE, 
                                                                c_var = c("X1", "X2", "X3", "X4"),
                                                                weights_var = 'weight')
  calibrate_always_treated_aggr <- aggregated_calibration_from_baseline(simdata, 
                                                                        var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), 
                                                                        censor = TRUE, 
                                                                        c_var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"),
                                                                        weights_var = 'weight')
  
  calibrate_always_treated_treat_only <- calibration_by_time_from_baseline(simdata, 
                                                                           var = c("X1", "X2", "X3", "X4"), 
                                                                           censor = FALSE, 
                                                                           c_var = c("X1", "X2", "X3", "X4"),
                                                                           weights_var = 'weight')
  calibrate_always_treated_aggr_treat_only <- aggregated_calibration_from_baseline(simdata, 
                                                                                   var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), 
                                                                                   censor = FALSE, 
                                                                                   c_var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"),
                                                                                   weights_var = 'weight')
  
  calibrate_always_treated_censor_only <- calibration_by_time_from_baseline(simdata, 
                                                                            calibrate_treat = FALSE,
                                                                            var = c("X1", "X2", "X3", "X4"), 
                                                                            censor = TRUE, 
                                                                            c_var = c("X1", "X2", "X3", "X4"),
                                                                            weights_var = 'weight')
  calibrate_always_treated_aggr_censor_only <- aggregated_calibration_from_baseline(simdata, 
                                                                                    calibrate_treat = FALSE,
                                                                                    var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), 
                                                                                    censor = TRUE, 
                                                                                    c_var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"),
                                                                                    weights_var = 'weight')
  
  simdata$Cweights <- calibrate_always_treated$data$Cweights
  simdata$Cweights_aggr <- calibrate_always_treated_aggr$data$Cweights
  
  simdata$Cweights_treat_only <- calibrate_always_treated_treat_only$data$Cweights
  simdata$Cweights_aggr_treat_only <- calibrate_always_treated_aggr_treat_only$data$Cweights
  
  simdata$Cweights_censor_only <- calibrate_always_treated_censor_only$data$Cweights
  simdata$Cweights_aggr_censor_only <- calibrate_always_treated_aggr_censor_only$data$Cweights
  
  simdata$RA <- 1
  simdata[simdata$t == 0 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 1 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 2 & !(simdata$CA == 0),]$RA <- 0
  
  calibrate_never_treated <- calibration_by_time_from_baseline(simdata, 
                                                               var = c("X1", "X2", "X3", "X4"), 
                                                               censor = TRUE, 
                                                               c_var = c("X1", "X2", "X3", "X4"),
                                                               weights_var = 'Cweights')
  calibrate_never_treated_aggr <- aggregated_calibration_from_baseline(simdata, 
                                                                       var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), 
                                                                       censor = TRUE, 
                                                                       c_var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"),
                                                                       weights_var = 'Cweights_aggr')
  
  calibrate_never_treated_treat_only <- calibration_by_time_from_baseline(simdata, 
                                                                          var = c("X1", "X2", "X3", "X4"), 
                                                                          censor = FALSE, 
                                                                          c_var = c("X1", "X2", "X3", "X4"),
                                                                          weights_var = 'Cweights_treat_only')
  calibrate_never_treated_aggr_treat_only <- aggregated_calibration_from_baseline(simdata, 
                                                                                  var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), 
                                                                                  censor = FALSE, 
                                                                                  c_var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"),
                                                                                  weights_var = 'Cweights_aggr_treat_only')
  
  calibrate_never_treated_censor_only <- calibration_by_time_from_baseline(simdata, 
                                                                           calibrate_treat = FALSE,
                                                                           var = c("X1", "X2", "X3", "X4"), 
                                                                           censor = TRUE, 
                                                                           c_var = c("X1", "X2", "X3", "X4"),
                                                                           weights_var = 'Cweights_censor_only')
  calibrate_never_treated_aggr_censor_only <- aggregated_calibration_from_baseline(simdata, 
                                                                                   calibrate_treat = FALSE,
                                                                                   var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), 
                                                                                   censor = TRUE, 
                                                                                   c_var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"),
                                                                                   weights_var = 'Cweights_aggr_censor_only')
  
  simdata$Cweights <- calibrate_never_treated$data$Cweights
  simdata$Cweights_aggr <- calibrate_never_treated_aggr$data$Cweights
  
  simdata$Cweights_treat_only <- calibrate_never_treated_treat_only$data$Cweights
  simdata$Cweights_aggr_treat_only <- calibrate_never_treated_aggr_treat_only$data$Cweights
  
  simdata$Cweights_censor_only <- calibrate_never_treated_censor_only$data$Cweights
  simdata$Cweights_aggr_censor_only <- calibrate_never_treated_aggr_censor_only$data$Cweights
  
  simdata$weights <- simdata$weight
  
  simdata$RA <- 1
  simdata[simdata$t == 1 & !(simdata$CA == 2 | simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 2 & !(simdata$CA == 3 | simdata$CA == 0),]$RA <- 0
 
  print(xtable(as.data.frame(t(sapply(simdata[simdata$RA ==1,c('weights', 'Cweights')], summary)))), type = "latex")
}

weight_dist_generator <- function(sample_size = 1000,  seeds = seeds,conf = 5,treat_prev_0 = 0,  treat_prev_d1_1 = 0.1, treat_prev_d0_1 = -5, treat_prev_d1_2 = -5, treat_prev_d0_2 = -5.1){
  set.seed(seeds[1])
  simdata<-DATA_GEN(ns = sample_size, conf = conf,treat_prev_0 = treat_prev_0, 
                    treat_prev_d1_1 = treat_prev_d1_1, 
                    treat_prev_d0_1 = treat_prev_d0_1, 
                    treat_prev_d1_2 =treat_prev_d1_2, 
                    treat_prev_d0_2 = treat_prev_d0_2)
 
  treatment_model_pooled <- glm(A~Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t !=0,], family = 'binomial')
  
  treat_model_A_0 <- glm(A~X1 + X2 + X3 + X4-1, data = simdata[simdata$t == 0,], family = 'binomial')
  treat_model_A_1 <- glm(A~ Ap+X1 + X2 + X3 + X4, data = simdata[simdata$t == 1,], family = 'binomial')
  treat_model_A_2 <- glm(A~ Ap+ X1 + X2 + X3 + X4, data = simdata[simdata$t == 2,], family = 'binomial')
  
  
  simdata$ps <- 1.0
  simdata[simdata$t == 0& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 1,], type = 'response'))
  simdata[simdata$t == 1& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 1,], type = 'response'))
  simdata[simdata$t == 2& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 1,], type = 'response'))
  
  simdata[simdata$t == 0& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 0,], type = 'response'))
  simdata[simdata$t == 1& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 0,], type = 'response'))
  simdata[simdata$t == 2& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 0,], type = 'response'))
  
  simdata$weight <- ave(simdata$ps, simdata$ID, FUN = function(X) 1/cumprod(X))
  simdata$tall <- simdata$t
  simdata$sub <- simdata$ID
  simdata$RA <- 1
  simdata[simdata$t == 0 & !(simdata$CA == 1),]$RA <- 0
  simdata[simdata$t == 1 & !(simdata$CA == 2),]$RA <- 0
  simdata[simdata$t == 2 & !(simdata$CA == 3),]$RA <- 0
  simdata$tX1 <- simdata$t*simdata$X1
  simdata$tX2 <- simdata$t*simdata$X2
  simdata$tX3 <- simdata$t*simdata$X3
  simdata$tX4 <- simdata$t*simdata$X4
  
  calibrate_always_treated <- calibration_by_time_from_baseline(simdata, var = c("X1", "X2", "X3", "X4"), weights_var = 'weight')
  calibrate_always_treated_aggr <- aggregated_calibration_from_baseline(simdata, var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), weights_var = 'weight')
  
  simdata$Cweights <- calibrate_always_treated$data$Cweights
  simdata$Cweights_aggr <- calibrate_always_treated_aggr$data$Cweights
  
  simdata$RA <- 1
  simdata[simdata$t == 0 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 1 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 2 & !(simdata$CA == 0),]$RA <- 0
  
  calibrate_never_treated <- calibration_by_time_from_baseline(simdata, var = c("X1", "X2", "X3", "X4"), weights_var = 'Cweights')
  calibrate_never_treated_aggr <- aggregated_calibration_from_baseline(simdata, var = c("X1", "X2", "X3", "X4", "t", "tX1", "tX2", "tX3", "tX4"), weights_var = 'Cweights_aggr')
  
  simdata$Cweights <- calibrate_never_treated$data$Cweights
  simdata$Cweights_aggr <- calibrate_never_treated_aggr$data$Cweights
  
  simdata$weights <- simdata$weight
  
  simdata$RA <- 1
  simdata[simdata$t == 1 & !(simdata$CA == 2 | simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 2 & !(simdata$CA == 3 | simdata$CA == 0),]$RA <- 0
  
  print(xtable(as.data.frame(t(sapply(simdata[simdata$RA ==1,c('weights', 'Cweights')], summary)))), type = "latex")
}

cat('NO CENSORING \n')
cat('Weak confounding \n')
weight_dist_generator(sample_size = 300,seeds = seeds,conf = 0.2, 
                               treat_prev_0 = 0, 
                               treat_prev_d1_1 = 1, 
                               treat_prev_d0_1 = -1.25, 
                               treat_prev_d1_2 =0.8, 
                               treat_prev_d0_2 = -1.25)
weight_dist_generator(sample_size = 500,seeds = seeds,conf = 0.2, 
                               treat_prev_0 = 0, 
                               treat_prev_d1_1 = 1, 
                               treat_prev_d0_1 = -1.25, 
                               treat_prev_d1_2 =0.8, 
                               treat_prev_d0_2 = -1.25)
weight_dist_generator(sample_size = 1000,seeds = seeds,conf = 0.2, 
                               treat_prev_0 = 0, 
                               treat_prev_d1_1 = 1, 
                               treat_prev_d0_1 = -1.25, 
                               treat_prev_d1_2 =0.8, 
                               treat_prev_d0_2 = -1.25)
weight_dist_generator(sample_size = 2500,seeds = seeds,conf = 0.2, 
                               treat_prev_0 = 0, 
                               treat_prev_d1_1 = 1, 
                               treat_prev_d0_1 = -1.25, 
                               treat_prev_d1_2 =0.8, 
                               treat_prev_d0_2 = -1.25)

cat('-------------------------------------')
cat('Strong confounding \n')

weight_dist_generator(sample_size = 300,  
                               seeds = seeds,conf = 5, 
                               treat_prev_d1_1 = 0.1, 
                               treat_prev_d0_1 = -5, 
                               treat_prev_d1_2 = -5, 
                               treat_prev_d0_2 = -5.1)
weight_dist_generator(sample_size = 500,  
                               seeds = seeds,conf = 5, 
                               treat_prev_d1_1 = 0.1, 
                               treat_prev_d0_1 = -5, 
                               treat_prev_d1_2 = -5, 
                               treat_prev_d0_2 = -5.1)
weight_dist_generator(sample_size = 1000,  
                               seeds = seeds,conf = 5, 
                               treat_prev_d1_1 = 0.1, 
                               treat_prev_d0_1 = -5, 
                               treat_prev_d1_2 = -5, 
                               treat_prev_d0_2 = -5.1)
weight_dist_generator(sample_size = 2500,  
                               seeds = seeds,conf = 5, 
                               treat_prev_d1_1 = 0.1, 
                               treat_prev_d0_1 = -5, 
                               treat_prev_d1_2 = -5, 
                               treat_prev_d0_2 = -5.1)


cat('-------------------------------------------------------------------------')


cat('CENSORING \n')
cat('Weak confounding \n')
weight_dist_generator_censored(sample_size = 300,seeds = seeds,conf = 0.2, 
                               treat_prev_0 = 0, 
                               treat_prev_d1_1 = 1, 
                               treat_prev_d0_1 = -1.25, 
                               treat_prev_d1_2 =0.8, 
                               treat_prev_d0_2 = -1.25)
weight_dist_generator_censored(sample_size = 500,seeds = seeds,conf = 0.2, 
                               treat_prev_0 = 0, 
                               treat_prev_d1_1 = 1, 
                               treat_prev_d0_1 = -1.25, 
                               treat_prev_d1_2 =0.8, 
                               treat_prev_d0_2 = -1.25)
weight_dist_generator_censored(sample_size = 1000,seeds = seeds,conf = 0.2, 
                      treat_prev_0 = 0, 
                      treat_prev_d1_1 = 1, 
                      treat_prev_d0_1 = -1.25, 
                      treat_prev_d1_2 =0.8, 
                      treat_prev_d0_2 = -1.25)
weight_dist_generator_censored(sample_size = 2500,seeds = seeds,conf = 0.2, 
                               treat_prev_0 = 0, 
                               treat_prev_d1_1 = 1, 
                               treat_prev_d0_1 = -1.25, 
                               treat_prev_d1_2 =0.8, 
                               treat_prev_d0_2 = -1.25)

cat('-------------------------------------')
cat('Strong confounding \n')

weight_dist_generator_censored(sample_size = 300,  
                      seeds = seeds,conf = 5, 
                      treat_prev_d1_1 = 0.1, 
                      treat_prev_d0_1 = -5, 
                      treat_prev_d1_2 = -5, 
                      treat_prev_d0_2 = -5.1)
weight_dist_generator_censored(sample_size = 500,  
                               seeds = seeds,conf = 5, 
                               treat_prev_d1_1 = 0.1, 
                               treat_prev_d0_1 = -5, 
                               treat_prev_d1_2 = -5, 
                               treat_prev_d0_2 = -5.1)
weight_dist_generator_censored(sample_size = 1000,  
                               seeds = seeds,conf = 5, 
                               treat_prev_d1_1 = 0.1, 
                               treat_prev_d0_1 = -5, 
                               treat_prev_d1_2 = -5, 
                               treat_prev_d0_2 = -5.1)
weight_dist_generator_censored(sample_size = 2500,  
                               seeds = seeds,conf = 5, 
                               treat_prev_d1_1 = 0.1, 
                               treat_prev_d0_1 = -5, 
                               treat_prev_d1_2 = -5, 
                               treat_prev_d0_2 = -5.1)


