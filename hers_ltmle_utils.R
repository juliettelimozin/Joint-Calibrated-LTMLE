h_function <- function(data, treat_models_n,censor_models_n){
  data$hs <- 1.0
  for(m in 0:4){
    if (m!= 0){
      data[data$t == m& data$A == 1,]$hs = as.vector(predict.glm(treat_models_n[[m+1]],
                                                                 newdata = data[data$t == m& data$A == 1,], 
                                                                 type = 'response')*(1-predict.glm(censor_models_n[[m]],
                                                                                                   newdata = data[data$t == m& data$A == 1,],
                                                                                                   type = 'response')))
      data[data$t == m& data$A == 0,]$hs = as.vector((1-predict.glm(treat_models_n[[m+1]],
                                                                    newdata = data[data$t == m& data$A == 0,], 
                                                                    type = 'response'))*(1-predict.glm(censor_models_n[[m]],
                                                                                                       newdata = data[data$t == m& data$A == 0,],
                                                                                                       type = 'response')))
      
    } else {
      data[data$t == m& data$A == 1,]$hs = as.vector(predict.glm(treat_models_n[[m+1]],
                                                                 newdata = data[data$t == m& data$A == 1,], 
                                                                 type = 'response'))
      data[data$t == m& data$A == 0,]$hs = as.vector(1-predict.glm(treat_models_n[[m+1]],
                                                                   newdata = data[data$t == m& data$A == 0,], 
                                                                   type = 'response'))
    }
  }
  data$h_weight <- 1.0
  data[data$A == 1,]$h_weight<- ave(data[data$A == 1,]$hs, data[data$A == 1,]$ID, FUN = function(X) cumprod(X))
  data[data$A == 0,]$h_weight<- ave(data[data$A == 0,]$hs, data[data$A == 0,]$ID, FUN = function(X) cumprod(X))
  
  return(data)
}

initial_weight_calculation <- function(simdata){
  treat_model_A_0 <- glm(A~ sqrtCD4_1 + sqrtCD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 0,], family = 'binomial')
  treat_model_A_1 <- glm(A~ Ap+sqrtCD4_1 + sqrtCD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 1,], family = 'binomial')
  treat_model_A_2 <- glm(A~ Ap+sqrtCD4_1 + sqrtCD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 2,], family = 'binomial')
  treat_model_A_3 <- glm(A~ Ap+sqrtCD4_1 + sqrtCD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 3,], family = 'binomial')
  treat_model_A_4 <- glm(A~ Ap+sqrtCD4_1 + sqrtCD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 4,], family = 'binomial')
  
  treat_model_A_0_num <- glm(A~ D, data = simdata[simdata$t == 0,], family = 'binomial')
  treat_model_A_1_num <- glm(A~ Ap + D, data = simdata[simdata$t == 1,], family = 'binomial')
  treat_model_A_2_num <- glm(A~ Ap + D, data = simdata[simdata$t == 2,], family = 'binomial')
  treat_model_A_3_num <- glm(A~ Ap + D, data = simdata[simdata$t == 3,], family = 'binomial')
  treat_model_A_4_num <- glm(A~ Ap + D, data = simdata[simdata$t == 4,], family = 'binomial')
  
  treat_models_d <- list(treat_model_A_0,treat_model_A_1,treat_model_A_2,treat_model_A_3,treat_model_A_4)
  treat_models_n <- list(treat_model_A_0_num,treat_model_A_1_num,treat_model_A_2_num,treat_model_A_3_num,treat_model_A_4_num)
  
  censor_model_0 <- glm(C ~ sqrtCD4_1 + sqrtCD4  + viral + viral_1  + HIVsym + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 0,], family = 'binomial')
  censor_model_1 <- glm(C ~ sqrtCD4_1 + sqrtCD4  + viral + viral_1  + HIVsym + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 1,], family = 'binomial')
  censor_model_2 <- glm(C ~ sqrtCD4_1 + sqrtCD4  + viral + viral_1  + HIVsym + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 2,], family = 'binomial')
  censor_model_3 <- glm(C ~ sqrtCD4_1 + sqrtCD4  + viral + viral_1  + HIVsym + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 3,], family = 'binomial')
  
  censor_model_0_num <- glm(C ~ D, data = simdata[simdata$t == 0,], family = 'binomial')
  censor_model_1_num <- glm(C ~ D, data = simdata[simdata$t == 1,], family = 'binomial')
  censor_model_2_num <- glm(C ~ D, data = simdata[simdata$t == 2,], family = 'binomial')
  censor_model_3_num <- glm(C ~ D, data = simdata[simdata$t == 3,], family = 'binomial')
  
  censor_models_d <- list(censor_model_0, censor_model_1, censor_model_2, censor_model_3)
  censor_models_n <- list(censor_model_0_num, censor_model_1_num, censor_model_2_num, censor_model_3_num)
  
  simdata$lagsqrtCD4_1 <- simdata$sqrtCD4_2
  simdata$lagsqrtCD4 <- simdata$sqrtCD4_1
  simdata$lagviral <- simdata$viral_1
  simdata$lagviral_1 <- simdata$viral_2
  simdata$lagHIVsym <- simdata$HIVsym_1
  simdata$lagHIVsym_1 <- simdata$HIVsym_2
  simdata$lagSITE2 <- simdata$SITE2
  simdata$lagSITE3 <- simdata$SITE3
  simdata$lagWHITE <- simdata$WHITE
  simdata$lagOTHER <- simdata$OTHER
  
  simdata$ps <- 1.0
  simdata$pr_c <- 1.0
  
  for(m in 0:4){
    
    simdata[simdata$t == m& simdata$A == 1,]$ps = as.vector(predict.glm(treat_models_d[[m+1]],
                                                                        newdata = simdata[simdata$t == m& simdata$A == 1,], 
                                                                        type = 'response')/predict.glm(treat_models_n[[m+1]], 
                                                                                                       newdata = simdata[simdata$t == m& simdata$A == 1,], 
                                                                                                       type = 'response'))
    simdata[simdata$t == m& simdata$A == 0,]$ps = as.vector((1-predict.glm(treat_models_d[[m+1]], 
                                                                           newdata = simdata[simdata$t == m& simdata$A == 0,], 
                                                                           type = 'response'))/(1-predict.glm(treat_models_n[[m+1]], 
                                                                                                              newdata = simdata[simdata$t == m& simdata$A == 0,], 
                                                                                                              type = 'response')))
    if (m!= 0){
      simdata[simdata$t == m,]$pr_c <- as.vector((1-predict.glm(censor_models_d[[m]], 
                                                                newdata = data.frame(sqrtCD4_1 = simdata[simdata$t == m,]$lagsqrtCD4_1,
                                                                                     sqrtCD4 = simdata[simdata$t == m,]$lagsqrtCD4,
                                                                                     viral = simdata[simdata$t == m,]$lagviral,
                                                                                     viral_1 = simdata[simdata$t == m,]$lagviral_1,
                                                                                     HIVsym = simdata[simdata$t == m,]$lagHIVsym,
                                                                                     HIVsym_1 = simdata[simdata$t == m,]$lagHIVsym_1,
                                                                                     SITE2 = simdata[simdata$t == m,]$SITE2,
                                                                                     SITE3 = simdata[simdata$t == m,]$SITE3,
                                                                                     WHITE = simdata[simdata$t == m,]$WHITE,
                                                                                     OTHER = simdata[simdata$t == m,]$OTHER),
                                                                type = 'response'))/(1-predict.glm(censor_models_n[[m]], 
                                                                                                   newdata = simdata[simdata$t == m,],
                                                                                                   type = 'response')))
    }
  }
  
  simdata$treat_weight <- ave(simdata$ps, simdata$ID, FUN = function(X) 1/cumprod(X))
  simdata$censor_weight <- ave(simdata$pr_c, simdata$ID, FUN = function(X) 1/cumprod(X))
  simdata$weight <- simdata$treat_weight*simdata$censor_weight
  
  return(list(data = simdata, treat_models_n = treat_models_n, censor_models_n = censor_models_n))
}

weight_calibration <- function(simdata){

  simdata$lagsqrtCD4_1 <- simdata$sqrtCD4_2
  simdata$lagsqrtCD4 <- simdata$sqrtCD4_1
  simdata$lagviral <- simdata$viral_1
  simdata$lagviral_1 <- simdata$viral_2
  simdata$lagHIVsym <- simdata$HIVsym_1
  simdata$lagHIVsym_1 <- simdata$HIVsym_2
  simdata$lagSITE2 <- simdata$SITE2
  simdata$lagSITE3 <- simdata$SITE3
  simdata$lagWHITE <- simdata$WHITE
  simdata$lagOTHER <- simdata$OTHER
  
  simdata$tall <- simdata$t
  simdata$sub <- simdata$ID
  simdata$RA <- 1
  simdata[simdata$t == 0 & !(simdata$CA == 1),]$RA <- 0
  simdata[simdata$t == 1 & !(simdata$CA == 2),]$RA <- 0
  simdata[simdata$t == 2 & !(simdata$CA == 3),]$RA <- 0
  simdata[simdata$t == 3 & !(simdata$CA == 4),]$RA <- 0
  simdata[simdata$t == 4 & !(simdata$CA == 5),]$RA <- 0
  
  simdata$tsqrtCD4_1     <- simdata$t * simdata$sqrtCD4_1
  simdata$tsqrtCD4_2     <- simdata$t * simdata$sqrtCD4_2
  simdata$tsqrtCD4       <- simdata$t * simdata$sqrtCD4
  simdata$tviral     <- simdata$t * simdata$viral
  simdata$tviral_1   <- simdata$t * simdata$viral_1
  simdata$tviral_2   <- simdata$t * simdata$viral_2
  simdata$tHIVsym    <- simdata$t * simdata$HIVsym
  simdata$tHIVsym_1  <- simdata$t * simdata$HIVsym_1
  simdata$tSITE2     <- simdata$t * simdata$SITE2
  simdata$tSITE3     <- simdata$t * simdata$SITE3
  simdata$tWHITE     <- simdata$t * simdata$WHITE
  simdata$tOTHER     <- simdata$t * simdata$OTHER
  
  
  calibrate_always_treated <- calibration_by_time_from_baseline(simdata, 
                                                                var = c("sqrtCD4_1","sqrtCD4_2","viral_1","viral_2","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"), 
                                                                censor = TRUE, 
                                                                c_var = c("sqrtCD4_1","sqrtCD4","viral","viral_1","HIVsym","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"
                                                                ),
                                                                weights_var = 'weight')
  calibrate_always_treated_aggr <- aggregated_calibration_from_baseline(simdata, 
                                                                        var = c(
                                                                          "sqrtCD4_1", "sqrtCD4_2", "viral_1", "viral_2", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER", 
                                                                          "tsqrtCD4_1", "tsqrtCD4_2", "tviral_1", "tviral_2", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                        ), 
                                                                        censor = TRUE, 
                                                                        c_var = c(
                                                                          "sqrtCD4_1", "sqrtCD4", "viral", "viral_1", "HIVsym", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER",
                                                                          "tsqrtCD4_1", "tsqrtCD4", "tviral", "tviral_1", "tHIVsym", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                        ),
                                                                        weights_var = 'weight')
  
  calibrate_always_treated_treat_only <- calibration_by_time_from_baseline(simdata, 
                                                                           var = c("sqrtCD4_1","sqrtCD4_2","viral_1","viral_2","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"), 
                                                                           censor = FALSE, 
                                                                           c_var = c("sqrtCD4_1","sqrtCD4","viral","viral_1","HIVsym","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"
                                                                           ),
                                                                           weights_var = 'weight')
  calibrate_always_treated_aggr_treat_only <- aggregated_calibration_from_baseline(simdata, 
                                                                                   var = c(
                                                                                     "sqrtCD4_1", "sqrtCD4_2", "viral_1", "viral_2", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER", 
                                                                                     "tsqrtCD4_1", "tsqrtCD4_2", "tviral_1", "tviral_2", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                                   ), 
                                                                                   censor = FALSE, 
                                                                                   c_var = c(
                                                                                     "sqrtCD4_1", "sqrtCD4", "viral", "viral_1", "HIVsym", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER",
                                                                                     "tsqrtCD4_1", "tsqrtCD4", "tviral", "tviral_1", "tHIVsym", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                                   ),
                                                                                   weights_var = 'weight')
  
  
  simdata$Cweights <- calibrate_always_treated$data$Cweights
  simdata$Cweights_aggr <- calibrate_always_treated_aggr$data$Cweights
  
  simdata$Cweights_treat_only <- calibrate_always_treated_treat_only$data$Cweights
  simdata$Cweights_aggr_treat_only <- calibrate_always_treated_aggr_treat_only$data$Cweights
  
  simdata$RA <- 1
  simdata[simdata$t == 0 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 1 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 2 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 3 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 4 & !(simdata$CA == 0),]$RA <- 0
  
  calibrate_never_treated <- calibration_by_time_from_baseline(simdata, 
                                                               var = c("sqrtCD4_1","sqrtCD4_2","viral_1","viral_2","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"), 
                                                               censor = TRUE, 
                                                               c_var = c("sqrtCD4_1","sqrtCD4","viral","viral_1","HIVsym","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"
                                                               ),
                                                               weights_var = 'Cweights')
  calibrate_never_treated_aggr <- aggregated_calibration_from_baseline(simdata, 
                                                                       var = c(
                                                                         "sqrtCD4_1", "sqrtCD4_2", "viral_1", "viral_2", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER", 
                                                                         "tsqrtCD4_1", "tsqrtCD4_2", "tviral_1", "tviral_2", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                       ), 
                                                                       censor = TRUE, 
                                                                       c_var = c(
                                                                         "sqrtCD4_1", "sqrtCD4", "viral", "viral_1", "HIVsym", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER", 
                                                                         "tsqrtCD4_1", "tsqrtCD4", "tviral", "tviral_1", "tHIVsym", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                       ),
                                                                       weights_var = 'Cweights_aggr')
  
  calibrate_never_treated_treat_only <- calibration_by_time_from_baseline(simdata, 
                                                                          var = c("sqrtCD4_1","sqrtCD4_2","viral_1","viral_2","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"),
                                                                          censor = FALSE, 
                                                                          c_var = c("sqrtCD4_1","sqrtCD4","viral","viral_1","HIVsym","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"
                                                                          ),
                                                                          weights_var = 'Cweights_treat_only')
  calibrate_never_treated_aggr_treat_only <- aggregated_calibration_from_baseline(simdata, 
                                                                                  var = c(
                                                                                    "sqrtCD4_1", "sqrtCD4_2", "viral_1", "viral_2", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER", 
                                                                                    "tsqrtCD4_1", "tsqrtCD4_2", "tviral_1", "tviral_2", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                                  ), 
                                                                                  censor = FALSE, 
                                                                                  c_var = c(
                                                                                    "sqrtCD4_1", "sqrtCD4", "viral", "viral_1", "HIVsym", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER", 
                                                                                    "tsqrtCD4_1", "tsqrtCD4", "tviral", "tviral_1", "tHIVsym", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                                  ),
                                                                                  weights_var = 'Cweights_aggr_treat_only')
  
  
  
  simdata$Cweights <- calibrate_never_treated$data$Cweights
  simdata$Cweights_aggr <- calibrate_never_treated_aggr$data$Cweights
  
  simdata$Cweights_treat_only <- calibrate_never_treated_treat_only$data$Cweights
  simdata$Cweights_aggr_treat_only <- calibrate_never_treated_aggr_treat_only$data$Cweights
  
  
  simdata$weights <- simdata$weight
  return(simdata)
}

fit_msm <- function(wideSimdata, 
                    Qstars, 
                    Qstars_cali, 
                    Qstars_cali_aggr,
                    Qstars_cali_treat_only,
                    Qstars_cali_aggr_treat_only,
                    treat_models_n,
                    censor_models_n){

  t <- c(rep(0,dim(wideSimdata)[1]*2), 
         rep(1,dim(wideSimdata)[1]*2), 
         rep(2,dim(wideSimdata)[1]*2),
         rep(3,dim(wideSimdata)[1]*2),
         rep(4,dim(wideSimdata)[1]*2))
  
  msm_fitting_data_transformed <- data.frame(ID = rep(wideSimdata$ID,10), 
                                             t  = t,
                                             Y =  unlist(Qstars[, 1]), 
                                             Y_cali =  unlist(Qstars_cali[, 1]),
                                             Y_cali_aggr =  unlist(Qstars_cali_aggr[, 1]),
                                             Y_cali_treat_only =  unlist(Qstars_cali_treat_only[, 1]),
                                             Y_cali_aggr_treat_only =  unlist(Qstars_cali_aggr_treat_only[, 1]),
                                             A = c(rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                   rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                   rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                   rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                   rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1])),
                                             Ap = c(rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                    rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                    rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                    rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                    rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1])),
                                             D0 = as.numeric(rep(wideSimdata$D_0, 10) == 0),
                                             D1 = as.numeric(rep(wideSimdata$D_0, 10) == 1),
                                             D2 = as.numeric(rep(wideSimdata$D_0, 10) == 2),
                                             D = rep(wideSimdata$D_0, 10),
                                             cumA = c(rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                      rep(2,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                      rep(3,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                      rep(4,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                      rep(5,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1])),
                                             intercept_0 = as.numeric(t == 0),
                                             intercept_1 = as.numeric(t == 1),
                                             intercept_2 = as.numeric(t == 2),
                                             intercept_3 = as.numeric(t == 3),
                                             intercept_4 = as.numeric(t == 4))
  
  msm_fitting_data_transformed <- h_function(msm_fitting_data_transformed, treat_models_n = treat_models_n, censor_models_n = censor_models_n)
  
  msm_transformed <- glm(Y ~ intercept_0 + intercept_1 + 
                           intercept_2 + intercept_3 + intercept_4 + 
                           D1 + D2 + D0:cumA + D1:cumA + D2:cumA  - 1, data = msm_fitting_data_transformed)
  msm_transformed_cali <- glm(Y_cali ~ intercept_0 + intercept_1 + 
                                intercept_2 + intercept_3 + intercept_4 + 
                                D1 + D2 + D0:cumA + D1:cumA + D2:cumA   - 1, data = msm_fitting_data_transformed)
  msm_transformed_cali_aggr <- glm(Y_cali_aggr ~ intercept_0 + intercept_1 + 
                                     intercept_2 + intercept_3 + intercept_4 + 
                                     D1 + D2 + D0:cumA + D1:cumA + D2:cumA  - 1, data = msm_fitting_data_transformed)
  msm_transformed_cali_treat_only <- glm(Y_cali_treat_only ~ intercept_0 + intercept_1 + 
                                           intercept_2 + intercept_3 + intercept_4 + 
                                           D1 + D2 + D0:cumA + D1:cumA + D2:cumA   - 1, data = msm_fitting_data_transformed)
  msm_transformed_cali_aggr_treat_only <- glm(Y_cali_aggr_treat_only ~ intercept_0 + intercept_1 + 
                                                intercept_2 + intercept_3 + intercept_4 + 
                                                D1 + D2 + D0:cumA + D1:cumA + D2:cumA  - 1, data = msm_fitting_data_transformed)
  
  
  MSM_estimates <- rbind(msm_transformed$coefficients,
                         msm_transformed_cali_treat_only$coefficients,
                         msm_transformed_cali_aggr_treat_only$coefficients,
                         msm_transformed_cali$coefficients,
                         msm_transformed_cali_aggr$coefficients
  )
  
  rownames(MSM_estimates) <- c('MLE LTMLE',
                               'CMLE LTMLE treat. only',
                               'Aggr. CMLE LTMLE treat. only',
                               'CMLE LTMLE',
                               'Aggr. CMLE LTMLE')
  colnames(MSM_estimates) <- names(msm_transformed$coefficients)
  
  fitting_data_strata0 <- msm_fitting_data_transformed[msm_fitting_data_transformed$D1 == 0 & msm_fitting_data_transformed$D2 == 0,]
  fitting_data_strata1 <- msm_fitting_data_transformed[msm_fitting_data_transformed$D1 == 1,]
  fitting_data_strata2 <- msm_fitting_data_transformed[msm_fitting_data_transformed$D2 == 1,]
  
  EY_always_treated_strata0 <- EY_always_treated_strata1 <- EY_always_treated_strata2<- EY_never_treated_strata0 <- EY_never_treated_strata1 <- EY_never_treated_strata2<-array(,dim = c(5,5))
  
  rownames(EY_always_treated_strata0) <- rownames(EY_always_treated_strata1) <- rownames(EY_always_treated_strata2) <-  c('MLE LTMLE',
                                                                                                                          'CMLE LTMLE treat. only',
                                                                                                                          'Aggr. CMLE LTMLE treat. only',
                                                                                                                          'CMLE LTMLE',
                                                                                                                          'Aggr. CMLE LTMLE')
  rownames(EY_never_treated_strata0) <- rownames(EY_never_treated_strata1) <- rownames(EY_never_treated_strata2) <- c('MLE LTMLE',
                                                                                                                      'CMLE LTMLE treat. only',
                                                                                                                      'Aggr. CMLE LTMLE treat. only',
                                                                                                                      'CMLE LTMLE',
                                                                                                                      'Aggr. CMLE LTMLE')
  
  EY_always_treated_strata0[1,] <- tapply(predict.glm(msm_transformed, newdata = fitting_data_strata0[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata0[4,] <- tapply(predict.glm(msm_transformed_cali, newdata = fitting_data_strata0[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata0[5,] <- tapply(predict.glm(msm_transformed_cali_aggr, newdata = fitting_data_strata0[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata0[2,] <- tapply(predict.glm(msm_transformed_cali_treat_only, newdata = fitting_data_strata0[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata0[3,] <- tapply(predict.glm(msm_transformed_cali_aggr_treat_only, newdata = fitting_data_strata0[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA !=0,]$t, mean)
  
  EY_always_treated_strata1[1,] <- tapply(predict.glm(msm_transformed, newdata = fitting_data_strata1[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata1[4,] <- tapply(predict.glm(msm_transformed_cali, newdata = fitting_data_strata1[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata1[5,] <- tapply(predict.glm(msm_transformed_cali_aggr, newdata = fitting_data_strata1[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata1[2,] <- tapply(predict.glm(msm_transformed_cali_treat_only, newdata = fitting_data_strata1[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata1[3,] <- tapply(predict.glm(msm_transformed_cali_aggr_treat_only, newdata = fitting_data_strata1[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA !=0,]$t, mean)
  
  EY_always_treated_strata2[1,] <- tapply(predict.glm(msm_transformed, newdata = fitting_data_strata2[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata2[4,] <- tapply(predict.glm(msm_transformed_cali, newdata = fitting_data_strata2[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata2[5,] <- tapply(predict.glm(msm_transformed_cali_aggr, newdata = fitting_data_strata2[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata2[2,] <- tapply(predict.glm(msm_transformed_cali_treat_only, newdata = fitting_data_strata2[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA !=0,]$t, mean)
  EY_always_treated_strata2[3,] <- tapply(predict.glm(msm_transformed_cali_aggr_treat_only, newdata = fitting_data_strata2[fitting_data_strata0$cumA !=0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA !=0,]$t, mean)
  
  EY_never_treated_strata0[1,] <- tapply(predict.glm(msm_transformed, newdata = fitting_data_strata0[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata0[4,] <- tapply(predict.glm(msm_transformed_cali, newdata = fitting_data_strata0[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata0[5,] <- tapply(predict.glm(msm_transformed_cali_aggr, newdata = fitting_data_strata0[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata0[2,] <- tapply(predict.glm(msm_transformed_cali_treat_only, newdata = fitting_data_strata0[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata0[3,] <- tapply(predict.glm(msm_transformed_cali_aggr_treat_only, newdata = fitting_data_strata0[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata0[fitting_data_strata0$cumA ==0,]$t, mean)
  
  EY_never_treated_strata1[1,] <- tapply(predict.glm(msm_transformed, newdata = fitting_data_strata1[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata1[4,] <- tapply(predict.glm(msm_transformed_cali, newdata = fitting_data_strata1[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata1[5,] <- tapply(predict.glm(msm_transformed_cali_aggr, newdata = fitting_data_strata1[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata1[2,] <- tapply(predict.glm(msm_transformed_cali_treat_only, newdata = fitting_data_strata1[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata1[3,] <- tapply(predict.glm(msm_transformed_cali_aggr_treat_only, newdata = fitting_data_strata1[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata1[fitting_data_strata0$cumA ==0,]$t, mean)
  
  EY_never_treated_strata2[1,] <- tapply(predict.glm(msm_transformed, newdata = fitting_data_strata2[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata2[4,] <- tapply(predict.glm(msm_transformed_cali, newdata = fitting_data_strata2[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata2[5,] <- tapply(predict.glm(msm_transformed_cali_aggr, newdata = fitting_data_strata2[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata2[2,] <- tapply(predict.glm(msm_transformed_cali_treat_only, newdata = fitting_data_strata2[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA ==0,]$t, mean)
  EY_never_treated_strata2[3,] <- tapply(predict.glm(msm_transformed_cali_aggr_treat_only, newdata = fitting_data_strata2[fitting_data_strata0$cumA ==0,], type = 'response'), fitting_data_strata2[fitting_data_strata0$cumA ==0,]$t, mean)
  
  counterfactual_means <- list(EY_always_treated_strata0 = EY_always_treated_strata0,
                               EY_always_treated_strata1 = EY_always_treated_strata1,
                               EY_always_treated_strata2 = EY_always_treated_strata2,
                               EY_never_treated_strata0 = EY_never_treated_strata0,
                               EY_never_treated_strata1 = EY_never_treated_strata1,
                               EY_never_treated_strata2 = EY_never_treated_strata2)  
  
  return(list(MSM_estimates = MSM_estimates, counterfactual_means = counterfactual_means))
}

MyLtmleMSM_modified_hers <- function(simdata, initial_Q_t = NA, refit_weights = TRUE, treat_models_n = NA, censor_models_n = NA){
  if (refit_weights){
    calculate_weights <- initial_weight_calculation(simdata)
   simdata <- calculate_weights$data
   treat_models_n <- calculate_weights$treat_models_n
   censor_models_n <- calculate_weights$censor_models_n
   simdata <- weight_calibration(simdata)
  }
  
  wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, 
                                   value.var = c("A", "Ap", "App", "CD4_1","sqrtCD4_1", "sqrtCD4_2", "viral", "viral_1", "viral_2","HIVsym","HIVsym_1", "HIVsym_2", 
                                                 "SITE2", "SITE3", "WHITE", "OTHER", "CA", "CD4","D", "weights", "Cweights", "Cweights_aggr",
                                                 "Cweights_treat_only", "Cweights_aggr_treat_only"))
  
  a = min(simdata$CD4,na.rm = TRUE); b = max(simdata$CD4, na.rm = TRUE)
  
  wideSimdata$Y_4_scaled = (wideSimdata$CD4_4-a)/(b-a)
  wideSimdata$Y_3_scaled = (wideSimdata$CD4_3-a)/(b-a)
  wideSimdata$Y_2_scaled = (wideSimdata$CD4_2-a)/(b-a)
  wideSimdata$Y_1_scaled = (wideSimdata$CD4_1-a)/(b-a)
  wideSimdata$Y_0_scaled = (wideSimdata$CD4_0-a)/(b-a)
  
  if (all(is.na(initial_Q_t))){
    Q_fits <- list()
  }
  
  logitQforms <- vector("list", 5 * 5)
  dim(logitQforms) <- c(5, 5)
  
  Qstars <- vector("list", 5 * 5)
  dim(Qstars) <- c(5, 5)
  
  Qstars_cali <- vector("list", 5 * 5)
  dim(Qstars_cali) <- c(5, 5)
  
  Qstars_cali_aggr <- vector("list", 5 * 5)
  dim(Qstars_cali_aggr) <- c(5, 5)
  
  Qstars_cali_treat_only <- vector("list", 5 * 5)
  dim(Qstars_cali_treat_only) <- c(5, 5)
  
  Qstars_cali_aggr_treat_only <- vector("list", 5 * 5)
  dim(Qstars_cali_aggr_treat_only) <- c(5, 5)
  ########### Initial Q fits ##########
  # First loop: get the Q forms
  
  
  for (k in 0:4){
    for (j in k:0){
      if (k == j){
        formula_string <- paste0('CD4_',j,' ~ A_',j, '+ Ap_',j,' + App_',j,
                                 '+ sqrtCD4_1_',j,' + sqrtCD4_2_',j,
                                 ' + viral_',j,' + viral_1_',j,' + viral_2_',j,
                                 ' + HIVsym_',j,' + HIVsym_1_',j,' + HIVsym_2_',j,
                                 ' + SITE2_',j,' + SITE3_',j,' + WHITE_',j,' + OTHER_',j)
        
        predict_data <- data.frame(c(rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1])),c(rep(1, dim(wideSimdata)[1]), rep(0, dim(wideSimdata)[1])),c(rep(1, dim(wideSimdata)[1]), rep(0, dim(wideSimdata)[1])))
        colnames(predict_data) <- paste0(c('A_', 'Ap_', 'App_'),j)
        
        vars <- attr(terms(as.formula(formula_string)), "term.labels")
        
        # strip "_4" suffix
        vars_clean <- gsub(paste0("_",j,"$"), "", vars)
        for (v in vars_clean[-1:-3]) {
          predict_data[, paste0(v, "_", j)] <- rep(wideSimdata[[paste0(v, "_", j)]], 2)
        }
        
        if (all(is.na(initial_Q_t))){
          
          Qfit <- glm(data = wideSimdata, formula = as.formula(formula_string))
          Q_fits[[paste0('Q',k,"_",j,'_fit')]] <- Qfit
        } else {
          Qfit <- initial_Q_t[[paste0('Q',k,"_",j,'_fit')]]
        }
        
        logitQforms[[(k+1),(j+1)]] <- qlogis((predict.glm(Qfit, newdata = predict_data)-a)/(b-a))
        
        
      } else{ 
        fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                                   Qk_j1 = plogis(logitQforms[[(k+1),(j+2)]])*(b-a)+a, 
                                   A_j       = rep(wideSimdata[[paste0('A_', j)]], 2),
                                   Ap_j      = rep(wideSimdata[[paste0('Ap_', j)]], 2),
                                   App_j     = rep(wideSimdata[[paste0('App_', j)]], 2),
                                   sqrtCD4_1_j   = rep(wideSimdata[[paste0('sqrtCD4_1_', j)]], 2),
                                   sqrtCD4_2_j   = rep(wideSimdata[[paste0('sqrtCD4_2_', j)]], 2),
                                   viral_j   = rep(wideSimdata[[paste0('viral_', j)]], 2),
                                   viral_1_j = rep(wideSimdata[[paste0('viral_1_', j)]], 2),
                                   viral_2_j = rep(wideSimdata[[paste0('viral_2_', j)]], 2),
                                   HIVsym_j   = rep(wideSimdata[[paste0('HIVsym_', j)]], 2),
                                   HIVsym_1_j = rep(wideSimdata[[paste0('HIVsym_1_', j)]], 2),
                                   HIVsym_2_j = rep(wideSimdata[[paste0('HIVsym_2_', j)]], 2),
                                   SITE2_j   = rep(wideSimdata[[paste0('SITE2_', j)]], 2),
                                   SITE3_j   = rep(wideSimdata[[paste0('SITE3_', j)]], 2),
                                   WHITE_j   = rep(wideSimdata[[paste0('WHITE_', j)]], 2),
                                   OTHER_j   = rep(wideSimdata[[paste0('OTHER_', j)]], 2))
        
        Q_predict_data_d1 <- fitting_data[1:dim(wideSimdata)[1],]
        Q_predict_data_d1$A_j <- rep(1,dim(wideSimdata)[1])
        Q_predict_data_d1$Ap_j <- rep(1,dim(wideSimdata)[1])
        Q_predict_data_d1$App_j <- rep(1,dim(wideSimdata)[1])
        
        Q_predict_data_d0 <- fitting_data[1:dim(wideSimdata)[1],]
        Q_predict_data_d0$A_j <- rep(0,dim(wideSimdata)[1])
        Q_predict_data_d0$Ap_j <- rep(0,dim(wideSimdata)[1])
        Q_predict_data_d0$App_j <- rep(0,dim(wideSimdata)[1])
        
        
        if (all(is.na(initial_Q_t))){
          formula_string <- paste0("Qk_j1 ~ A_j + Ap_j + App_j +",
                                   "sqrtCD4_1_j + sqrtCD4_2_j +",
                                   "viral_j + viral_1_j + viral_2_j +",
                                   "HIVsym_j + HIVsym_1_j + HIVsym_2_j +",
                                   "SITE2_j + SITE3_j + WHITE_j + OTHER_j") 
          
          Qform_d1 <- glm(data = fitting_data[1:dim(wideSimdata)[1],], formula = as.formula(formula_string))
          Qform_d0 <- glm(data = fitting_data[(dim(wideSimdata)[1]+1):(dim(wideSimdata)[1]*2),], formula = as.formula(formula_string))
          
          Q_fits[[paste0('Q',k,"_",j,'_fit_d1')]]<- Qform_d1
          Q_fits[[paste0('Q',k,"_",j,'_fit_d0')]]<- Qform_d0
          
        } else {
          Qform_d1 <- initial_Q_t[[paste0('Q',k,"_",j,'_fit_d1')]] 
          Qform_d0 <- initial_Q_t[[paste0('Q',k,"_",j,'_fit_d0')]] 
        }
        
        logitQforms[[(k+1), (j+1)]] <- as.matrix(c(qlogis((predict.glm(Qform_d1, 
                                                                       newdata = Q_predict_data_d1)-a)/(b-a)),
                                                   qlogis((predict.glm(Qform_d0, 
                                                                       newdata = Q_predict_data_d0)-a)/(b-a))
        ))
      }
    }
  }
  
  ################# TARGETING STEP #######################
  for (k in 0:4){
    for (j in k:0){
      regimen_ind <- c(as.numeric(wideSimdata[[paste0('CA_',j)]] ==j+1), as.numeric(wideSimdata[[paste0('CA_',j)]] ==0))
      regimen_ind[is.na(regimen_ind)] <- 0
      weight <- rep(wideSimdata[[paste0('weights_',j)]],2)
      Cweight <- rep(wideSimdata[[paste0('Cweights_',j)]],2)
      Cweight_aggr <- rep(wideSimdata[[paste0('Cweights_aggr_',j)]],2)
      Cweight_treat_only <- rep(wideSimdata[[paste0('Cweights_treat_only_',j)]],2)
      Cweight_aggr_treat_only <- rep(wideSimdata[[paste0('Cweights_aggr_treat_only_',j)]],2)
      
      if (k == j){
        update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                  off = logitQforms[[(k+1),(j+1)]],
                                  intercept_0 = as.numeric(k == 0)*rep(1,dim(wideSimdata)[1]*2), 
                                  intercept_1 = as.numeric(k == 1)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_2 = as.numeric(k == 2)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_3 = as.numeric(k == 3)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_4 = as.numeric(k == 4)*rep(1,dim(wideSimdata)[1]*2),
                                  D0 = as.numeric(rep(wideSimdata$D_0, 2) == 0),
                                  D1 = as.numeric(rep(wideSimdata$D_0, 2) == 1),
                                  D2 = as.numeric(rep(wideSimdata$D_0, 2) == 2),
                                  cumA = c(rep(k+1,dim(wideSimdata)[1]), 
                                           rep(0,dim(wideSimdata)[1]))
                                 )
        
        update_data[[paste0('Y_',j,'_scaled')]] <-rep(wideSimdata[[paste0('Y_',j,'_scaled')]], 2)
        
        formula_string <- paste0(
          "Y_", j, "_scaled ~ ",
          "intercept_0 + intercept_1 + intercept_2 + intercept_3 + intercept_4 + D1 + D2 + ",
          "D0:cumA + D1:cumA + D2:cumA + ",
          "offset(off) - 1"
        )
        
        Qstar_fit <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                         data = update_data[regimen_ind == 1,],
                         weights = as.vector(weight[regimen_ind == 1]))
        
        Qstars[[(k+1),(j+1)]] <- predict.glm(Qstar_fit, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                              data = update_data[regimen_ind == 1,],
                              weights = as.vector(Cweight[regimen_ind == 1]))
        
        Qstars_cali[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali_aggr <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                   data = update_data[regimen_ind == 1,],
                                   weights = as.vector(Cweight_aggr[regimen_ind == 1]))
        
        Qstars_cali_aggr[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali_treat_only <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                         data = update_data[regimen_ind == 1,],
                                         weights = as.vector(Cweight_treat_only[regimen_ind == 1]))
        
        Qstars_cali_treat_only[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali_aggr_treat_only <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                              data = update_data[regimen_ind == 1,],
                                              weights = as.vector(Cweight_aggr_treat_only[regimen_ind == 1]))
        
        Qstars_cali_aggr_treat_only[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_aggr_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
      } else{ 
        update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                  Qk_j1star = (Qstars[[(k+1),(j+2)]]-a)/(b-a), 
                                  Qk_j1star_cali = (Qstars_cali[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_aggr = (Qstars_cali_aggr[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_treat_only = (Qstars_cali_treat_only[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_aggr_treat_only = (Qstars_cali_aggr_treat_only[[(k+1),(j+2)]]-a)/(b-a),
                                  off = logitQforms[[(k+1),(j+1)]],
                                  intercept_0 = as.numeric(k == 0)*rep(1,dim(wideSimdata)[1]*2), 
                                  intercept_1 = as.numeric(k == 1)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_2 = as.numeric(k == 2)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_3 = as.numeric(k == 3)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_4 = as.numeric(k == 4)*rep(1,dim(wideSimdata)[1]*2),
                                  D0 = as.numeric(rep(wideSimdata$D_0, 2) == 0),
                                  D1 = as.numeric(rep(wideSimdata$D_0, 2) == 1),
                                  D2 = as.numeric(rep(wideSimdata$D_0, 2) == 2),
                                  cumA = c(rep(k+1,dim(wideSimdata)[1]), 
                                           rep(0,dim(wideSimdata)[1])))
        
        Qk_j_star_fit <- glm(Qk_j1star ~ intercept_0 + intercept_1 + 
                               intercept_2 + intercept_3 + intercept_4 + 
                               D1 + D2 + D0:cumA + D1:cumA + D2:cumA + 
                               offset(off) - 1, 
                             family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(weight[regimen_ind == 1]))
        
        Qstars[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali <- glm(Qk_j1star_cali ~ intercept_0 + intercept_1 + 
                                    intercept_2 + intercept_3 + intercept_4 + 
                                    D1 + D2 + D0:cumA + D1:cumA + D2:cumA + 
                                    offset(off) - 1,
                                  family = 'quasibinomial', 
                                  data = update_data[regimen_ind == 1,],
                                  weights = as.vector(Cweight[regimen_ind == 1]))
        
        Qstars_cali[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_aggr <- glm(Qk_j1star_cali_aggr ~ intercept_0 + intercept_1 + 
                                         intercept_2 + intercept_3 + intercept_4 + 
                                         D1 + D2 + D0:cumA + D1:cumA + D2:cumA + 
                                         offset(off) - 1,
                                       family = 'quasibinomial', 
                                       data = update_data[regimen_ind == 1,],
                                       weights = as.vector(Cweight_aggr[regimen_ind == 1]))
        
        Qstars_cali_aggr[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_treat_only <- glm(Qk_j1star_cali_treat_only ~ intercept_0 + intercept_1 + 
                                               intercept_2 + intercept_3 + intercept_4 + 
                                               D1 + D2 + D0:cumA + D1:cumA + D2:cumA + 
                                               offset(off) - 1,
                                             family = 'quasibinomial', 
                                             data = update_data[regimen_ind == 1,],
                                             weights = as.vector(Cweight_treat_only[regimen_ind == 1]))
        
        Qstars_cali_treat_only[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_aggr_treat_only <- glm(Qk_j1star_cali_aggr_treat_only ~ intercept_0 + intercept_1 + 
                                                    intercept_2 + intercept_3 + intercept_4 + 
                                                    D1 + D2 + D0:cumA + D1:cumA + D2:cumA + 
                                                    offset(off) - 1,
                                                  family = 'quasibinomial', 
                                                  data = update_data[regimen_ind == 1,],
                                                  weights = as.vector(Cweight_aggr_treat_only[regimen_ind == 1]))
        
        Qstars_cali_aggr_treat_only[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_aggr_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
      }
    }
  }
  
  ################# Fit MSM ##########################
  msm_fit <- fit_msm(wideSimdata, 
                      Qstars, 
                      Qstars_cali, 
                      Qstars_cali_aggr,
                      Qstars_cali_treat_only,
                      Qstars_cali_aggr_treat_only,
                      treat_models_n,
                      censor_models_n)
  
  if (!all(is.na(initial_Q_t))){
    Q_fits <- initial_Q_t
  }
  list(MSM_estimates = msm_fit$MSM_estimates,
       counterfactual_means = msm_fit$counterfactual_means,
       Q_fits = Q_fits,
       data_with_weights = simdata,
       treat_models_n = treat_models_n, censor_models_n = censor_models_n)
}

MyLtmleMSM_hers <- function(simdata){
  
  calculate_weights <- initial_weight_calculation(simdata)
  simdata <- calculate_weights$data
  treat_models_n <- calculate_weights$treat_models_n
  censor_models_n <- calculate_weights$censor_models_n
  simdata <- weight_calibration(simdata)
  
  wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, 
                                   value.var = c("A", "Ap", "App", "CD4_1","sqrtCD4_1", "sqrtCD4_2", "viral", "viral_1", "viral_2","HIVsym","HIVsym_1", "HIVsym_2", 
                                                 "SITE2", "SITE3", "WHITE", "OTHER", "CA", "CD4","D", "weights", "Cweights", "Cweights_aggr",
                                                 "Cweights_treat_only", "Cweights_aggr_treat_only"))
  
  a = min(simdata$CD4,na.rm = TRUE); b = max(simdata$CD4, na.rm = TRUE)
  
  wideSimdata$Y_4_scaled = (wideSimdata$CD4_4-a)/(b-a)
  wideSimdata$Y_3_scaled = (wideSimdata$CD4_3-a)/(b-a)
  wideSimdata$Y_2_scaled = (wideSimdata$CD4_2-a)/(b-a)
  wideSimdata$Y_1_scaled = (wideSimdata$CD4_1-a)/(b-a)
  wideSimdata$Y_0_scaled = (wideSimdata$CD4_0-a)/(b-a)
  
  
  logitQforms <- vector("list", 5 * 5)
  dim(logitQforms) <- c(5, 5)
  
  logitQforms_cali <- vector("list", 5 * 5)
  dim(logitQforms_cali) <- c(5, 5)
  
  logitQforms_cali_aggr <- vector("list", 5 * 5)
  dim(logitQforms_cali_aggr) <- c(5, 5)
  
  logitQforms_cali_treat_only <- vector("list", 5 * 5)
  dim(logitQforms_cali_treat_only) <- c(5, 5)
  
  logitQforms_cali_aggr_treat_only <- vector("list", 5 * 5)
  dim(logitQforms_cali_aggr_treat_only) <- c(5, 5)
  
  Qstars <- vector("list", 5 * 5)
  dim(Qstars) <- c(5, 5)
  
  Qstars_cali <- vector("list", 5 * 5)
  dim(Qstars_cali) <- c(5, 5)
  
  Qstars_cali_aggr <- vector("list", 5 * 5)
  dim(Qstars_cali_aggr) <- c(5, 5)
  
  Qstars_cali_treat_only <- vector("list", 5 * 5)
  dim(Qstars_cali_treat_only) <- c(5, 5)
  
  Qstars_cali_aggr_treat_only <- vector("list", 5 * 5)
  dim(Qstars_cali_aggr_treat_only) <- c(5, 5)
  
  
  ################# TMLE #######################
  for (k in 0:4){
    for (j in k:0){
      regimen_ind <- c(as.numeric(wideSimdata[[paste0('CA_',j)]] ==j+1), as.numeric(wideSimdata[[paste0('CA_',j)]] ==0))
      regimen_ind[is.na(regimen_ind)] <- 0
      weight <- rep(wideSimdata[[paste0('weights_',j)]],2)
      Cweight <- rep(wideSimdata[[paste0('Cweights_',j)]],2)
      Cweight_aggr <- rep(wideSimdata[[paste0('Cweights_aggr_',j)]],2)
      Cweight_treat_only <- rep(wideSimdata[[paste0('Cweights_treat_only_',j)]],2)
      Cweight_aggr_treat_only <- rep(wideSimdata[[paste0('Cweights_aggr_treat_only_',j)]],2)
      
      if (k == j){
        formula_string <- paste0('CD4_',j,' ~ A_',j, '+ Ap_',j,' + App_',j,
                                 '+ sqrtCD4_1_',j,' + sqrtCD4_2_',j,
                                 ' + viral_',j,' + viral_1_',j,' + viral_2_',j,
                                 ' + HIVsym_',j,' + HIVsym_1_',j,' + HIVsym_2_',j,
                                 ' + SITE2_',j,' + SITE3_',j,' + WHITE_',j,' + OTHER_',j)
        
        predict_data <- data.frame(c(rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1])),c(rep(1, dim(wideSimdata)[1]), rep(0, dim(wideSimdata)[1])),c(rep(1, dim(wideSimdata)[1]), rep(0, dim(wideSimdata)[1])))
        colnames(predict_data) <- paste0(c('A_', 'Ap_', 'App_'),j)
        
        vars <- attr(terms(as.formula(formula_string)), "term.labels")
        
        # strip "_4" suffix
        vars_clean <- gsub(paste0("_",j,"$"), "", vars)
        for (v in vars_clean[-1:-3]) {
          predict_data[, paste0(v, "_", j)] <- rep(wideSimdata[[paste0(v, "_", j)]], 2)
        }
        
        Qfit <- glm(data = wideSimdata, formula = as.formula(formula_string))
        
        logitQforms[[(k+1),(j+1)]] <- qlogis((predict.glm(Qfit, newdata = predict_data)-a)/(b-a))
        
        
        
        update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                  off = logitQforms[[(k+1),(j+1)]],
                                  intercept_0 = as.numeric(k == 0)*rep(1,dim(wideSimdata)[1]*2), 
                                  intercept_1 = as.numeric(k == 1)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_2 = as.numeric(k == 2)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_3 = as.numeric(k == 3)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_4 = as.numeric(k == 4)*rep(1,dim(wideSimdata)[1]*2),
                                  D0 = as.numeric(rep(wideSimdata$D_0, 2) == 0),
                                  D1 = as.numeric(rep(wideSimdata$D_0, 2) == 1),
                                  D2 = as.numeric(rep(wideSimdata$D_0, 2) == 2),
                                  cumA = c(rep(k+1,dim(wideSimdata)[1]), 
                                           rep(0,dim(wideSimdata)[1])))
        
        update_data[[paste0('Y_',j,'_scaled')]] <-rep(wideSimdata[[paste0('Y_',j,'_scaled')]], 2)
        
        formula_string <- paste0(
          "Y_", j, "_scaled ~ ",
          "intercept_0 + intercept_1 + intercept_2 + intercept_3 + intercept_4 + D1 + D2 + ",
          "D0:cumA + D1:cumA + D2:cumA + ",
          "offset(off) - 1"
        )
        
        Qstar_fit <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                         data = update_data[regimen_ind == 1,],
                         weights = as.vector(weight[regimen_ind == 1]))
        
        Qstars[[(k+1),(j+1)]] <- predict.glm(Qstar_fit, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                              data = update_data[regimen_ind == 1,],
                              weights = as.vector(Cweight[regimen_ind == 1]))
        
        Qstars_cali[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali_aggr <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                   data = update_data[regimen_ind == 1,],
                                   weights = as.vector(Cweight_aggr[regimen_ind == 1]))
        
        Qstars_cali_aggr[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali_treat_only <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                         data = update_data[regimen_ind == 1,],
                                         weights = as.vector(Cweight_treat_only[regimen_ind == 1]))
        
        Qstars_cali_treat_only[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali_aggr_treat_only <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                              data = update_data[regimen_ind == 1,],
                                              weights = as.vector(Cweight_aggr_treat_only[regimen_ind == 1]))
        
        Qstars_cali_aggr_treat_only[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_aggr_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
      } else{ 
        fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                                   Qk_j1star = Qstars[[(k+1),(j+2)]], 
                                   Qk_j1star_cali = Qstars_cali[[(k+1),(j+2)]],
                                   Qk_j1star_cali_aggr = Qstars_cali_aggr[[(k+1),(j+2)]], 
                                   Qk_j1star_cali_treat_only = Qstars_cali_treat_only[[(k+1),(j+2)]],
                                   Qk_j1star_cali_aggr_treat_only = Qstars_cali_aggr_treat_only[[(k+1),(j+2)]], 
                                   A_j       = rep(wideSimdata[[paste0('A_', j)]], 2),
                                   Ap_j      = rep(wideSimdata[[paste0('Ap_', j)]], 2),
                                   App_j     = rep(wideSimdata[[paste0('App_', j)]], 2),
                                   sqrtCD4_1_j   = rep(wideSimdata[[paste0('sqrtCD4_1_', j)]], 2),
                                   sqrtCD4_2_j   = rep(wideSimdata[[paste0('sqrtCD4_2_', j)]], 2),
                                   viral_j   = rep(wideSimdata[[paste0('viral_', j)]], 2),
                                   viral_1_j = rep(wideSimdata[[paste0('viral_1_', j)]], 2),
                                   viral_2_j = rep(wideSimdata[[paste0('viral_2_', j)]], 2),
                                   HIVsym_j   = rep(wideSimdata[[paste0('HIVsym_', j)]], 2),
                                   HIVsym_1_j = rep(wideSimdata[[paste0('HIVsym_1_', j)]], 2),
                                   HIVsym_2_j = rep(wideSimdata[[paste0('HIVsym_2_', j)]], 2),
                                   SITE2_j   = rep(wideSimdata[[paste0('SITE2_', j)]], 2),
                                   SITE3_j   = rep(wideSimdata[[paste0('SITE3_', j)]], 2),
                                   WHITE_j   = rep(wideSimdata[[paste0('WHITE_', j)]], 2),
                                   OTHER_j   = rep(wideSimdata[[paste0('OTHER_', j)]], 2))
        
        Q_predict_data_d1 <- fitting_data[1:dim(wideSimdata)[1],]
        Q_predict_data_d1$A_j <- rep(1,dim(wideSimdata)[1])
        Q_predict_data_d1$Ap_j <- rep(1,dim(wideSimdata)[1])
        Q_predict_data_d1$App_j <- rep(1,dim(wideSimdata)[1])
        
        Q_predict_data_d0 <- fitting_data[1:dim(wideSimdata)[1],]
        Q_predict_data_d0$A_j <- rep(0,dim(wideSimdata)[1])
        Q_predict_data_d0$Ap_j <- rep(0,dim(wideSimdata)[1])
        Q_predict_data_d0$App_j <- rep(0,dim(wideSimdata)[1])
        
        
        
        formula_string <- "Qk_j1star ~ A_j + Ap_j + App_j +
                   sqrtCD4_1_j + sqrtCD4_2_j +
                   viral_j + viral_1_j + viral_2_j +
                   HIVsym_j + HIVsym_1_j + HIVsym_2_j +
                   SITE2_j + SITE3_j + WHITE_j + OTHER_j"
        formula_string_cali <- "Qk_j1star_cali ~  A_j + Ap_j + App_j +
                   sqrtCD4_1_j + sqrtCD4_2_j +
                   viral_j + viral_1_j + viral_2_j +
                   HIVsym_j + HIVsym_1_j + HIVsym_2_j +
                   SITE2_j + SITE3_j + WHITE_j + OTHER_j"
        formula_string_cali_aggr <- "Qk_j1star_cali_aggr ~  A_j + Ap_j + App_j +
                   sqrtCD4_1_j + sqrtCD4_2_j +
                   viral_j + viral_1_j + viral_2_j +
                   HIVsym_j + HIVsym_1_j + HIVsym_2_j +
                   SITE2_j + SITE3_j + WHITE_j + OTHER_j"
        
        formula_string_cali_treat_only <- "Qk_j1star_cali_treat_only ~  A_j + Ap_j + App_j +
                   sqrtCD4_1_j + sqrtCD4_2_j +
                   viral_j + viral_1_j + viral_2_j +
                   HIVsym_j + HIVsym_1_j + HIVsym_2_j +
                   SITE2_j + SITE3_j + WHITE_j + OTHER_j"
        formula_string_cali_aggr_treat_only <- "Qk_j1star_cali_aggr_treat_only ~  A_j + Ap_j + App_j +
                   sqrtCD4_1_j + sqrtCD4_2_j +
                   viral_j + viral_1_j + viral_2_j +
                   HIVsym_j + HIVsym_1_j + HIVsym_2_j +
                   SITE2_j + SITE3_j + WHITE_j + OTHER_j"
        
        Qform_d1 <- glm(data = fitting_data[1:dim(wideSimdata)[1],], formula = as.formula(formula_string))
        Qform_d0 <- glm(data = fitting_data[(dim(wideSimdata)[1]+1):(dim(wideSimdata)[1]*2),], formula = as.formula(formula_string))
        
        
        logitQforms[[(k+1), (j+1)]] <- as.matrix(c(qlogis((predict.glm(Qform_d1, 
                                                                       newdata = Q_predict_data_d1)-a)/(b-a)),
                                                   qlogis((predict.glm(Qform_d0, 
                                                                       newdata = Q_predict_data_d0)-a)/(b-a))
        ))
        
        Qform_d1_cali <- glm(data = fitting_data[1:dim(wideSimdata)[1],], formula = as.formula(formula_string_cali))
        Qform_d0_cali <- glm(data = fitting_data[(dim(wideSimdata)[1]+1):(dim(wideSimdata)[1]*2),], formula = as.formula(formula_string_cali))
        
        logitQforms_cali[[(k+1), (j+1)]] <- as.matrix(c(qlogis((predict.glm(Qform_d1_cali, 
                                                                            newdata = Q_predict_data_d1)-a)/(b-a)),
                                                        qlogis((predict.glm(Qform_d0_cali, 
                                                                            newdata = Q_predict_data_d0)-a)/(b-a))
        ))
        
        Qform_d1_cali_aggr <- glm(data = fitting_data[1:dim(wideSimdata)[1],], formula = as.formula(formula_string_cali_aggr))
        Qform_d0_cali_aggr <- glm(data = fitting_data[(dim(wideSimdata)[1]+1):(dim(wideSimdata)[1]*2),], formula = as.formula(formula_string_cali_aggr))
        
        logitQforms_cali_aggr[[(k+1), (j+1)]] <- as.matrix(c(qlogis((predict.glm(Qform_d1_cali_aggr, 
                                                                                 newdata = Q_predict_data_d1)-a)/(b-a)),
                                                             qlogis((predict.glm(Qform_d0_cali_aggr, 
                                                                                 newdata = Q_predict_data_d0)-a)/(b-a))
        ))
        
        Qform_d1_cali_treat_only <- glm(data = fitting_data[1:dim(wideSimdata)[1],], formula = as.formula(formula_string_cali_treat_only))
        Qform_d0_cali_treat_only <- glm(data = fitting_data[(dim(wideSimdata)[1]+1):(dim(wideSimdata)[1]*2),], formula = as.formula(formula_string_cali_treat_only))
        
        logitQforms_cali_treat_only[[(k+1), (j+1)]] <- as.matrix(c(qlogis((predict.glm(Qform_d1_cali_treat_only, 
                                                                                       newdata = Q_predict_data_d1)-a)/(b-a)),
                                                                   qlogis((predict.glm(Qform_d0_cali_treat_only, 
                                                                                       newdata = Q_predict_data_d0)-a)/(b-a))
        ))
        
        Qform_d1_cali_aggr_treat_only <- glm(data = fitting_data[1:dim(wideSimdata)[1],], formula = as.formula(formula_string_cali_aggr_treat_only))
        Qform_d0_cali_aggr_treat_only <- glm(data = fitting_data[(dim(wideSimdata)[1]+1):(dim(wideSimdata)[1]*2),], formula = as.formula(formula_string_cali_aggr_treat_only))
        
        logitQforms_cali_aggr_treat_only[[(k+1), (j+1)]] <- as.matrix(c(qlogis((predict.glm(Qform_d1_cali_aggr_treat_only, 
                                                                                            newdata = Q_predict_data_d1)-a)/(b-a)),
                                                                        qlogis((predict.glm(Qform_d0_cali_aggr_treat_only, 
                                                                                            newdata = Q_predict_data_d0)-a)/(b-a))
        ))
        
        
        update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                  Qk_j1star = (Qstars[[(k+1),(j+2)]]-a)/(b-a), 
                                  Qk_j1star_cali = (Qstars_cali[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_aggr = (Qstars_cali_aggr[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_treat_only = (Qstars_cali_treat_only[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_aggr_treat_only = (Qstars_cali_aggr_treat_only[[(k+1),(j+2)]]-a)/(b-a),
                                  off = logitQforms[[(k+1),(j+1)]],
                                  off_cali = logitQforms_cali[[(k+1),(j+1)]],
                                  off_cali_aggr = logitQforms_cali_aggr[[(k+1),(j+1)]],
                                  off_cali_treat_only = logitQforms_cali_treat_only[[(k+1),(j+1)]],
                                  off_cali_aggr_treat_only = logitQforms_cali_aggr_treat_only[[(k+1),(j+1)]],
                                  intercept_0 = as.numeric(k == 0)*rep(1,dim(wideSimdata)[1]*2), 
                                  intercept_1 = as.numeric(k == 1)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_2 = as.numeric(k == 2)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_3 = as.numeric(k == 3)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_4 = as.numeric(k == 4)*rep(1,dim(wideSimdata)[1]*2),
                                  D0 = as.numeric(rep(wideSimdata$D_0, 2) == 0),
                                  D1 = as.numeric(rep(wideSimdata$D_0, 2) == 1),
                                  D2 = as.numeric(rep(wideSimdata$D_0, 2) == 2),
                                  cumA = c(rep(k+1,dim(wideSimdata)[1]), 
                                           rep(0,dim(wideSimdata)[1])))
        
        Qk_j_star_fit <- glm(Qk_j1star ~ intercept_0 + intercept_1 + 
                               intercept_2 + intercept_3 + intercept_4 + 
                               D1 + D2 + D0:cumA + D1:cumA + D2:cumA +
                               offset(off) - 1, 
                             family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(weight[regimen_ind == 1]))
        
        Qstars[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali <- glm(Qk_j1star_cali ~ intercept_0 + intercept_1 + 
                                    intercept_2 + intercept_3 + intercept_4 + 
                                    D1 + D2 + D0:cumA + D1:cumA + D2:cumA +
                                    offset(off_cali) - 1,
                                  family = 'quasibinomial', 
                                  data = update_data[regimen_ind == 1,],
                                  weights = as.vector(Cweight[regimen_ind == 1]))
        
        Qstars_cali[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_aggr <- glm(Qk_j1star_cali_aggr ~ intercept_0 + intercept_1 + 
                                         intercept_2 + intercept_3 + intercept_4 + 
                                         D1 + D2 + D0:cumA + D1:cumA + D2:cumA +
                                         offset(off_cali_aggr) - 1,
                                       family = 'quasibinomial', 
                                       data = update_data[regimen_ind == 1,],
                                       weights = as.vector(Cweight_aggr[regimen_ind == 1]))
        
        Qstars_cali_aggr[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_treat_only <- glm(Qk_j1star_cali_treat_only ~ intercept_0 + intercept_1 + 
                                               intercept_2 + intercept_3 + intercept_4 + 
                                               D1 + D2 + D0:cumA + D1:cumA + D2:cumA +
                                               offset(off_cali_treat_only) - 1,
                                             family = 'quasibinomial', 
                                             data = update_data[regimen_ind == 1,],
                                             weights = as.vector(Cweight_treat_only[regimen_ind == 1]))
        
        Qstars_cali_treat_only[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_aggr_treat_only <- glm(Qk_j1star_cali_aggr_treat_only ~ intercept_0 + intercept_1 + 
                                                    intercept_2 + intercept_3 + intercept_4 + 
                                                    D1 + D2 + D0:cumA + D1:cumA + D2:cumA +
                                                    offset(off_cali_aggr_treat_only) - 1,
                                                  family = 'quasibinomial', 
                                                  data = update_data[regimen_ind == 1,],
                                                  weights = as.vector(Cweight_aggr_treat_only[regimen_ind == 1]))
        
        Qstars_cali_aggr_treat_only[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_aggr_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
      }
    }
  }
  ################# Fit MSM ##########################
  msm_fit <- fit_msm(wideSimdata, 
                     Qstars, 
                     Qstars_cali, 
                     Qstars_cali_aggr,
                     Qstars_cali_treat_only,
                     Qstars_cali_aggr_treat_only,
                     treat_models_n,
                     censor_models_n)
  
  list(MSM_estimates = msm_fit$MSM_estimates,
       counterfactual_means = msm_fit$counterfactual_means)
}

MyLtmleMSM_cali_boot_hers <- function(simdata, initial_Q_t = NA, refit_weights = TRUE,treat_models_n = NA, censor_models_n = NA){
  if (refit_weights){
    calculate_weights <- initial_weight_calculation(simdata)
    simdata <- calculate_weights$data
    treat_models_n <- calculate_weights$treat_models_n
    censor_models_n <- calculate_weights$censor_models_n
  } 
  simdata <- weight_calibration(simdata)
  
  wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, 
                                   value.var = c("A", "Ap", "App", "CD4_1","sqrtCD4_1", "sqrtCD4_2", "viral", "viral_1", "viral_2","HIVsym","HIVsym_1", "HIVsym_2", 
                                                 "SITE2", "SITE3", "WHITE", "OTHER", "CA", "CD4","D", "weights", "Cweights", "Cweights_aggr",
                                                 "Cweights_treat_only", "Cweights_aggr_treat_only"))
  
  a = min(simdata$CD4,na.rm = TRUE); b = max(simdata$CD4, na.rm = TRUE)
  
  wideSimdata$Y_4_scaled = (wideSimdata$CD4_4-a)/(b-a)
  wideSimdata$Y_3_scaled = (wideSimdata$CD4_3-a)/(b-a)
  wideSimdata$Y_2_scaled = (wideSimdata$CD4_2-a)/(b-a)
  wideSimdata$Y_1_scaled = (wideSimdata$CD4_1-a)/(b-a)
  wideSimdata$Y_0_scaled = (wideSimdata$CD4_0-a)/(b-a)
  
  if (all(is.na(initial_Q_t))){
    Q_fits <- list()
  }
  
  logitQforms <- vector("list", 5 * 5)
  dim(logitQforms) <- c(5, 5)
  
  Qstars <- vector("list", 5 * 5)
  dim(Qstars) <- c(5, 5)
  
  Qstars_cali <- vector("list", 5 * 5)
  dim(Qstars_cali) <- c(5, 5)
  
  Qstars_cali_aggr <- vector("list", 5 * 5)
  dim(Qstars_cali_aggr) <- c(5, 5)
  
  Qstars_cali_treat_only <- vector("list", 5 * 5)
  dim(Qstars_cali_treat_only) <- c(5, 5)
  
  Qstars_cali_aggr_treat_only <- vector("list", 5 * 5)
  dim(Qstars_cali_aggr_treat_only) <- c(5, 5)
  ########### Initial Q fits ##########
  # First loop: get the Q forms
  
  
  for (k in 0:4){
    for (j in k:0){
      if (k == j){
        formula_string <- paste0('CD4_',j,' ~ A_',j, '+ Ap_',j,' + App_',j,
                                 '+ sqrtCD4_1_',j,' + sqrtCD4_2_',j,
                                 ' + viral_',j,' + viral_1_',j,' + viral_2_',j,
                                 ' + HIVsym_',j,' + HIVsym_1_',j,' + HIVsym_2_',j,
                                 ' + SITE2_',j,' + SITE3_',j,' + WHITE_',j,' + OTHER_',j)
        
        predict_data <- data.frame(c(rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1])),c(rep(1, dim(wideSimdata)[1]), rep(0, dim(wideSimdata)[1])),c(rep(1, dim(wideSimdata)[1]), rep(0, dim(wideSimdata)[1])))
        colnames(predict_data) <- paste0(c('A_', 'Ap_', 'App_'),j)
        
        vars <- attr(terms(as.formula(formula_string)), "term.labels")
        
        # strip "_4" suffix
        vars_clean <- gsub(paste0("_",j,"$"), "", vars)
        for (v in vars_clean[-1:-3]) {
          predict_data[, paste0(v, "_", j)] <- rep(wideSimdata[[paste0(v, "_", j)]], 2)
        }
        
        if (all(is.na(initial_Q_t))){
          
          Qfit <- glm(data = wideSimdata, formula = as.formula(formula_string))
          Q_fits[[paste0('Q',k,"_",j,'_fit')]] <- Qfit
        } else {
          Qfit <- initial_Q_t[[paste0('Q',k,"_",j,'_fit')]]
        }
        
        logitQforms[[(k+1),(j+1)]] <- qlogis((predict.glm(Qfit, newdata = predict_data)-a)/(b-a))
        
        
      } else{ 
        fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                                   Qk_j1 = plogis(logitQforms[[(k+1),(j+2)]])*(b-a)+a, 
                                   A_j       = rep(wideSimdata[[paste0('A_', j)]], 2),
                                   Ap_j      = rep(wideSimdata[[paste0('Ap_', j)]], 2),
                                   App_j     = rep(wideSimdata[[paste0('App_', j)]], 2),
                                   sqrtCD4_1_j   = rep(wideSimdata[[paste0('sqrtCD4_1_', j)]], 2),
                                   sqrtCD4_2_j   = rep(wideSimdata[[paste0('sqrtCD4_2_', j)]], 2),
                                   viral_j   = rep(wideSimdata[[paste0('viral_', j)]], 2),
                                   viral_1_j = rep(wideSimdata[[paste0('viral_1_', j)]], 2),
                                   viral_2_j = rep(wideSimdata[[paste0('viral_2_', j)]], 2),
                                   HIVsym_j   = rep(wideSimdata[[paste0('HIVsym_', j)]], 2),
                                   HIVsym_1_j = rep(wideSimdata[[paste0('HIVsym_1_', j)]], 2),
                                   HIVsym_2_j = rep(wideSimdata[[paste0('HIVsym_2_', j)]], 2),
                                   SITE2_j   = rep(wideSimdata[[paste0('SITE2_', j)]], 2),
                                   SITE3_j   = rep(wideSimdata[[paste0('SITE3_', j)]], 2),
                                   WHITE_j   = rep(wideSimdata[[paste0('WHITE_', j)]], 2),
                                   OTHER_j   = rep(wideSimdata[[paste0('OTHER_', j)]], 2))
        
        Q_predict_data_d1 <- fitting_data[1:dim(wideSimdata)[1],]
        Q_predict_data_d1$A_j <- rep(1,dim(wideSimdata)[1])
        Q_predict_data_d1$Ap_j <- rep(1,dim(wideSimdata)[1])
        Q_predict_data_d1$App_j <- rep(1,dim(wideSimdata)[1])
        
        Q_predict_data_d0 <- fitting_data[1:dim(wideSimdata)[1],]
        Q_predict_data_d0$A_j <- rep(0,dim(wideSimdata)[1])
        Q_predict_data_d0$Ap_j <- rep(0,dim(wideSimdata)[1])
        Q_predict_data_d0$App_j <- rep(0,dim(wideSimdata)[1])
        
        
        if (all(is.na(initial_Q_t))){
          formula_string <- paste0("Qk_j1 ~ A_j + Ap_j + App_j +",
                                   "sqrtCD4_1_j + sqrtCD4_2_j +",
                                   "viral_j + viral_1_j + viral_2_j +",
                                   "HIVsym_j + HIVsym_1_j + HIVsym_2_j +",
                                   "SITE2_j + SITE3_j + WHITE_j + OTHER_j") 
          
          Qform_d1 <- glm(data = fitting_data[1:dim(wideSimdata)[1],], formula = as.formula(formula_string))
          Qform_d0 <- glm(data = fitting_data[(dim(wideSimdata)[1]+1):(dim(wideSimdata)[1]*2),], formula = as.formula(formula_string))
          
          Q_fits[[paste0('Q',k,"_",j,'_fit_d1')]]<- Qform_d1
          Q_fits[[paste0('Q',k,"_",j,'_fit_d0')]]<- Qform_d0
          
        } else {
          Qform_d1 <- initial_Q_t[[paste0('Q',k,"_",j,'_fit_d1')]] 
          Qform_d0 <- initial_Q_t[[paste0('Q',k,"_",j,'_fit_d0')]] 
        }
        
        logitQforms[[(k+1), (j+1)]] <- as.matrix(c(qlogis((predict.glm(Qform_d1, 
                                                                       newdata = Q_predict_data_d1)-a)/(b-a)),
                                                   qlogis((predict.glm(Qform_d0, 
                                                                       newdata = Q_predict_data_d0)-a)/(b-a))
        ))
      }
    }
  }
  
  ################# TARGETING STEP #######################
  for (k in 0:4){
    for (j in k:0){
      regimen_ind <- c(as.numeric(wideSimdata[[paste0('CA_',j)]] ==j+1), as.numeric(wideSimdata[[paste0('CA_',j)]] ==0))
      regimen_ind[is.na(regimen_ind)] <- 0
      weight <- rep(wideSimdata[[paste0('weights_',j)]],2)
      Cweight <- rep(wideSimdata[[paste0('Cweights_',j)]],2)
      Cweight_aggr <- rep(wideSimdata[[paste0('Cweights_aggr_',j)]],2)
      Cweight_treat_only <- rep(wideSimdata[[paste0('Cweights_treat_only_',j)]],2)
      Cweight_aggr_treat_only <- rep(wideSimdata[[paste0('Cweights_aggr_treat_only_',j)]],2)
      
      if (k == j){
        update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                  off = logitQforms[[(k+1),(j+1)]],
                                  intercept_0 = as.numeric(k == 0)*rep(1,dim(wideSimdata)[1]*2), 
                                  intercept_1 = as.numeric(k == 1)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_2 = as.numeric(k == 2)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_3 = as.numeric(k == 3)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_4 = as.numeric(k == 4)*rep(1,dim(wideSimdata)[1]*2),
                                  D0 = as.numeric(rep(wideSimdata$D_0, 2) == 0),
                                  D1 = as.numeric(rep(wideSimdata$D_0, 2) == 1),
                                  D2 = as.numeric(rep(wideSimdata$D_0, 2) == 2),
                                  cumA = c(rep(k+1,dim(wideSimdata)[1]), 
                                           rep(0,dim(wideSimdata)[1])))
        
        update_data[[paste0('Y_',j,'_scaled')]] <-rep(wideSimdata[[paste0('Y_',j,'_scaled')]], 2)
        
        formula_string <- paste0(
          "Y_", j, "_scaled ~ ",
          "intercept_0 + intercept_1 + intercept_2 + intercept_3 + intercept_4 + D1 + D2 + ",
          "D0:cumA + D1:cumA + D2:cumA + ",
          "offset(off) - 1"
        )
        
        Qstar_fit <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                         data = update_data[regimen_ind == 1,],
                         weights = as.vector(weight[regimen_ind == 1]))
        
        Qstars[[(k+1),(j+1)]] <- predict.glm(Qstar_fit, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                              data = update_data[regimen_ind == 1,],
                              weights = as.vector(Cweight[regimen_ind == 1]))
        
        Qstars_cali[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali_aggr <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                   data = update_data[regimen_ind == 1,],
                                   weights = as.vector(Cweight_aggr[regimen_ind == 1]))
        
        Qstars_cali_aggr[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali_treat_only <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                         data = update_data[regimen_ind == 1,],
                                         weights = as.vector(Cweight_treat_only[regimen_ind == 1]))
        
        Qstars_cali_treat_only[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
        Qstar_fit_cali_aggr_treat_only <- glm(formula = as.formula(formula_string), family = 'quasibinomial', 
                                              data = update_data[regimen_ind == 1,],
                                              weights = as.vector(Cweight_aggr_treat_only[regimen_ind == 1]))
        
        Qstars_cali_aggr_treat_only[[(k+1),(j+1)]] <- predict.glm(Qstar_fit_cali_aggr_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
      } else{ 
        update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                  Qk_j1star = (Qstars[[(k+1),(j+2)]]-a)/(b-a), 
                                  Qk_j1star_cali = (Qstars_cali[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_aggr = (Qstars_cali_aggr[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_treat_only = (Qstars_cali_treat_only[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_aggr_treat_only = (Qstars_cali_aggr_treat_only[[(k+1),(j+2)]]-a)/(b-a),
                                  off = logitQforms[[(k+1),(j+1)]],
                                  intercept_0 = as.numeric(k == 0)*rep(1,dim(wideSimdata)[1]*2), 
                                  intercept_1 = as.numeric(k == 1)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_2 = as.numeric(k == 2)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_3 = as.numeric(k == 3)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_4 = as.numeric(k == 4)*rep(1,dim(wideSimdata)[1]*2),
                                  D0 = as.numeric(rep(wideSimdata$D_0, 2) == 0),
                                  D1 = as.numeric(rep(wideSimdata$D_0, 2) == 1),
                                  D2 = as.numeric(rep(wideSimdata$D_0, 2) == 2),
                                  cumA = c(rep(k+1,dim(wideSimdata)[1]), 
                                           rep(0,dim(wideSimdata)[1])))
        
        Qk_j_star_fit <- glm(Qk_j1star ~ intercept_0 + intercept_1 + 
                               intercept_2 + intercept_3 + intercept_4 + 
                               D1 + D2 + D0:cumA + D1:cumA + D2:cumA +
                               offset(off) - 1, 
                             family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(weight[regimen_ind == 1]))
        
        Qstars[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali <- glm(Qk_j1star_cali ~ intercept_0 + intercept_1 + 
                                    intercept_2 + intercept_3 + intercept_4 + 
                                    D1 + D2 + D0:cumA + D1:cumA + D2:cumA + 
                                    offset(off) - 1,
                                  family = 'quasibinomial', 
                                  data = update_data[regimen_ind == 1,],
                                  weights = as.vector(Cweight[regimen_ind == 1]))
        
        Qstars_cali[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_aggr <- glm(Qk_j1star_cali_aggr ~ intercept_0 + intercept_1 + 
                                         intercept_2 + intercept_3 + intercept_4 + 
                                         D1 + D2 + D0:cumA + D1:cumA + D2:cumA +
                                         offset(off) - 1,
                                       family = 'quasibinomial', 
                                       data = update_data[regimen_ind == 1,],
                                       weights = as.vector(Cweight_aggr[regimen_ind == 1]))
        
        Qstars_cali_aggr[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_treat_only <- glm(Qk_j1star_cali_treat_only ~ intercept_0 + intercept_1 + 
                                               intercept_2 + intercept_3 + intercept_4 + 
                                               D1 + D2 + D0:cumA + D1:cumA + D2:cumA +
                                               offset(off) - 1,
                                             family = 'quasibinomial', 
                                             data = update_data[regimen_ind == 1,],
                                             weights = as.vector(Cweight_treat_only[regimen_ind == 1]))
        
        Qstars_cali_treat_only[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_aggr_treat_only <- glm(Qk_j1star_cali_aggr_treat_only ~intercept_0 + intercept_1 + 
                                                    intercept_2 + intercept_3 + intercept_4 + 
                                                    D1 + D2 + D0:cumA + D1:cumA + D2:cumA +
                                                    offset(off) - 1,
                                                  family = 'quasibinomial', 
                                                  data = update_data[regimen_ind == 1,],
                                                  weights = as.vector(Cweight_aggr_treat_only[regimen_ind == 1]))
        
        Qstars_cali_aggr_treat_only[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_aggr_treat_only, newdata = update_data, type = 'response')*(b-a) +a
        
      }
    }
  }
  ################# Fit MSM ##########################
  msm_fit <- fit_msm(wideSimdata, 
                     Qstars, 
                     Qstars_cali, 
                     Qstars_cali_aggr,
                     Qstars_cali_treat_only,
                     Qstars_cali_aggr_treat_only,
                     treat_models_n,
                     censor_models_n)
  
  if (!all(is.na(initial_Q_t))){
    Q_fits <- initial_Q_t
  }
  list(MSM_estimates = msm_fit$MSM_estimates,
       counterfactual_means = msm_fit$counterfactual_means,
       Q_fits = Q_fits,
       data_with_weights = simdata,
       treat_models_n = treat_models_n, censor_models_n = censor_models_n)
}