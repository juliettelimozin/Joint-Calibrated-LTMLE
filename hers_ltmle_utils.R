MyLtmleMSM_modified_hers <- function(simdata, initial_Q_t = NA, refit_weights = TRUE){
  
  treat_model_A_0 <- glm(A~ CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER-1, data = simdata[simdata$t == 0,], family = 'binomial')
  treat_model_A_1 <- glm(A~ Ap+CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 1,], family = 'binomial')
  treat_model_A_2 <- glm(A~ Ap+CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 2,], family = 'binomial')
  treat_model_A_3 <- glm(A~ Ap+CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 3,], family = 'binomial')
  treat_model_A_4 <- glm(A~ Ap+CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 4,], family = 'binomial')
  
  
  simdata$ps <- 1.0
  simdata[simdata$t == 0& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 1,], type = 'response'))
  simdata[simdata$t == 1& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 1,], type = 'response'))
  simdata[simdata$t == 2& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 1,], type = 'response'))
  simdata[simdata$t == 3& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_3, newdata = simdata[simdata$t == 3& simdata$A == 1,], type = 'response'))
  simdata[simdata$t == 4& simdata$A == 1,]$ps = as.vector(predict.glm(treat_model_A_4, newdata = simdata[simdata$t == 4& simdata$A == 1,], type = 'response'))
  
  simdata[simdata$t == 0& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_0, newdata = simdata[simdata$t == 0& simdata$A == 0,], type = 'response'))
  simdata[simdata$t == 1& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_1, newdata = simdata[simdata$t == 1& simdata$A == 0,], type = 'response'))
  simdata[simdata$t == 2& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_2, newdata = simdata[simdata$t == 2& simdata$A == 0,], type = 'response'))
  simdata[simdata$t == 3& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_3, newdata = simdata[simdata$t == 3& simdata$A == 0,], type = 'response'))
  simdata[simdata$t == 4& simdata$A == 0,]$ps = as.vector(1-predict.glm(treat_model_A_4, newdata = simdata[simdata$t == 4& simdata$A == 0,], type = 'response'))
  
  censor_model_0 <- glm(C ~ CD4_1 + CD4  + viral + viral_1  + HIVsym + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 0,], family = 'binomial')
  censor_model_1 <- glm(C ~ CD4_1 + CD4  + viral + viral_1  + HIVsym + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 1,], family = 'binomial')
  censor_model_2 <- glm(C ~ CD4_1 + CD4  + viral + viral_1  + HIVsym + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 2,], family = 'binomial')
  censor_model_3 <- glm(C ~ CD4_1 + CD4  + viral + viral_1  + HIVsym + HIVsym_1+ SITE2 + SITE3 + WHITE + OTHER, data = simdata[simdata$t == 3,], family = 'binomial')
  
  simdata$lagCD4_1 <- simdata$CD4_2
  simdata$lagCD4 <- simdata$CD4_1
  simdata$lagviral <- simdata$viral_1
  simdata$lagviral_1 <- simdata$viral_2
  simdata$lagHIVsym <- simdata$HIVsym_1
  simdata$lagHIVsym_1 <- simdata$HIVsym_2
  simdata$lagSITE2 <- simdata$SITE2
  simdata$lagSITE3 <- simdata$SITE3
  simdata$lagWHITE <- simdata$WHITE
  simdata$lagOTHER <- simdata$OTHER
  
  
  simdata$pr_c <- 1.0
  simdata[simdata$t == 1,]$pr_c <- as.vector(1-predict.glm(censor_model_0, 
                                                           newdata = data.frame(CD4_1 = simdata[simdata$t == 1,]$lagCD4_1,
                                                                                CD4 = simdata[simdata$t == 1,]$lagCD4,
                                                                                viral = simdata[simdata$t == 1,]$lagviral,
                                                                                viral_1 = simdata[simdata$t == 1,]$lagviral_1,
                                                                                HIVsym = simdata[simdata$t == 1,]$lagHIVsym,
                                                                                HIVsym_1 = simdata[simdata$t == 1,]$lagHIVsym_1,
                                                                                SITE2 = simdata[simdata$t == 1,]$SITE2,
                                                                                SITE3 = simdata[simdata$t == 1,]$SITE3,
                                                                                WHITE = simdata[simdata$t == 1,]$WHITE,
                                                                                OTHER = simdata[simdata$t == 1,]$OTHER),
                                                           type = 'response'))
  simdata[simdata$t == 2,]$pr_c <- as.vector(1-predict.glm(censor_model_1,
                                                           newdata = data.frame(CD4_1 = simdata[simdata$t == 2,]$lagCD4_1,
                                                                                CD4 = simdata[simdata$t == 2,]$lagCD4,
                                                                                viral = simdata[simdata$t == 2,]$lagviral,
                                                                                viral_1 = simdata[simdata$t == 2,]$lagviral_1,
                                                                                HIVsym = simdata[simdata$t == 2,]$lagHIVsym,
                                                                                HIVsym_1 = simdata[simdata$t == 2,]$lagHIVsym_1,
                                                                                SITE2 = simdata[simdata$t == 2,]$SITE2,
                                                                                SITE3 = simdata[simdata$t == 2,]$SITE3,
                                                                                WHITE = simdata[simdata$t == 2,]$WHITE,
                                                                                OTHER = simdata[simdata$t == 2,]$OTHER),
                                                           type = 'response'))
  simdata[simdata$t == 3,]$pr_c <- as.vector(1-predict.glm(censor_model_2,
                                                           newdata = data.frame(CD4_1 = simdata[simdata$t == 3,]$lagCD4_1,
                                                                                CD4 = simdata[simdata$t == 3,]$lagCD4,
                                                                                viral = simdata[simdata$t == 3,]$lagviral,
                                                                                viral_1 = simdata[simdata$t == 3,]$lagviral_1,
                                                                                HIVsym = simdata[simdata$t == 3,]$lagHIVsym,
                                                                                HIVsym_1 = simdata[simdata$t == 3,]$lagHIVsym_1,
                                                                                SITE2 = simdata[simdata$t == 3,]$SITE2,
                                                                                SITE3 = simdata[simdata$t == 3,]$SITE3,
                                                                                WHITE = simdata[simdata$t == 3,]$WHITE,
                                                                                OTHER = simdata[simdata$t == 3,]$OTHER),
                                                           type = 'response'))
  simdata[simdata$t == 4,]$pr_c <- as.vector(1-predict.glm(censor_model_3,
                                                           newdata = data.frame(CD4_1 = simdata[simdata$t == 4,]$lagCD4_1,
                                                                                CD4 = simdata[simdata$t == 4,]$lagCD4,
                                                                                viral = simdata[simdata$t == 4,]$lagviral,
                                                                                viral_1 = simdata[simdata$t == 4,]$lagviral_1,
                                                                                HIVsym = simdata[simdata$t == 4,]$lagHIVsym,
                                                                                HIVsym_1 = simdata[simdata$t == 4,]$lagHIVsym_1,
                                                                                SITE2 = simdata[simdata$t == 4,]$SITE2,
                                                                                SITE3 = simdata[simdata$t == 4,]$SITE3,
                                                                                WHITE = simdata[simdata$t == 4,]$WHITE,
                                                                                OTHER = simdata[simdata$t == 4,]$OTHER),
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
  simdata[simdata$t == 3 & !(simdata$CA == 4),]$RA <- 0
  simdata[simdata$t == 4 & !(simdata$CA == 5),]$RA <- 0
  
  simdata$tCD4_1     <- simdata$t * simdata$CD4_1
  simdata$tCD4_2     <- simdata$t * simdata$CD4_2
  simdata$tCD4       <- simdata$t * simdata$CD4
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
                                                                var = c("CD4_1","CD4_2","viral_1","viral_2","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"), 
                                                                censor = TRUE, 
                                                                c_var = c("CD4_1","CD4","viral","viral_1","HIVsym","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"
                                                                ),
                                                                weights_var = 'weight')
  calibrate_always_treated_aggr <- aggregated_calibration_from_baseline(simdata, 
                                                                        var = c(
                                                                          "CD4_1", "CD4_2", "viral_1", "viral_2", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER", 
                                                                          "tCD4_1", "tCD4_2", "tviral_1", "tviral_2", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                        ), 
                                                                        censor = TRUE, 
                                                                        c_var = c(
                                                                          "CD4_1", "CD4", "viral", "viral_1", "HIVsym", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER",
                                                                          "tCD4_1", "tCD4", "tviral", "tviral_1", "tHIVsym", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                        ),
                                                                        weights_var = 'weight')
  
  
  simdata$Cweights <- calibrate_always_treated$data$Cweights
  simdata$Cweights_aggr <- calibrate_always_treated_aggr$data$Cweights
  
  simdata$RA <- 1
  simdata[simdata$t == 0 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 1 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 2 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 3 & !(simdata$CA == 0),]$RA <- 0
  simdata[simdata$t == 4 & !(simdata$CA == 0),]$RA <- 0
  
  calibrate_never_treated <- calibration_by_time_from_baseline(simdata, 
                                                               var = c("CD4_1","CD4_2","viral_1","viral_2","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"), 
                                                               censor = TRUE, 
                                                               c_var = c("CD4_1","CD4","viral","viral_1","HIVsym","HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER"
                                                               ),
                                                               weights_var = 'Cweights')
  calibrate_never_treated_aggr <- aggregated_calibration_from_baseline(simdata, 
                                                                       var = c(
                                                                         "CD4_1", "CD4_2", "viral_1", "viral_2", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER", 
                                                                         "tCD4_1", "tCD4_2", "tviral_1", "tviral_2", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                       ), 
                                                                       censor = TRUE, 
                                                                       c_var = c(
                                                                         "CD4_1", "CD4", "viral", "viral_1", "HIVsym", "HIVsym_1", "SITE2", "SITE3", "WHITE", "OTHER", 
                                                                         "tCD4_1", "tCD4", "tviral", "tviral_1", "tHIVsym", "tHIVsym_1","tSITE2", "tSITE3", "tWHITE", "tOTHER"
                                                                       ),
                                                                       weights_var = 'Cweights_aggr')
  
  
  simdata$Cweights <- calibrate_never_treated$data$Cweights
  simdata$Cweights_aggr <- calibrate_never_treated_aggr$data$Cweights
  
  simdata$weights <- simdata$weight
  
  wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, 
                                   value.var = c("A", "Ap", "App", "CD4_1", "CD4_2", "viral", "viral_1", "viral_2","HIVsym","HIVsym_1", "HIVsym_2", 
                                                 "SITE2", "SITE3", "WHITE", "OTHER", "CA", "CD4", "weights", "Cweights", "Cweights_aggr"))
  
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
  ########### Initial Q fits ##########
  # First loop: get the Q forms
  
  
  for (k in 0:4){
    for (j in k:0){
      if (k == j){
        formula_string <- paste0('CD4_',j,' ~ A_',j, '+ Ap_',j,' + App_',j,
                                 '+ CD4_1_',j,' + CD4_2_',j,
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
                                   CD4_1_j   = rep(wideSimdata[[paste0('CD4_1_', j)]], 2),
                                   CD4_2_j   = rep(wideSimdata[[paste0('CD4_2_', j)]], 2),
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
          formula_string <- "Qk_j1 ~ A_j + Ap_j + App_j +
                   CD4_1_j + CD4_2_j +
                   viral_j + viral_1_j + viral_2_j +
                   HIVsym_j + HIVsym_1_j + HIVsym_2_j +
                   SITE2_j + SITE3_j + WHITE_j + OTHER_j"
          
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
      
      if (k == j){
        update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                  off = logitQforms[[(k+1),(j+1)]],
                                  intercept_0 = as.numeric(k == 0)*rep(1,dim(wideSimdata)[1]*2), 
                                  intercept_1 = as.numeric(k == 1)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_2 = as.numeric(k == 2)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_3 = as.numeric(k == 3)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_4 = as.numeric(k == 4)*rep(1,dim(wideSimdata)[1]*2),
                                  strata_1 = rep(as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58),2),
                                  strata_2 = rep(as.numeric(wideSimdata$CD4_1_0 > 0.58),2),
                                  viral_1 = rep(wideSimdata$viral_1_0,2),
                                  HIVsym_1 = rep(wideSimdata$HIVsym_1_0,2),
                                  Ap = rep(wideSimdata$Ap_0,2),
                                  strata_0cumA = c(as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(k+1,dim(wideSimdata)[1]), 
                                                   as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(0,dim(wideSimdata)[1])),
                                  strata_1cumA = c(as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58)*rep(k+1,dim(wideSimdata)[1]),
                                                   as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58)*rep(0,dim(wideSimdata)[1])),
                                  strata_2cumA = c(as.numeric(wideSimdata$CD4_1_0 > 0.58)*rep(k+1,dim(wideSimdata)[1]), 
                                                   as.numeric(wideSimdata$CD4_1_0 > 0.58)*rep(0,dim(wideSimdata)[1])))
        
        update_data[[paste0('Y_',j,'_scaled')]] <-rep(wideSimdata[[paste0('Y_',j,'_scaled')]], 2)
        
        formula_string <- paste0(
          "Y_", j, "_scaled ~ ",
          "intercept_0 + intercept_1 + intercept_2 + intercept_3 + intercept_4 + strata_1 + strata_2 + viral_1 + HIVsym_1 + Ap + ",
          "strata_0cumA + strata_1cumA + strata_2cumA + ",
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
        
      } else{ 
        update_data <- data.frame(id = rep(wideSimdata$ID,2),
                                  Qk_j1star = (Qstars[[(k+1),(j+2)]]-a)/(b-a), 
                                  Qk_j1star_cali = (Qstars_cali[[(k+1),(j+2)]]-a)/(b-a),
                                  Qk_j1star_cali_aggr = (Qstars_cali_aggr[[(k+1),(j+2)]]-a)/(b-a),
                                  off = logitQforms[[(k+1),(j+1)]],
                                  intercept_0 = as.numeric(k == 0)*rep(1,dim(wideSimdata)[1]*2), 
                                  intercept_1 = as.numeric(k == 1)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_2 = as.numeric(k == 2)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_3 = as.numeric(k == 3)*rep(1,dim(wideSimdata)[1]*2),
                                  intercept_4 = as.numeric(k == 4)*rep(1,dim(wideSimdata)[1]*2),
                                  strata_1 = rep(as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58),2),
                                  strata_2 = rep(as.numeric(wideSimdata$CD4_1_0 > 0.58),2),
                                  viral_1 = rep(wideSimdata$viral_1_0,2),
                                  HIVsym_1 = rep(wideSimdata$HIVsym_1_0,2),
                                  Ap = rep(wideSimdata$Ap_0,2),
                                  strata_0cumA = c(as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(k+1,dim(wideSimdata)[1]), 
                                                   as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(0,dim(wideSimdata)[1])),
                                  strata_1cumA = c(as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58)*rep(k+1,dim(wideSimdata)[1]),
                                                   as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58)*rep(0,dim(wideSimdata)[1])),
                                  strata_2cumA = c(as.numeric(wideSimdata$CD4_1_0 > 0.58)*rep(k+1,dim(wideSimdata)[1]), 
                                                   as.numeric(wideSimdata$CD4_1_0 > 0.58)*rep(0,dim(wideSimdata)[1])))
        
        Qk_j_star_fit <- glm(Qk_j1star ~ intercept_0 + intercept_1 + 
                               intercept_2 + intercept_3 + intercept_4 + 
                               strata_1 + strata_2 + viral_1 + HIVsym_1 + 
                               Ap + strata_0cumA + strata_1cumA + strata_2cumA + 
                               offset(off) - 1, 
                             family = 'quasibinomial', 
                             data = update_data[regimen_ind == 1,],
                             weights = as.vector(weight[regimen_ind == 1]))
        
        Qstars[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali <- glm(Qk_j1star_cali ~ intercept_0 + intercept_1 + 
                                    intercept_2 + intercept_3 + intercept_4 + 
                                    strata_1 + strata_2 + viral_1 + HIVsym_1 + 
                                    Ap + strata_0cumA + strata_1cumA + strata_2cumA + 
                                    offset(off) - 1,
                                  family = 'quasibinomial', 
                                  data = update_data[regimen_ind == 1,],
                                  weights = as.vector(Cweight[regimen_ind == 1]))
        
        Qstars_cali[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali, newdata = update_data, type = 'response')*(b-a) +a
        
        Qk_j_star_fit_cali_aggr <- glm(Qk_j1star_cali_aggr ~ intercept_0 + intercept_1 + 
                                         intercept_2 + intercept_3 + intercept_4 + 
                                         strata_1 + strata_2 + viral_1 + HIVsym_1 + 
                                         Ap + strata_0cumA + strata_1cumA + strata_2cumA + 
                                         offset(off) - 1,
                                       family = 'quasibinomial', 
                                       data = update_data[regimen_ind == 1,],
                                       weights = as.vector(Cweight_aggr[regimen_ind == 1]))
        
        Qstars_cali_aggr[[(k+1),(j+1)]] <- predict.glm(Qk_j_star_fit_cali_aggr, newdata = update_data, type = 'response')*(b-a) +a
        
      }
    }
  }
  ################# Fit MSM ##########################
  
  msm_fitting_data_transformed <- data.frame(id = rep(wideSimdata$ID,10), 
                                             t = c(rep(0,dim(wideSimdata)[1]*2), 
                                                   rep(1,dim(wideSimdata)[1]*2), 
                                                   rep(2,dim(wideSimdata)[1]*2),
                                                   rep(3,dim(wideSimdata)[1]*2),
                                                   rep(4,dim(wideSimdata)[1]*2)),
                                             Y =  unlist(Qstars[, 1]), 
                                             Y_cali =  unlist(Qstars_cali[, 1]),
                                             Y_cali_aggr =  unlist(Qstars_cali_aggr[, 1]),
                                             cumA = c(rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                      rep(2,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                      rep(3,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                      rep(4,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                      rep(5,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1])),
                                             intercept_0 = as.numeric(t == 0),
                                             intercept_1 = as.numeric(t == 1),
                                             intercept_2 = as.numeric(t == 2),
                                             intercept_3 = as.numeric(t == 3),
                                             intercept_4 = as.numeric(t == 4),
                                             strata_1 = rep(as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58),10),
                                             strata_2 = rep(as.numeric(wideSimdata$CD4_1_0 > 0.58),10),
                                             viral_1 = rep(wideSimdata$viral_1_0,10),
                                             HIVsym_1 = rep(wideSimdata$HIVsym_1_0,10),
                                             Ap = rep(wideSimdata$Ap_0,10),
                                             strata_0cumA = c(as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(1,dim(wideSimdata)[1]), as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(0,dim(wideSimdata)[1]), 
                                                              as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(2,dim(wideSimdata)[1]), as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(0,dim(wideSimdata)[1]), 
                                                              as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(3,dim(wideSimdata)[1]), as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(0,dim(wideSimdata)[1]),
                                                              as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(4,dim(wideSimdata)[1]), as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(0,dim(wideSimdata)[1]),
                                                              as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(5,dim(wideSimdata)[1]), as.numeric(-0.61 > wideSimdata$CD4_1_0)*rep(0,dim(wideSimdata)[1])),
                                             strata_1cumA =  c(as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58)*rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                               as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58)*rep(2,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                               as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58)*rep(3,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                               as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58)*rep(4,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                               as.numeric(-0.61 <= wideSimdata$CD4_1_0 & wideSimdata$CD4_1_0 <= 0.58)*rep(5,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1])),
                                             strata_2cumA = c(as.numeric(wideSimdata$CD4_1_0 > 0.58)*rep(1,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                              as.numeric(wideSimdata$CD4_1_0 > 0.58)*rep(2,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]), 
                                                              as.numeric(wideSimdata$CD4_1_0 > 0.58)*rep(3,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                              as.numeric(wideSimdata$CD4_1_0 > 0.58)*rep(4,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1]),
                                                              as.numeric(wideSimdata$CD4_1_0 > 0.58)*rep(5,dim(wideSimdata)[1]), rep(0,dim(wideSimdata)[1])))
  
  msm_transformed <- glm(Y ~ intercept_0 + intercept_1 + 
                           intercept_2 + intercept_3 + intercept_4 + 
                           strata_1 + strata_2 + viral_1 + HIVsym_1 + 
                           Ap + strata_0cumA + strata_1cumA + strata_2cumA  - 1, data = msm_fitting_data_transformed)
  msm_transformed_cali <- glm(Y_cali ~ intercept_0 + intercept_1 + 
                                intercept_2 + intercept_3 + intercept_4 + 
                                strata_1 + strata_2 + viral_1 + HIVsym_1 + 
                                Ap + strata_0cumA + strata_1cumA + strata_2cumA  - 1, data = msm_fitting_data_transformed)
  msm_transformed_cali_aggr <- glm(Y_cali_aggr ~ intercept_0 + intercept_1 + 
                                     intercept_2 + intercept_3 + intercept_4 + 
                                     strata_1 + strata_2 + viral_1 + HIVsym_1 + 
                                     Ap + strata_0cumA + strata_1cumA + strata_2cumA  - 1, data = msm_fitting_data_transformed)
  
  
  estimates <- c(msm_transformed$coefficients,
                 msm_transformed_cali$coefficients,
                 msm_transformed_cali_aggr$coefficients)
  
  dim(estimates) <- c(13,3)
  
  names(estimates) <- c('MLE LTML',
                        'CMLE LTMLE',
                        'Aggr. CMLE LTMLE')
  
  
  list(estimates = estimates,
       Q_fits = Q_fits)
}

