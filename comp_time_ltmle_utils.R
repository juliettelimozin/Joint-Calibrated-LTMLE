MyLtmleMSM <- function(simdata){
  
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
  
  # censor_model_0 <- glm(C ~ X1 + X2 + X3 + X4, data = simdata[simdata$t == 0,], family = 'binomial')
  # censor_model_1 <- glm(C ~ X1 + X2 + X3 + X4, data = simdata[simdata$t == 1,], family = 'binomial')
  # 
  # simdata$lagX1 <- ave(simdata$X1, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
  # simdata$lagX2 <- ave(simdata$X2, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
  # simdata$lagX3 <- ave(simdata$X3, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
  # simdata$lagX4 <- ave(simdata$X4, simdata$ID, FUN = function(x) c(NA, head(x, -1)))
  
  simdata$pr_c <- 1.0
  # simdata[simdata$t == 1,]$pr_c <- as.vector(1-predict.glm(censor_model_0, 
  #                                                          newdata = data.frame(X1 = simdata[simdata$t == 1,]$lagX1,
  #                                                                               X2 = simdata[simdata$t == 1,]$lagX2,
  #                                                                               X3 = simdata[simdata$t == 1,]$lagX3,
  #                                                                               X4 = simdata[simdata$t == 1,]$lagX4),
  #                                                          type = 'response'))
  # simdata[simdata$t == 2,]$pr_c <- as.vector(1-predict.glm(censor_model_1,
  #                                                          newdata = data.frame(X1 = simdata[simdata$t == 2,]$lagX1,
  #                                                                               X2 = simdata[simdata$t == 2,]$lagX2,
  #                                                                               X3 = simdata[simdata$t == 2,]$lagX3,
  #                                                                               X4 = simdata[simdata$t == 2,]$lagX4),
  #                                                          type = 'response'))
  
  simdata$treat_weight <- ave(simdata$ps, simdata$ID, FUN = function(X) 1/cumprod(X))
  simdata$censor_weight <- ave(simdata$pr_c, simdata$ID, FUN = function(X) 1/cumprod(X))
  simdata$weight <- simdata$treat_weight*simdata$censor_weight
  
  simdata$weights <- simdata$weight
  
  
  wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, value.var = c("A", "X1", "X2","X3", "X4", "CA", "Y","C", 
                                                                         "weights"))
  
  a = min(simdata$Y,na.rm = TRUE); b = max(simdata$Y, na.rm = TRUE)
  
  wideSimdata$Y_2_scaled = (wideSimdata$Y_2-a)/(b-a)
  wideSimdata$Y_1_scaled = (wideSimdata$Y_1-a)/(b-a)
  wideSimdata$Y_0_scaled = (wideSimdata$Y_0-a)/(b-a)
  
  wideSimdata$RC_0 <- 1.0
  wideSimdata$RC_1 <- 1-wideSimdata$C_0
  wideSimdata$RC_2 <- ifelse(wideSimdata$RC_1 == 0, 0, 1-wideSimdata$C_1 )
  
  
  
  ########### t = 2 ##########
  
  Q2_2_fit <- glm(data = wideSimdata[wideSimdata$RC_2 ==1,], formula = Y_2_scaled ~ A_2 + A_1 + X1_2 + X2_2 + X3_2 +X4_2 + X1_1 + X2_1 + X3_1 +X4_1, family = 'quasibinomial')
  
  logitQ2_2 <- predict.glm(Q2_2_fit, newdata = data.frame(A_2 = c(rep(1,sample_size), rep(0,sample_size)),
                                                          A_1 = c(rep(1, sample_size), rep(0, sample_size)),
                                                          X1_2 = rep(wideSimdata$X1_2, 2), 
                                                          X2_2 = rep(wideSimdata$X2_2, 2),
                                                          X3_2 = rep(wideSimdata$X3_2, 2), 
                                                          X4_2 = rep(wideSimdata$X4_2, 2),
                                                          X1_1 = rep(wideSimdata$X1_1, 2), 
                                                          X2_1 = rep(wideSimdata$X2_1, 2),
                                                          X3_1 = rep(wideSimdata$X3_1, 2), 
                                                          X4_1 = rep(wideSimdata$X4_1, 2)),type = 'link')
  
  #------------- update Q2_2-----------------------
  regimen_ind <- c(as.numeric(wideSimdata$CA_2 ==3), as.numeric(wideSimdata$CA_2 ==0))
  regimen_ind[is.na(regimen_ind)] <- 0
  weight <- rep(wideSimdata$weights_2,2)

  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Y_2 = rep(wideSimdata$Y_2_scaled, 2), 
                            off = logitQ2_2,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(3,sample_size), rep(0,sample_size)))
  
  Q2_2star_fit <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q2_2star <- predict.glm(Q2_2star_fit, newdata = update_data, type = 'response')
  
  
  ########### t = 1 ##########
  #------------- Get Q2_1, Q1_1 ---------------------
  fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                             Q2_2star = Q2_2star, 
                             CA_1 = rep(wideSimdata$CA_1, 2),
                             A_1 = rep(wideSimdata$A_1, 2),
                             A_0 = rep(wideSimdata$A_0, 2),
                             X1_1 = rep(wideSimdata$X1_1, 2), 
                             X2_1 = rep(wideSimdata$X2_1, 2),
                             X3_1 = rep(wideSimdata$X3_1, 2), 
                             X4_1 = rep(wideSimdata$X4_1, 2))
  
  Q2_1_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_2star ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
  Q2_1_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2star ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
  
  logitQ2_1 <- as.matrix(c(predict.glm(Q2_1_fit_d1, 
                                       newdata = data.frame(CA_1 = rep(2,sample_size),
                                                            X1_1 = wideSimdata$X1_1, 
                                                            X2_1 = wideSimdata$X2_1,
                                                            X3_1 = wideSimdata$X3_1, 
                                                            X4_1 = wideSimdata$X4_1), 
                                       type = 'link'),
                           predict.glm(Q2_1_fit_d0, 
                                       newdata = data.frame(CA_1 = rep(0,sample_size),
                                                            X1_1 = wideSimdata$X1_1, 
                                                            X2_1 = wideSimdata$X2_1,
                                                            X3_1 = wideSimdata$X3_1, 
                                                            X4_1 = wideSimdata$X4_1), 
                                       type = 'link')))
  
  
  Q1_1_fit <- glm(data = wideSimdata, formula = Y_1_scaled ~ A_1 + A_0 + X1_1 + X2_1 + X3_1 + X4_1 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
  
  logitQ1_1 <- predict.glm(Q1_1_fit, newdata = data.frame(A_1 = c(rep(1,sample_size), rep(0,sample_size)),
                                                          A_0 = c(rep(1,sample_size), rep(0,sample_size)),
                                                          X1_1 = rep(wideSimdata$X1_1, 2), 
                                                          X2_1 = rep(wideSimdata$X2_1, 2),
                                                          X3_1 = rep(wideSimdata$X3_1, 2), 
                                                          X4_1 = rep(wideSimdata$X4_1, 2),
                                                          X1_0 = rep(wideSimdata$X1_0, 2), 
                                                          X2_0 = rep(wideSimdata$X2_0, 2),
                                                          X3_0 = rep(wideSimdata$X3_0, 2), 
                                                          X4_0 = rep(wideSimdata$X4_0, 2)), type = 'link')
  
  
  #------------- update Q1_1-----------------------
  regimen_ind = c(as.numeric(wideSimdata$CA_1 ==2), as.numeric(wideSimdata$CA_1 ==0))
  regimen_ind[is.na(regimen_ind)] <- 0
  weight <- rep(wideSimdata$weights_1,2)
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Y_1 = rep(wideSimdata$Y_1_scaled, 2), 
                            off = logitQ1_1,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(2,sample_size), rep(0,sample_size)))
  
  Q1_1star_fit <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q1_1star <- predict.glm(Q1_1star_fit, newdata = update_data, type = 'response')
  
  #------------- update Q2_1-----------------------
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Q2_2star = Q2_2star,
                            off = logitQ2_1,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(3,sample_size), rep(0,sample_size)))
  
  Q2_1star_fit <- glm(Q2_2star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q2_1star <- predict.glm(Q2_1star_fit, newdata = update_data, type = 'response')
  
  ########### t = 0 ##########
  #------------- Get Q4_0, Q3_0, Q2_0, Q1_0, Q0_0 ---------------------
  fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                             Q2_1star = Q2_1star, 
                             Q1_1star = Q1_1star, 
                             A_0 = rep(wideSimdata$A_0, 2),
                             X1_0 = rep(wideSimdata$X1_0, 2), 
                             X2_0 = rep(wideSimdata$X2_0, 2),
                             X3_0 = rep(wideSimdata$X3_0, 2), 
                             X4_0 = rep(wideSimdata$X4_0, 2))
  
  
  Q2_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_1star ~ A_0, family = 'quasibinomial')
  Q2_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1star ~ A_0, family = 'quasibinomial')
  
  logitQ2_0 <- as.matrix(c(predict.glm(Q2_0_fit_d1, 
                                       newdata = data.frame(A_0 = rep(1,sample_size)),
                                       type = 'link'),
                           predict.glm(Q2_0_fit_d0, 
                                       newdata = data.frame(A_0 = rep(0,sample_size)),
                                       type = 'link')))
  
 
  Q1_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q1_1star ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
  Q1_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1star ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
  
  logitQ1_0 <- as.matrix(c(predict.glm(Q1_0_fit_d1, 
                                       newdata = data.frame(A_0 = rep(1,sample_size),
                                                            X1_0 = wideSimdata$X1_0, 
                                                            X2_0 = wideSimdata$X2_0,
                                                            X3_0 = wideSimdata$X3_0, 
                                                            X4_0 = wideSimdata$X4_0),
                                       type = 'link'),
                           predict.glm(Q1_0_fit_d0, 
                                       newdata = data.frame(A_0 = rep(0,sample_size),
                                                            X1_0 = wideSimdata$X1_0, 
                                                            X2_0 = wideSimdata$X2_0,
                                                            X3_0 = wideSimdata$X3_0, 
                                                            X4_0 = wideSimdata$X4_0),
                                       type = 'link')))
  
  
  Q0_0_fit <- glm(data = wideSimdata, formula = Y_0_scaled ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
  
  logitQ0_0 <- predict.glm(Q0_0_fit, newdata = data.frame(A_0 = c(rep(1,sample_size), rep(0,sample_size)), 
                                                          X1_0 = rep(wideSimdata$X1_0, 2), 
                                                          X2_0 = rep(wideSimdata$X2_0, 2),
                                                          X3_0 = rep(wideSimdata$X3_0, 2), 
                                                          X4_0 = rep(wideSimdata$X4_0, 2)), type = 'link')
  
  
  
  
  #------------- update Q0_0-----------------------
  regimen_ind = c(as.numeric(wideSimdata$CA_0 ==1), as.numeric(wideSimdata$CA_0 ==0))
  regimen_ind[is.na(regimen_ind)] <- 0
  weight <- rep(wideSimdata$weights_0,2)
 
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Y_0 = rep(wideSimdata$Y_0_scaled, 2), 
                            off = logitQ0_0,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(1,sample_size), rep(0,sample_size)))
  
  Q0_0star_fit <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q0_0star <- predict.glm(Q0_0star_fit, newdata = update_data, type = 'response')
  
  
  #------------- update Q2_0-----------------------
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Q2_1star = Q2_1star,
                            off = logitQ2_0,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(3,sample_size), rep(0,sample_size)))
  
  Q2_0star_fit <- glm(Q2_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q2_0star <- predict.glm(Q2_0star_fit, newdata = update_data, type = 'response')
  
  #------------- update Q1_0-----------------------
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Q1_1star = Q1_1star,
                            off = logitQ1_0,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(2,sample_size), rep(0,sample_size)))
  
  Q1_0star_fit <- glm(Q1_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q1_0star <- predict.glm(Q1_0star_fit, newdata = update_data, type = 'response')
  
  
  ################# Fit MSM ##########################
  
  msm_fitting_data_transformed <- data.frame(id = rep(1:sample_size,6), 
                                             t = c(rep(0,sample_size*2), 
                                                   rep(1,sample_size*2), 
                                                   rep(2,sample_size*2)),
                                             Y = c(Q0_0star, Q1_0star, Q2_0star)*(b-a) + a,
                                             cumA = c(rep(1,sample_size), rep(0,sample_size), 
                                                      rep(2,sample_size), rep(0,sample_size), 
                                                      rep(3,sample_size), rep(0,sample_size)))
  
  msm_transformed <- glm(Y ~ cumA, data = msm_fitting_data_transformed)

  
  
  estimates <- c(msm_transformed$coefficients[1],
                 msm_transformed$coefficients[2],
                 msm_transformed$coefficients[1] + 3* msm_transformed$coefficients[2]
  )
  names(estimates) <- c('MLE LTMLE beta_0',
                        'MLE LTMLE beta_1',
                        'MLE LTMLE E(Y^1_2)')
  
  Q_t <- list(
    logitQ0_0        = logitQ0_0,
    logitQ1_0        = logitQ1_0,
    logitQ2_0        = logitQ2_0,
    logitQ1_1        = logitQ1_1,
    logitQ2_1        = logitQ2_1,
    logitQ2_2        = logitQ2_2
  )
  
  list(estimates = estimates,
       Q_t = Q_t)
}


# Modified TMLE -- the Q fits are all done first for t= 0,...,T
#
#
#
MyLtmleMSM_modified <- function(simdata, initial_Q_t = NA, refit_weights = TRUE){
  
  if (refit_weights){
    
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
    simdata$weights <- simdata$weight
    
  }
  
  wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, value.var = c("A", "X1", "X2","X3", "X4", "CA", "Y","C", 
                                                                         "weights"))
  
  a = min(simdata$Y,na.rm = TRUE); b = max(simdata$Y, na.rm = TRUE)
  
  wideSimdata$Y_2_scaled = (wideSimdata$Y_2-a)/(b-a)
  wideSimdata$Y_1_scaled = (wideSimdata$Y_1-a)/(b-a)
  wideSimdata$Y_0_scaled = (wideSimdata$Y_0-a)/(b-a)
  
  wideSimdata$RC_0 <- 1.0
  wideSimdata$RC_1 <- 1-wideSimdata$C_0
  wideSimdata$RC_2 <- ifelse(wideSimdata$RC_1 == 0, 0, 1-wideSimdata$C_1 )
  
  
  
  ########### Initial Q fits ##########
  if (all(is.na(initial_Q_t))){
    Q2_2_fit <- glm(data = wideSimdata[wideSimdata$RC_2 ==1,], formula = Y_2_scaled ~ A_2 + A_1 + X1_2 + X2_2 + X3_2 +X4_2 + X1_1 + X2_1 + X3_1 +X4_1, family = 'quasibinomial')
  } else {
    Q2_2_fit <- initial_Q_t$Q2_2_fit
  }
  logitQ2_2 <- predict.glm(Q2_2_fit, newdata = data.frame(A_2 = c(rep(1,sample_size), rep(0,sample_size)),
                                                          A_1 = c(rep(1, sample_size), rep(0, sample_size)),
                                                          X1_2 = rep(wideSimdata$X1_2, 2), 
                                                          X2_2 = rep(wideSimdata$X2_2, 2),
                                                          X3_2 = rep(wideSimdata$X3_2, 2), 
                                                          X4_2 = rep(wideSimdata$X4_2, 2),
                                                          X1_1 = rep(wideSimdata$X1_1, 2), 
                                                          X2_1 = rep(wideSimdata$X2_1, 2),
                                                          X3_1 = rep(wideSimdata$X3_1, 2), 
                                                          X4_1 = rep(wideSimdata$X4_1, 2)),type = 'link')
  
  if (all(is.na(initial_Q_t))){
    fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                               Q2_2 = plogis(logitQ2_2), 
                               CA_1 = rep(wideSimdata$CA_1, 2),
                               A_1 = rep(wideSimdata$A_1, 2),
                               A_0 = rep(wideSimdata$A_0, 2),
                               X1_1 = rep(wideSimdata$X1_1, 2), 
                               X2_1 = rep(wideSimdata$X2_1, 2),
                               X3_1 = rep(wideSimdata$X3_1, 2), 
                               X4_1 = rep(wideSimdata$X4_1, 2))
    
    Q2_1_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_2 ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
    Q2_1_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_2 ~ CA_1 + X1_1 + X2_1 + X3_1 + X4_1, family = 'quasibinomial')
    
    Q1_1_fit <- glm(data = wideSimdata, formula = Y_1_scaled ~ A_1 + A_0 + X1_1 + X2_1 + X3_1 + X4_1 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
    
  } else {
    Q2_1_fit_d1 <- initial_Q_t$Q2_1_fit_d1
    Q2_1_fit_d0 <- initial_Q_t$Q2_1_fit_d0
    Q1_1_fit <- initial_Q_t$Q1_1_fit
  }
  logitQ2_1 <- as.matrix(c(predict.glm(Q2_1_fit_d1, 
                                       newdata = data.frame(CA_1 = rep(2,sample_size),
                                                            X1_1 = wideSimdata$X1_1, 
                                                            X2_1 = wideSimdata$X2_1,
                                                            X3_1 = wideSimdata$X3_1, 
                                                            X4_1 = wideSimdata$X4_1), 
                                       type = 'link'),
                           predict.glm(Q2_1_fit_d0, 
                                       newdata = data.frame(CA_1 = rep(0,sample_size),
                                                            X1_1 = wideSimdata$X1_1, 
                                                            X2_1 = wideSimdata$X2_1,
                                                            X3_1 = wideSimdata$X3_1, 
                                                            X4_1 = wideSimdata$X4_1), 
                                       type = 'link')))
  
  
  
  logitQ1_1 <- predict.glm(Q1_1_fit, newdata = data.frame(A_1 = c(rep(1,sample_size), rep(0,sample_size)),
                                                          A_0 = c(rep(1,sample_size), rep(0,sample_size)),
                                                          X1_1 = rep(wideSimdata$X1_1, 2), 
                                                          X2_1 = rep(wideSimdata$X2_1, 2),
                                                          X3_1 = rep(wideSimdata$X3_1, 2), 
                                                          X4_1 = rep(wideSimdata$X4_1, 2),
                                                          X1_0 = rep(wideSimdata$X1_0, 2), 
                                                          X2_0 = rep(wideSimdata$X2_0, 2),
                                                          X3_0 = rep(wideSimdata$X3_0, 2), 
                                                          X4_0 = rep(wideSimdata$X4_0, 2)), type = 'link')
  
  
  if (all(is.na(initial_Q_t))){
    fitting_data <- data.frame(ID = rep(wideSimdata$ID,2), 
                               Q2_1 = plogis(logitQ2_1), 
                               Q1_1 = plogis(logitQ1_1), 
                               A_0 = rep(wideSimdata$A_0, 2),
                               X1_0 = rep(wideSimdata$X1_0, 2), 
                               X2_0 = rep(wideSimdata$X2_0, 2),
                               X3_0 = rep(wideSimdata$X3_0, 2), 
                               X4_0 = rep(wideSimdata$X4_0, 2))
    
    
    Q2_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q2_1 ~ A_0, family = 'quasibinomial')
    Q2_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q2_1 ~ A_0, family = 'quasibinomial')
    
    Q1_0_fit_d1 <- glm(data = fitting_data[1:sample_size,], formula = Q1_1 ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
    Q1_0_fit_d0 <- glm(data = fitting_data[(sample_size+1):(sample_size*2),], formula = Q1_1 ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
    
    Q0_0_fit <- glm(data = wideSimdata, formula = Y_0_scaled ~ A_0 + X1_0 + X2_0 + X3_0 + X4_0, family = 'quasibinomial')
    
  } else{
    Q2_0_fit_d1 <- initial_Q_t$Q2_0_fit_d1
    Q2_0_fit_d0 <- initial_Q_t$Q2_0_fit_d0
    Q1_0_fit_d1 <- initial_Q_t$Q1_0_fit_d1
    Q1_0_fit_d0 <- initial_Q_t$Q1_0_fit_d0
    Q0_0_fit <- initial_Q_t$Q0_0_fit
  }
  
  logitQ2_0 <- as.matrix(c(predict.glm(Q2_0_fit_d1, 
                                       newdata = data.frame(A_0 = rep(1,sample_size)),
                                       type = 'link'),
                           predict.glm(Q2_0_fit_d0, 
                                       newdata = data.frame(A_0 = rep(0,sample_size)),
                                       type = 'link')))
  logitQ1_0 <- as.matrix(c(predict.glm(Q1_0_fit_d1, 
                                       newdata = data.frame(A_0 = rep(1,sample_size),
                                                            X1_0 = wideSimdata$X1_0, 
                                                            X2_0 = wideSimdata$X2_0,
                                                            X3_0 = wideSimdata$X3_0, 
                                                            X4_0 = wideSimdata$X4_0),
                                       type = 'link'),
                           predict.glm(Q1_0_fit_d0, 
                                       newdata = data.frame(A_0 = rep(0,sample_size),
                                                            X1_0 = wideSimdata$X1_0, 
                                                            X2_0 = wideSimdata$X2_0,
                                                            X3_0 = wideSimdata$X3_0, 
                                                            X4_0 = wideSimdata$X4_0),
                                       type = 'link')))
  
  logitQ0_0 <- predict.glm(Q0_0_fit, newdata = data.frame(A_0 = c(rep(1,sample_size), rep(0,sample_size)), 
                                                          X1_0 = rep(wideSimdata$X1_0, 2), 
                                                          X2_0 = rep(wideSimdata$X2_0, 2),
                                                          X3_0 = rep(wideSimdata$X3_0, 2), 
                                                          X4_0 = rep(wideSimdata$X4_0, 2)), type = 'link')
  
  # Targrting step
  #------------- update Q2_2-----------------------
  regimen_ind <- c(as.numeric(wideSimdata$CA_2 ==3), as.numeric(wideSimdata$CA_2 ==0))
  regimen_ind[is.na(regimen_ind)] <- 0
  weight <- rep(wideSimdata$weights_2,2)
 
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Y_2 = rep(wideSimdata$Y_2_scaled, 2), 
                            off = logitQ2_2,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(3,sample_size), rep(0,sample_size)))
  
  Q2_2star_fit <- glm(Y_2 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q2_2star <- predict.glm(Q2_2star_fit, newdata = update_data, type = 'response')
  
  
  ########### t = 1 ##########
  #------------- Get Q2_1, Q1_1 ---------------------
  
  #------------- update Q1_1-----------------------
  regimen_ind = c(as.numeric(wideSimdata$CA_1 ==2), as.numeric(wideSimdata$CA_1 ==0))
  regimen_ind[is.na(regimen_ind)] <- 0
  weight <- rep(wideSimdata$weights_1,2)
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Y_1 = rep(wideSimdata$Y_1_scaled, 2), 
                            off = logitQ1_1,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(2,sample_size), rep(0,sample_size)))
  
  Q1_1star_fit <- glm(Y_1 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q1_1star <- predict.glm(Q1_1star_fit, newdata = update_data, type = 'response')
  
  
  #------------- update Q2_1-----------------------
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Q2_2star = Q2_2star,
                            off = logitQ2_1,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(3,sample_size), rep(0,sample_size)))
  
  Q2_1star_fit <- glm(Q2_2star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q2_1star <- predict.glm(Q2_1star_fit, newdata = update_data, type = 'response')

  
  ########### t = 0 ##########
  #------------- Get Q4_0, Q3_0, Q2_0, Q1_0, Q0_0 ---------------------
  
  
  
  
  #------------- update Q0_0-----------------------
  regimen_ind = c(as.numeric(wideSimdata$CA_0 ==1), as.numeric(wideSimdata$CA_0 ==0))
  regimen_ind[is.na(regimen_ind)] <- 0
  weight <- rep(wideSimdata$weights_0,2)
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Y_0 = rep(wideSimdata$Y_0_scaled, 2), 
                            off = logitQ0_0,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(1,sample_size), rep(0,sample_size)))
  
  Q0_0star_fit <- glm(Y_0 ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q0_0star <- predict.glm(Q0_0star_fit, newdata = update_data, type = 'response')
  
  
  #------------- update Q2_0-----------------------
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Q2_1star = Q2_1star,
                            off = logitQ2_0,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(3,sample_size), rep(0,sample_size)))
  
  Q2_0star_fit <- glm(Q2_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q2_0star <- predict.glm(Q2_0star_fit, newdata = update_data, type = 'response')
  
  #------------- update Q1_0-----------------------
  
  update_data <- data.frame(id = rep(wideSimdata$ID,2),
                            Q1_1star = Q1_1star,
                            off = logitQ1_0,
                            intercept = rep(1,sample_size*2), 
                            cumA = c(rep(2,sample_size), rep(0,sample_size)))
  
  Q1_0star_fit <- glm(Q1_1star ~ intercept + cumA + offset(off) - 1, family = 'quasibinomial', 
                      data = update_data[regimen_ind == 1,],
                      weights = as.vector(weight[regimen_ind == 1]))
  
  Q1_0star <- predict.glm(Q1_0star_fit, newdata = update_data, type = 'response')
  
  ################# Fit MSM ##########################
  
  msm_fitting_data_transformed <- data.frame(id = rep(1:sample_size,6), 
                                             t = c(rep(0,sample_size*2), 
                                                   rep(1,sample_size*2), 
                                                   rep(2,sample_size*2)),
                                             Y = c(Q0_0star, Q1_0star, Q2_0star)*(b-a) + a, 
                                             cumA = c(rep(1,sample_size), rep(0,sample_size), 
                                                      rep(2,sample_size), rep(0,sample_size), 
                                                      rep(3,sample_size), rep(0,sample_size)))
  
  msm_transformed <- glm(Y ~ cumA, data = msm_fitting_data_transformed)

  
  estimates <- c(msm_transformed$coefficients[1],
                 msm_transformed$coefficients[2],
                 msm_transformed$coefficients[1] + 3* msm_transformed$coefficients[2]
  )
  names(estimates) <- c('MLE LTMLE beta_0',
                        'MLE LTMLE beta_1',
                        'MLE LTMLE E(Y^1_2)')
  
  Q_fits <- list(
    Q2_2_fit      = Q2_2_fit,
    Q1_1_fit      = Q1_1_fit,       # if you really want to include the duplicate
    Q2_1_fit_d1   = Q2_1_fit_d1,
    Q2_1_fit_d0   = Q2_1_fit_d0,  # duplicate
    Q2_0_fit_d1   = Q2_0_fit_d1,
    Q2_0_fit_d0   = Q2_0_fit_d0,
    Q1_0_fit_d1   = Q1_0_fit_d1,
    Q1_0_fit_d0   = Q1_0_fit_d0,
    Q0_0_fit      = Q0_0_fit
  )
  
  list(estimates = estimates,
       Q_t = Q_fits,
       simdata = simdata)
}
