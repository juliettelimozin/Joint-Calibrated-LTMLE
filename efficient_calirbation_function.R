eff_calibration_by_time_from_baseline <- function(simdatafinal,
                                              calibrate_treat = TRUE,
                                              var = c("X1", "X2"),
                                              subset = FALSE,
                                              subset_t = NULL,
                                              subset_var = NULL,
                                              censor = FALSE,
                                              c_var = NULL,
                                              subset_c_var = NULL,
                                              weights_var = "weights") {
  simdatafinal$weights <- simdatafinal[, weights_var]
  T <- max(simdatafinal$tall)
  
  if (calibrate_treat) {
    data1 <- simdatafinal[simdatafinal$tall == 0, ]
    data1 <- cbind(data1$S1, data1$S1 * data1[, var],data1$S0, data1$S0 * data1[, var]) # RC_0*RA_0*X_0 = RA_0*X_0
    Tdata1 <- t(data1)
    data0 <- cbind(1, simdatafinal[simdatafinal$tall == 0, var],1, simdatafinal[simdatafinal$tall == 0, var]) # RC_0*RA_-1*W_-1*X_0 = X_0
    RHS <- colSums(data0)
    restrictions_weight1 <- function(w) {
      ## Objective function
      lp3 <- colSums(Tdata1 * w)
      we <- simdatafinal$weights[simdatafinal$tall == 0] * exp(lp3)
      m3 <- (colSums(data1 * we) - RHS) / nrow(data1) ## Restrictions (6)
      c(m3)
    }
    
    ## Hessian matrix
    
    N1MATAR <- Tdata1
    N2MATAR <- data1
    
    DgAR1 <- function(w) {
      lp3 <- colSums(Tdata1 * w)
      we <- simdatafinal$weights[simdatafinal$tall == 0] * exp(lp3)
      MAT <- N2MATAR * we
      MAT1 <- diag(length(w))
      for (k in 1:length(w)) {
        MAT1[, k] <- colMeans(N1MATAR[k, ] * MAT)
      }
      MAT1
    }
    
    wei1optAR <- nleqslv::nleqslv(rep(0, dim(data1)[2]), restrictions_weight1, DgAR1,
                                  method = "Broyden",
                                  control = list(maxit = 10000, ftol = 10^(-16), xtol = 10^(-16)), jacobian = T
    )
    
    weconsAR <- function(w) {
      lp3 <- colSums(Tdata1 * w)
      we <- simdatafinal$weights[simdatafinal$tall == 0] * exp(lp3)
      return(we)
    }
    
    CALW1 <- weconsAR(wei1optAR$x)
    
    simdatafinal$Cweights <- simdatafinal$weights
    simdatafinal$Cweights[simdatafinal$tall == 0] <- CALW1
  } else {
    simdatafinal$Cweights <- simdatafinal$weights
  }
  
  if (censor == TRUE) {
    c_var <- paste0("lag", c_var) # Create H_k-1
    if (subset == TRUE) {
      subset_c_var <- paste0("lag", subset_c_var) # Create H_k-1
    }
    simdatafinal$lagS1 <- ave(simdatafinal$S1, simdatafinal$sub, FUN = function(x) c(0.0, head(x, -1))) # Create RA_k-1
    simdatafinal$lagS0 <- ave(simdatafinal$S0, simdatafinal$sub, FUN = function(x) c(0.0, head(x, -1))) # Create RA_k-1
  }
  for (vars in c(var, c_var)) {
    simdatafinal[[paste0(vars, "_lead")]] <- ave(simdatafinal[[vars]], simdatafinal$sub, FUN = function(x) c(tail(x, -1), 0.0)) # Create H_k, X_k+1
  }
  
  var_lead <- paste0(var, "_lead")
  if (subset == TRUE) {
    subset_var_lead <- paste0(subset_var, "_lead")
  }
  if (censor == TRUE) {
    c_var_lead <- paste0(c_var, "_lead")
    if (subset == TRUE) {
      subset_c_var_lead <- paste0(subset_c_var, "_lead")
    }
  }
  for (k in 1:T) {
    data1RA <- simdatafinal[simdatafinal$tall == k, ]
    
    if (isTRUE(k >= subset_t)) {
      data1RA <- cbind(data1RA$S1, data1RA$S1 * data1RA[, subset_var],data1RA$S0, data1RA$S0 * data1RA[, subset_var])
    } else {
      data1RA <- cbind(data1RA$S1, data1RA$S1 * data1RA[, var],data1RA$S0, data1RA$S0 * data1RA[, var]) # RC_k*RA_k*X_k
    }
    
    if (censor == TRUE) {
      if (isTRUE(k >= subset_t)) {
        data1RC <- cbind(simdatafinal[simdatafinal$tall == k, ]$lagS1, simdatafinal[simdatafinal$tall == k, ]$lagS1 * simdatafinal[simdatafinal$tall == k, subset_c_var],
                         simdatafinal[simdatafinal$tall == k, ]$lagS0, simdatafinal[simdatafinal$tall == k, ]$lagS0 * simdatafinal[simdatafinal$tall == k, subset_c_var])
      } else {
        data1RC <- cbind(simdatafinal[simdatafinal$tall == k, ]$lagS1, simdatafinal[simdatafinal$tall == k, ]$lagS1 * simdatafinal[simdatafinal$tall == k, c_var],
                         simdatafinal[simdatafinal$tall == k, ]$lagS0, simdatafinal[simdatafinal$tall == k, ]$lagS0 * simdatafinal[simdatafinal$tall == k, c_var])
        # RA_k-1*RC_k*H_k-1
      }
      
      if (calibrate_treat) {
        Tdata1 <- t(cbind(data1RA, data1RC))
      } else {
        Tdata1 <- t(data1RC)
      }
    } else {
      Tdata1 <- t(data1RA)
    }
    data1 <- t(Tdata1)
    
    
    data0RA <- simdatafinal[simdatafinal$tall == k - 1, ]
    if (isTRUE(k >= subset_t)) {
      data0RA <- cbind((1 - data0RA$C) * data0RA$S1 * data0RA$Cweights, (1 - data0RA$C) * data0RA$S1 * data0RA$Cweights * data0RA[, subset_var_lead],
                       (1 - data0RA$C) * data0RA$S0 * data0RA$Cweights, (1 - data0RA$C) * data0RA$S0 * data0RA$Cweights * data0RA[, subset_var_lead])
    } else {
      data0RA <- cbind((1 - data0RA$C) * data0RA$S1 * data0RA$Cweights, (1 - data0RA$C) * data0RA$S1 * data0RA$Cweights * data0RA[, var_lead],
                       (1 - data0RA$C) * data0RA$S0 * data0RA$Cweights, (1 - data0RA$C) * data0RA$S0 * data0RA$Cweights * data0RA[, var_lead])
      # RC_k*RA_k-1 * W_k-1 *X_k = (1-C_k-1)*RA_k-1 * W_k-1 *X_k
    }
    
    RHS_RA <- colSums(data0RA)
    if (censor == TRUE) {
      data0RC <- simdatafinal[simdatafinal$tall == k - 1, ]
      
      if (isTRUE(k >= subset_t)) {
        data0RC <- cbind(data0RC$S1 * data0RC$Cweights, data0RC$S1 * data0RC$Cweights * data0RC[, subset_c_var_lead],
                         data0RC$S0 * data0RC$Cweights, data0RC$S0 * data0RC$Cweights * data0RC[, subset_c_var_lead]) # RA_k-1*RC_k-1*W_k-1*H_k-1
      } else {
        data0RC <- cbind(data0RC$S1 * data0RC$Cweights, data0RC$S1 * data0RC$Cweights * data0RC[, c_var_lead],
                         data0RC$S0 * data0RC$Cweights, data0RC$S0 * data0RC$Cweights * data0RC[, c_var_lead]) # RA_k-1*RC_k-1*W_k-1*H_k-1
      }
      RHS_RC <- colSums(data0RC)
      if (calibrate_treat) {
        RHS <- c(RHS_RA, RHS_RC)
      } else {
        RHS <- RHS_RC
      }
    } else {
      RHS <- RHS_RA
    }
    
    
    restrictions_weightk <- function(w) {
      ## Objective function
      lp3 <- colSums(Tdata1 * w)
      we <- simdatafinal$weights[simdatafinal$tall == k] * exp(lp3)
      m3 <- (colSums(data1 * we) - RHS) / nrow(data1) ## Restrictions (6)
      c(m3)
    }
    
    ## Hessian matrix
    
    N1MATAR <- Tdata1
    N2MATAR <- data1
    
    DgAR1 <- function(w) {
      lp3 <- colSums(Tdata1 * w)
      we <- simdatafinal$weights[simdatafinal$tall == k] * exp(lp3)
      MAT <- N2MATAR * we
      MAT1 <- diag(length(w))
      for (j in 1:length(w)) {
        MAT1[, j] <- colMeans(N1MATAR[j, ] * MAT)
      }
      MAT1
    }
    
    weikoptAR <- nleqslv::nleqslv(rep(0, dim(data1)[2]), restrictions_weightk, DgAR1,
                                  method = "Broyden",
                                  control = list(maxit = 10000, ftol = 10^(-16), xtol = 10^(-16)), jacobian = T
    )
    
    weconsAR <- function(w) {
      lp3 <- colSums(Tdata1 * w)
      we <- simdatafinal$weights[simdatafinal$tall == k] * exp(lp3)
      return(we)
    }
    
    CALW1 <- weconsAR(weikoptAR$x)
    
    simdatafinal$Cweights[simdatafinal$tall == k] <- CALW1
  }
  list(
    data = simdatafinal,
    objective.IPW = NA,
    objective.Cali = NA
  )
}


#!/usr/bin R
###############################
# This script runs simulations with Q functions on logit form
###############################
# setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
load("/home/juliette/hers.Rdata")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/calibration_func_trials.R")
source("/home/juliette/Calibrated-weights-sequential-trial-emulation/hers_ltmle_utils.R")
library(dplyr)
library(MASS)
library(foreach)
library(doParallel)
library(doRNG)
library(nleqslv)
library(data.table)
library(matrixStats)
library(xtable)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(cobalt)
set.seed(160625)

# Data preparation
simdata <- HERS

simdata$C <- simdata$C + simdata$Y
simdata$t <- simdata$visit - 8
simdata$SITE1 <- as.integer(simdata$SITE1)
simdata$SITE2 <- as.integer(simdata$SITE2)
simdata$SITE3 <- as.integer(simdata$SITE3)
simdata$WHITE <- as.integer(simdata$WHITE)
simdata$OTHER <- as.integer(simdata$OTHER)

simdata <- simdata %>%
  dplyr::arrange(id, t) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(CAp = cumsum(as.numeric(haart_1 == 1 | haart_2 == 1))) %>%
  dplyr::ungroup()

simdata$A <- simdata$haart
simdata$Ap <- simdata$haart_1
simdata$App <- simdata$haart_2
simdata[, "ID"] <- simdata$id

simdata$sqrtCD4 <- (sqrt(as.numeric(simdata$CD4)) - mean(sqrt(simdata$CD4))) / sd(sqrt(simdata$CD4))
simdata$sqrtCD4_1 <- (sqrt(as.numeric(simdata$CD4_1)) - mean(sqrt(simdata$CD4_1))) / sd(sqrt(simdata$CD4_1))
simdata$sqrtCD4_2 <- (sqrt(as.numeric(simdata$CD4_2)) - mean(sqrt(simdata$CD4_2))) / sd(sqrt(simdata$CD4_2))

simdata <- simdata %>%
  dplyr::select(
    ID, t, A, Ap, App, CAp, CD4, CD4_1, CD4_2, sqrtCD4, sqrtCD4_1, sqrtCD4_2,
    viral, viral_1, viral_2, HIVsym, HIVsym_1, HIVsym_2,
    SITE1, SITE2, SITE3, WHITE, OTHER, Y, C
  )
simdata$eligible <- as.numeric(simdata$CAp == 0 & simdata$t == 0)
simdata$CA <- ave(simdata$A, simdata$ID, FUN = cumsum)
simdata <- as.data.frame(simdata)
simdata$eligible <- ave(simdata$eligible, simdata$ID, FUN = cumsum)


simdata <- simdata[simdata$eligible == 1, ]

simdata <- simdata %>%
  group_by(ID) %>%
  mutate(D = factor(
    case_when(
      CD4_1[t == 0][1] < 200 ~ 0,
      CD4_1[t == 0][1] >= 200 & CD4_1[t == 0][1] <= 500 ~ 1,
      CD4_1[t == 0][1] > 500 ~ 2
    ),
    levels = c(0, 1, 2)
  )) %>%
  ungroup()

simdata$viral <- (log10(simdata$viral) - mean(log10(simdata$viral))) / sd(log10(simdata$viral))
simdata$viral_1 <- (log10(simdata$viral_1) - mean(log10(simdata$viral_1))) / sd(log10(simdata$viral_1))
simdata$viral_2 <- (log10(simdata$viral_2) - mean(log10(simdata$viral_2))) / sd(log10(simdata$viral_2))

simdata <- as.data.frame(simdata)

calculate_weights <- initial_weight_calculation(simdata)
simdata <- calculate_weights$data
treat_models_n <- calculate_weights$treat_models_n
censor_models_n <- calculate_weights$censor_models_n
simdata <- weight_calibration(simdata)

simdata_old_cali <- simdata


simdata$S1 <- 1
simdata[simdata$t == 0 & !(simdata$CA == 1), ]$S1 <- 0
simdata[simdata$t == 1 & !(simdata$CA == 2), ]$S1 <- 0
simdata[simdata$t == 2 & !(simdata$CA == 3), ]$S1 <- 0
simdata[simdata$t == 3 & !(simdata$CA == 4), ]$S1 <- 0
simdata[simdata$t == 4 & !(simdata$CA == 5), ]$S1 <- 0

simdata$S0 <- 1
simdata[simdata$t == 0 & !(simdata$CA == 0), ]$S0 <- 0
simdata[simdata$t == 1 & !(simdata$CA == 0), ]$S0 <- 0
simdata[simdata$t == 2 & !(simdata$CA == 0), ]$S0 <- 0
simdata[simdata$t == 3 & !(simdata$CA == 0), ]$S0 <- 0
simdata[simdata$t == 4 & !(simdata$CA == 0), ]$S0 <- 0



calibrate_always_treated <- eff_calibration_by_time_from_baseline(simdata,
                                                              var = c("sqrtCD4_1", "sqrtCD4_2", "viral_1", "viral_2", "HIVsym_1", "SITE1", "SITE2", "SITE3", "WHITE", "OTHER"),
                                                              censor = TRUE,
                                                              c_var = c("sqrtCD4_1", "sqrtCD4", "viral", "viral_1", "HIVsym", "HIVsym_1", "SITE1", "SITE2", "SITE3", "WHITE", "OTHER"),
                                                              weights_var = "weight"
)
simdata$Cweights <- calibrate_always_treated$data$Cweights
all(simdata$Cweights == simdata_old_cali$Cweights)

simdata_old_cali$S1 <- simdata$S1
simdata_old_cali$S0 <- simdata$S0
mean(simdata[simdata$S1 == 1,]$Cweights - simdata_old_cali[simdata$S1 == 1,]$Cweights)
mean(simdata[simdata$S0 == 1,]$Cweights - simdata_old_cali[simdata$S0 == 1,]$Cweights)


