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
load('hers.Rdata')

############# DATA PREPARATION ##################
HERS$Y <- as.factor(HERS$Y)
HERS$t <- HERS$visit - 8
HERS$SITE1 <- as.factor(HERS$SITE1)
HERS$SITE2 <- as.factor(HERS$SITE2)
HERS$SITE3 <- as.factor(HERS$SITE3)
HERS$WHITE <- as.factor(HERS$WHITE)
HERS$OTHER <- as.factor(HERS$OTHER)

HERS$CD4 <- (sqrt(as.numeric(HERS$CD4)) - mean(sqrt(HERS$CD4)))/sd(sqrt(HERS$CD4))
HERS$CD4_1 <- (sqrt(as.numeric(HERS$CD4_1)) - mean(sqrt(HERS$CD4_1)))/sd(sqrt(HERS$CD4_1))
HERS$CD4_2 <- (sqrt(as.numeric(HERS$CD4_2)) - mean(sqrt(HERS$CD4_2)))/sd(sqrt(HERS$CD4_2))

HERS$viral <- (log10(HERS$viral) - mean(log10(HERS$viral)))/sd(log10(HERS$viral))
HERS$viral_1 <- (log10(HERS$viral_1) - mean(log10(HERS$viral_1)))/sd(log10(HERS$viral_1))
HERS$viral_2 <- (log10(HERS$viral_2) - mean(log10(HERS$viral_2)))/sd(log10(HERS$viral_2))

HERS <- HERS %>% 
  dplyr::arrange(id,t) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(CAp = ifelse(t == 0, haart_1 + haart_2, cumsum(haart_1))) %>% 
  dplyr::ungroup()

HERS$A <- HERS$haart
HERS$Ap <- HERS$haart_1
HERS$App <- HERS$haart_2
HERS[,'ID'] <- HERS$id
HERS <- HERS %>% 
  dplyr::mutate(SITE = as.factor(ifelse(SITE1 == 1, 1, ifelse(SITE2 == 1,2,3)))) %>% 
  dplyr::select(ID, t, A, Ap, App,CAp, CD4, CD4_1,CD4_2,
                viral,viral_1,viral_2,HIVsym,HIVsym_1,HIVsym_2,
                SITE, WHITE, OTHER, Y, C) %>% 
  dplyr::arrange(ID, t) %>% 
  dplyr:: group_by(ID) %>% 
  dplyr::mutate(cumA = cumsum(A),
                A_0 = first(A),
                protocol = ifelse(A_0 == 1, t+1, 0)) %>% 
  dplyr::ungroup()
HERS$eligible <- as.numeric(HERS$CAp == 0)

 
############## SEQUENTIAL G COMPUTATION ###################
library(data.table)
wideHERS <- data.table::dcast(setDT(HERS), ID + WHITE + OTHER ~ t, value.var = c("viral", "HIVsym", "SITE", "CD4", "A", "C"))

wideHERS$y4pred <- wideHERS$CD4_4

q12 <- glm(y4pred ~ WHITE + OTHER + A_0 + A_1 + A_2 + A_3 + A_4 + HIVsym_2 + HIVsym_3 + HIVsym_4 + CD4_2 + CD4_3,
           data = wideHERS[wideHERS$C_3 == 0,], family = 'gaussian')

wideHERS$y4pred <- predict.glm(q12,type = "response", newdata = wideHERS %>% mutate(A_1 = A_0, A_2 = A_0, A_3 = A_0, A_4 = A_0))

q11 <- glm(y4pred ~ WHITE + OTHER + A_0 + A_1 + A_2 + A_3 + HIVsym_1 + HIVsym_2 + HIVsym_3 + CD4_1 + CD4_2,
           data = wideHERS[wideHERS$C_2 == 0,], family = 'gaussian')

wideHERS$y4pred <- predict.glm(q11,type = "response", newdata = wideHERS %>% mutate(A_1 = A_0, A_2 = A_0, A_3 = A_0))

q10 <- glm(y4pred ~ WHITE + OTHER + A_0 + A_1 + A_2 + HIVsym_0 + HIVsym_1 + HIVsym_2 + CD4_0 + CD4_1,
           data = wideHERS[wideHERS$C_1 == 0,], family = 'gaussian')

wideHERS$y4pred <- predict.glm(q10,type = "response", newdata = wideHERS %>% mutate(A_1 = A_0, A_2 = A_0))

q9 <- glm(y4pred ~ WHITE + OTHER + A_0 + A_1 + HIVsym_0 + HIVsym_1  + CD4_0,
           data = wideHERS[wideHERS$C_0 == 0,], family = 'gaussian')

wideHERS$y4pred <- predict.glm(q9,type = "response", newdata = wideHERS %>% mutate(A_1 = A_0))

q8 <- glm(y4pred ~ WHITE + OTHER + A_0 + HIVsym_0,
          data = wideHERS, family = 'gaussian')

wideHERS$y4pred <- predict.glm(q8,type = "response", newdata = wideHERS)

wideHERS$Y0 <- predict.glm(q8, type = "response", newdata = wideHERS)
wideHERS$Y1 <- predict.glm(q9, type = "response", newdata = wideHERS)
wideHERS$Y2 <- predict.glm(q10, type = "response", newdata = wideHERS)
wideHERS$Y3 <- predict.glm(q11, type = "response", newdata = wideHERS)
wideHERS$Y4 <- predict.glm(q12, type = "response", newdata = wideHERS)

HERS_imputed <- pivot_longer(wideHERS %>% dplyr::select(ID, Y0, Y1, Y2, Y3, Y4),
                             cols = Y0:Y4, names_to = 't', names_prefix = "Y", values_to = "Yhat") %>% 
  mutate(t = as.numeric(t))

######################### CALIBRATION #################################################
HERS <- HERS %>% 
  merge(HERS_imputed, by = c("ID", "t")) %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1)) %>% 
  dplyr::mutate(CRA = cumsum(RA),
                RA = ifelse(CRA == t+1,1,0),
                RC = ifelse(lag(C) == 0,1,0),
                A1 = A_0,
                A0 = 1-A_0,
                tA1 = t*A_0,
                tA0 = t*(1-A_0),
                A1HIVsym = A_0*HIVsym,
                A0HIVsym = (1-A_0)*HIVsym,
                tA1HIVsym =t* A_0*HIVsym,
                tA0HIVsym = t*(1-A_0)*HIVsym,
                A1HIVsym_1 = A_0*HIVsym_1,
                A0HIVsym_1 = (1-A_0)*HIVsym_1,
                tA1HIVsym_1 =t* A_0*HIVsym_1,
                tA0HIVsym_1 = t*(1-A_0)*HIVsym_1,
                A1HIVsym_2 = A_0*HIVsym_2,
                A0HIVsym_2 = (1-A_0)*HIVsym_2,
                tA1HIVsym_2 =t* A_0*HIVsym_2,
                tA0HIVsym_2 = t*(1-A_0)*HIVsym_2,
                A1CD4_1 = A_0*CD4_1,
                A0CD4_1 = (1-A_0)*CD4_1,
                tA1CD4_1 =t* A_0*CD4_1,
                tA0CD4_1 = t*(1-A_0)*CD4_1,
                A1CD4_2 = A_0*CD4_2,
                A0CD4_2 = (1-A_0)*CD4_2,
                tA1CD4_2 =t* A_0*CD4_2,
                tA0CD4_2 = t*(1-A_0)*CD4_2,
                A1viral = A_0*viral,
                A0viral = (1-A_0)*viral,
                tA1viral =t*A_0*viral,
                tA0viral = t*(1-A_0)*viral,
                A1viral_1 = A_0*viral_1,
                A0viral_1 = (1-A_0)*viral_1,
                tA1viral_1 =t* A_0*viral_1,
                tA0viral_1 = t*(1-A_0)*viral_1,
                A1viral_2 = A_0*viral_2,
                A0viral_2 = (1-A_0)*viral_2,
                tA1viral_2 =t*A_0*viral_2,
                tA0viral_2 = t*(1-A_0)*viral_2,
                A1SITE = A_0*as.numeric(SITE),
                A0SITE = (1-A_0)*as.numeric(SITE),
                tA1SITE =t*A_0*as.numeric(SITE),
                tA0SITE = t*(1-A_0)*as.numeric(SITE),
                A1WHITE = A_0*as.numeric(WHITE),
                A0WHITE = (1-A_0)*as.numeric(WHITE),
                tA1WHITE =t* A_0*as.numeric(WHITE),
                tA0WHITE = t*(1-A_0)*as.numeric(WHITE),
                A1OTHER = A_0*as.numeric(OTHER),
                A0OTHER = (1-A_0)*as.numeric(OTHER),
                tA1OTHER =t*A_0*as.numeric(OTHER),
                tA0OTHER = t*(1-A_0)*as.numeric(OTHER),
                A1Yhat = A_0*as.numeric(Yhat),
                A0Yhat = (1-A_0)*as.numeric(Yhat),
                tA1Yhat =t*A_0*as.numeric(Yhat),
                tA0Yhat = t*(1-A_0)*as.numeric(Yhat),
                sub = ID,
                tall = t,
                One = 1.0) %>% 
  dplyr::arrange(ID, t) 

#weight_training_data <- HERS %>% 
 # group_by(ID) %>% 
  #filter(cumsum(1-RA) <= 1, t!= 0)

weight_model <- glm(data = HERS, 
                    formula = A ~ Ap + CD4_1 + CD4_2 + viral_1 + viral_2 + HIVsym_1 + HIVsym_2 + SITE + WHITE + OTHER, family = 'binomial')
summary(weight_model)
HERS$p_1 <- 1.0
HERS[HERS$t != 0,]$p_1 <- predict.glm(weight_model, HERS[HERS$t != 0,], type = 'response')
HERS <- HERS %>% 
  arrange(ID, t) %>% 
  group_by(ID) %>% 
  dplyr::mutate(
    wt = ifelse( t == 0, 1.0,ifelse(A == 1, 1/p_1, 1/(1-p_1))),
    wtprod = cumprod(wt),
    weights = wtprod)

################### Calibration aggregated#######################
simdatafinal1 <- calibration(simdatafinal = HERS, 
                             var = c('A1','A1CD4_1', 'A1CD4_2', 'A1viral_1', 'A1viral_2', 'A1HIVsym_1', 'A1HIVsym_2', 'A1SITE', 'A1WHITE', 'A1OTHER','A1Yhat',
                                     'A0','A0CD4_1', 'A0CD4_2', 'A0viral_1', 'A0viral_2', 'A0HIVsym_1', 'A0HIVsym_2', 'A0SITE', 'A0WHITE', 'A0OTHER','A0Yhat'))


result$objectiveIPW <- simdatafinal1$objective.IPW
result$objectiveCali <- simdatafinal1$objective.Cali

################### Calibration by time #######################
simdatafinal2 <- calibration_by_time(simdatafinal = HERS, 
                                     var = c('A1','A1CD4_1', 'A1CD4_2', 'A1viral_1', 'A1viral_2', 'A1HIVsym_1', 'A1HIVsym_2', 'A1SITE', 'A1WHITE', 'A1OTHER','A1Yhat',
                                             'A0','A0CD4_1', 'A0CD4_2', 'A0viral_1', 'A0viral_2', 'A0HIVsym_1', 'A0HIVsym_2', 'A0SITE', 'A0WHITE', 'A0OTHER','A0Yhat'))


result$objectiveIPWseq <- simdatafinal2$objective.IPW
result$objectiveCaliseq <- simdatafinal2$objective.Cali

################## Calibration aggregated with time interaction ###########################
simdatafinal3 <- calibration(simdatafinal = HERS, 
                             var =c('tA1','tA1CD4_1', 'tA1CD4_2', 'tA1viral_1', 'tA1viral_2', 'tA1HIVsym_1', 'tA1HIVsym_2', 'tA1SITE', 'tA1WHITE', 'tA1OTHER','A1Yhat',
                                    'tA0','tA0CD4_1', 'tA0CD4_2', 'tA0viral_1', 'tA0viral_2', 'tA0HIVsym_1', 'tA0HIVsym_2', 'tA0SITE', 'tA0WHITE', 'tA0OTHER','A0Yhat'))
                             


result$objectiveIPWCaliT <- simdatafinal3$objective.IPW
result$objectiveCaliT <- simdatafinal3$objective.Cali


HERS$weight <- simdatafinal1$data$weights
HERS$Cweights <- simdatafinal1$data$Cweights
HERS$Cweights_sequential <- simdatafinal2$data$Cweights
HERS$CweightsT <- simdatafinal3$data$Cweights


summary(HERS[HERS$RA == 1 & HERS$A0 == 0,] %>% dplyr::select(weights, Cweights, Cweights_sequential, CweightsT))

########################## COUNTERFACTUAL OUTCOME MODELLING ####################
PP_naive <- glm(data = HERS[HERS$RA == 1 & HERS$t == 4,],
                formula = CD4 ~ cumA + WHITE + OTHER + SITE1 + STE2 + SITE3,
                weights = NULL, family = 'gaussian')
summary(PP_naive)

PP_IPW <- glm(data = HERS[HERS$RA == 1 & HERS$t == 4,],
              formula = CD4 ~ cumA + WHITE + OTHER + SITE1 + STE2 + SITE3,
              weights = weight, family = 'gaussian')
summary(PP_IPW)
PP_calibrated <- glm(data = HERS[HERS$RA == 1 & HERS$t == 4,],
                     formula = CD4 ~ cumA + WHITE + OTHER + SITE1 + STE2 + SITE3,
                     weights = Cweights, family = 'gaussian')
summary(PP_calibrated)
PP_calibrated_sequential <- glm(data = HERS[HERS$RA == 1 & HERS$t == 4,],
                                formula = CD4 ~ cumA + WHITE + OTHER + SITE1 + STE2 + SITE3,
                                weights = Cweights_sequential, family = 'gaussian')
summary(PP_calibrated_sequential)
PP_calibrated_T <- glm(data = HERS[HERS$RA == 1 & HERS$t == 4,],
                       formula = CD4 ~ cumA,
                       weights = CweightsT, family = 'gaussian')
summary(PP_calibrated_T)

ATE_naive <- predict.glm(PP_naive,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 5), type = 'response')[1] -  predict.glm(PP_naive,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 0), type = 'response')[1]
ATE_IPW <- predict.glm(PP_IPW,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 5), type = 'response')[1] -  predict.glm(PP_IPW,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 0), type = 'response')[1]
ATE_calibrated <- predict.glm(PP_calibrated,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 5), type = 'response')[1] -  predict.glm(PP_calibrated,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 0), type = 'response')[1]
ATE_calibrated_sequential <- predict.glm(PP_calibrated_sequential,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 5), type = 'response')[1] -  predict.glm(PP_calibrated_sequential,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 0), type = 'response')[1]
ATE_calibrated_T <- predict.glm(PP_calibrated_T,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 5), type = 'response')[1] -  predict.glm(PP_calibrated_T,newdata = HERS[HERS$RA == 1 & HERS$t == 4,] %>% dplyr::mutate(cumA = 0), type = 'response')[1]



