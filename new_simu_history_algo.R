library(dplyr)
library(tidyr)
#setwd("~/rds/hpc-work/Calibrated-weights-sequential-trial-emulation")
library(TrialEmulation)
library(MASS)
library(sandwich)
library(foreach)
library(doParallel)
library(doRNG)
library(nleqslv)
source('calibration_func_trials.R')
set.seed(NULL)
library(SuperLearner)
## simulate data for testing TrialEmulation package, using the algorithm in Young and Tchetgen Tchetgen (2014) 


DATA_GEN_cont_outcome_treatment_switch_history<-function(ns, nv, conf = 0.5, treat_prev = 0,
                                                                  outcome_prev = -3.8, all_treat = FALSE,
                                                                  all_control = FALSE, censor = TRUE){   
  # ns= number of subjects, nv=no of visits including baseline visit
  
  
  nvisit<-nv+1
  
  X1<-rep(0,nvisit*ns)          ## place holders for time-varying covariates
  Z1<-rnorm(nvisit*ns,0,0.5)
  X2<-rep(0,nvisit*ns)    
  Z2<-rnorm(nvisit*ns,0,0.5)
  
  X1p <- rep(0,nvisit*ns)
  X2p <- rep(0,nvisit*ns)
  
  A<-rep(0,nvisit*ns) ##place holders for current  treatments
  Ap<-rep(0,nvisit*ns) ##place holders for  previous treatments
  CAp<-rep(0,nvisit*ns)

  Y<-rep(0,nvisit*ns)     ##place holders for outcome vector

  ##Fill in initial values
  seq1<-seq(1,nvisit*ns-nv,nvisit)
  
  X1[seq1]<-0
  X2[seq1]<-0
  X1p[seq1]<-0
  X2p[seq1]<-0
  P1<-list() ##list of treatment probabilities 
  P1[[1]]<-rep(0, ns) 
  seqlist<-list()                              
  seqlist[[1]]<-seq1
  
  for (k in 2:nvisit){  
    ## update covariates
    
    
    seqlist[[k]]<-seqlist[[k-1]]+1
    Ap[seqlist[[k]]]<-A[seqlist[[k-1]]]
    CAp[seqlist[[k]]]<-CAp[seqlist[[k-1]]] + A[seqlist[[k-1]]]
    X1p[seqlist[[k]]]<- X1[seqlist[[k-1]]]
    X2p[seqlist[[k]]]<- X2[seqlist[[k-1]]]
    

    X1[seqlist[[k]]]<-Z1[seqlist[[k]]]-Ap[seqlist[[k]]]-0.5*Ap[seqlist[[k-1]]] + 1.5*X1p[seqlist[[k]]]+ X1p[seqlist[[k-1]]] ## continuous time-varying confounder 
    X2[seqlist[[k]]]<-Z2[seqlist[[k]]]-Ap[seqlist[[k]]]-0.5*Ap[seqlist[[k-1]]] + 1.5*X2p[seqlist[[k]]]+ X2p[seqlist[[k-1]]] ## continuous time-varying confounder 
    
    ## update treatment
    
    lpp<- 5*Ap[seqlist[[k]]] -3*(1-Ap[seqlist[[k]]]) + 0.25*(X1[seqlist[[k]]]+X2[seqlist[[k]]]) +0.3*X1[seqlist[[k]]]*X2[seqlist[[k]]] + 0.25*(X1p[seqlist[[k]]]+X2p[seqlist[[k]]]) 
    P1[[k]]<-1/(1+exp(-lpp))
    
    if (all_treat == TRUE){
      A[seqlist[[k]]]<- 1.0
    } else{ if (all_control == TRUE){
      A[seqlist[[k]]]<- 0.0
    } else{ if (k == 2){
      A[seqlist[[k]]]<-rbinom(ns,1,as.numeric(treat_prev))
    } else{
      A[seqlist[[k]]]<-rbinom(ns,1,P1[[k]])}
    }
    }
    ##Generate outcome
  
    
    Y[seqlist[[k]]]<- 300 + 2*(X1[seqlist[[k]]] + X2[seqlist[[k]]]) + 0.5*(X1p[seqlist[[k]]] + X2p[seqlist[[k]]]) -0.1*X1[seqlist[[k]]]*X2[seqlist[[k]]] -3*A[seqlist[[k]]] -1.5*Ap[seqlist[[k]]] + rnorm(ns,0,5)
    
    
    
  }
  
  ##Make data frame
  
  ID<-rep(1:ns,each=nv)
  
  ##Align data by removing values 
  NSEQ<-seq1
  
  X1<-X1[-NSEQ]
  X2<-X2[-NSEQ]
  X1p<-X1p[-NSEQ]
  X2p<-X2p[-NSEQ]
  A<-A[-NSEQ]
  Ap<-Ap[-NSEQ]
  CAp<-CAp[-NSEQ]
  Y<-Y[-seq(1,nvisit*ns-nv,nvisit)]

  ##Create data frame
  
  DATA<-data.frame(ID,t=rep(c(0:(nv-1)),ns),A,Ap,CAp,X1,X2,X1p,X2p, Y)
  DATA$eligible<-as.numeric(CAp==0)  ## eligibility criteria: age>=18, had no treatment so far, no event so far
  
  ##censoring
  if (censor == T){
    Dprob<-1/(1+exp(2.5 + Ap-0.5*X1-0.2*X2)) ##Probability of dropout
    
    DATA$C<-rbinom(nv*ns,1,Dprob) ##C=0 is remain in the study
    
    
    
    indfun<-function(n){
      if (sum(n)==0) {rep(0,nv)}
      else{k<-min(which(n==1))
      c(rep(0,k),rep(1,nv-k))}}
    
    RL<-ave(DATA$C,DATA$ID,FUN=indfun)
    
    
    eligCum<-ave(DATA$eligible,DATA$ID,FUN=cumsum)
    
    DATA[RL==0 & eligCum>0,] #remove observations after event occurrence and censoring, and not eligible
  } else {
    DATA$C <- 0
    eligCum<-ave(DATA$eligible,DATA$ID,FUN=cumsum)
    
    DATA[eligCum>0,] 
  }
  
}


############################### CHECK DATA SPARSITY ############################
simdata<-DATA_GEN_cont_outcome_treatment_switch_history(ns = as.numeric(size[l]),nv = 5,treat_prev =  0.5,
                                                                 conf =  1.5,
                                                                 censor = F)
simdata <- simdata %>% 
  mutate(switch = ifelse(t == 0, 0,ifelse(A!=Ap,1,0)),
         missX1 = log(abs(X1))/4,
         missX2 = sqrt(abs(X2))/3) %>% 
  group_by(ID) %>% 
  dplyr::mutate(missCX1 = cumsum(missX1),
                missCX2 = cumsum(missX2),
                CA = cumsum(A),
                p = 1/(1+exp(-(2.5*Ap - 2.5*(1-Ap) + 0.5*(X1 + X2 + X1p + X2p))))) %>% 
  dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1)) %>% 
  dplyr::mutate(CRA = cumsum(RA),
                RA = ifelse(CRA == t+1,1,0))
plot(simdata$X2, simdata$p)
con4<-xtabs(~t + A, data=simdata[simdata$RA == 1,])
ftable(con4)

con4<-xtabs(~t + switch, data=simdata)
ftable(con4)

simdata_alltreat <- DATA_GEN_cont_outcome_treatment_switch_history(ns = 1000000,nv = 5,treat_prev =  0.5,
                                                                   conf =  1.5,
                                                                   censor = F, all_treat = TRUE)
simdata_allcontrol <- DATA_GEN_cont_outcome_treatment_switch_history(ns = 1000000,nv = 5,treat_prev =  0.5,
                                                                     conf =  1.5,
                                                                     censor = F, all_control = TRUE)
true_ATE <- mean(simdata_alltreat[simdata_alltreat$t == 4,]$Y)-mean(simdata_allcontrol[simdata_allcontrol$t == 4,]$Y)
true_EY1 <- mean(simdata_alltreat[simdata_alltreat$t == 4,]$Y)
true_EY0 <- mean(simdata_allcontrol[simdata_allcontrol$t == 4,]$Y)
min(simdata$Y)

############################### SIMULATION ####################################

#Simulation setup up:
# Correct specification: treatment or outcome model includes an interaction term
# Misspecification: No interaction term
# Functional misspecification: transformed covariates U1 = log(|X1|)/4, U2 = sqrt(|X2|)/3 (interaction included)

iters <- 500
#l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
size <- c(200,500,1000,2500,5000,10000)

# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 10)
multiResultClass <- function(EY1 = NULL,
                             EY0=NULL,predict_estimates = NULL)
{
  me <- list(
    EY1 = EY1,
    EY0 = EY0,
    predict_estimates = predict_estimates
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}
for (l in 1:6){
  oper <- foreach(i = 1:iters,.combine=cbind) %dopar% {
    tryCatch({
      result <- multiResultClass()
      
      simdata<-DATA_GEN_cont_outcome_treatment_switch_history(ns = as.numeric(size[l]),nv = 5,treat_prev =  0.5,
                                                                       conf =  1.5,
                                                                       censor = F)
      simdata <- simdata %>% 
        mutate(switch = ifelse(t == 0, 0,ifelse(A!=Ap,1,0)),
               missX1 = log(abs(X1))/4,
               missX2 = sqrt(abs(X2))/3,
               missX1p = log(abs(X1p))/4,
               missX2p = sqrt(abs(X2p))/3) %>% 
        mutate(missX1p= ifelse(is.infinite(missX1p), 0, missX1p),
               missX2p= ifelse(is.infinite(missX2p), 0, missX2p)) %>% 
        group_by(ID) %>% 
        dplyr::mutate(CA = cumsum(A))
      #con4<-xtabs(~t + switch, data=simdata)
      #ftable(con4)
      
      
      
      ######### Correctly specified ############### 
      PP_prep <- TrialEmulation::data_preparation(simdata, id='ID', period='t', treatment='A', outcome='Y', cense = 'C',
                                                  eligible ='eligible',
                                                  estimand_type = 'PP',
                                                  switch_d_cov = ~ X1*X2,
                                                  switch_n_cov = ~ -1,
                                                  outcome_cov = ~ X1+X2 + X1p + X2p, model_var = c('assigned_treatment'),
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
      wideSimdata <- data.table::dcast(setDT(simdata), ID ~ t, value.var = c("A", "X1", "X2", "missX1", "missX2", "CA", "Y"))
      
      wideSimdata$y4pred <- wideSimdata$Y_4
      
      
      q4 <- glm(y4pred ~ X1_4*X2_4 + X1_3 + X2_3 + A_4 + A_3,
                data = wideSimdata, family = 'gaussian')
      q4miss <- glm(y4pred ~ X1_4 + X2_4 + X1_3 + X2_3 + A_4 + A_3,
                    data = wideSimdata, family = 'gaussian')
      q4funcmiss <- glm(y4pred ~ missX1_4*missX2_4 + missX1_3 + missX2_3 + A_4 + A_3,
                        data = wideSimdata, family = 'gaussian')
      
      q4ad1 <- glm(y4pred ~ X1_4*X2_4 + X1_3 + X2_3,
                 data = wideSimdata[wideSimdata$A_4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      q4ad0 <- glm(y4pred ~ X1_4*X2_4 + X1_3 + X2_3,
                   data = wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      
      q4admiss1 <- glm(y4pred ~ X1_4+X2_4 + X1_3 + X2_3,
                   data = wideSimdata[wideSimdata$A_4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      q4admiss0 <- glm(y4pred ~ X1_4+X2_4 + X1_3 + X2_3,
                   data = wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      q4adfuncmiss1 <- glm(y4pred ~ missX1_4*missX2_4 + missX1_3 + missX2_3,
                   data = wideSimdata[wideSimdata$A_4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      q4adfuncmiss0 <- glm(y4pred ~ missX1_4*missX2_4 + missX1_3 + missX2_3,
                   data = wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      
      q4SL <- SuperLearner(Y =  wideSimdata$y4pred,
                          X = wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                             "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                             "A_0", "A_1", "A_2", "A_3", "A_4")],
                          family = gaussian(),
                          SL.library = c("SL.glm", "SL.glm.interaction",
                                         "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      q4funcmissSL <- SuperLearner(Y =  wideSimdata$y4pred,
                           X = wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                              "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                              "A_0", "A_1", "A_2", "A_3", "A_4")],
                           family = gaussian(),
                           SL.library = c("SL.glm", "SL.glm.interaction",
                                          "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      
      wideSimdata$y4pred1 <- predict.glm(q4,type = "response", newdata = wideSimdata %>% mutate(A_4 = 1, A_3 = 1))
      wideSimdata$y4predmiss1 <- predict.glm(q4miss,type = "response", newdata = wideSimdata %>% mutate(A_4 = 1, A_3 = 1))
      wideSimdata$y4predfuncmiss1 <- predict.glm(q4funcmiss,type = "response", newdata = wideSimdata %>% mutate(A_4 = 1, A_3 = 1))
      wideSimdata$y4predad1 <- 0.0
      wideSimdata$y4predadmiss1 <- 0.0
      wideSimdata$y4predadfuncmiss1 <- 0.0
      wideSimdata[wideSimdata$A_4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predad1 <- q4ad1$fitted.values
      wideSimdata[wideSimdata$A_4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadmiss1 <- q4admiss1$fitted.values
      wideSimdata[wideSimdata$A_4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadfuncmiss1 <- q4adfuncmiss1$fitted.values

      wideSimdata$y4predSL1 <- predict(q4SL, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                            "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                            "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL1 <- predict(q4funcmissSL, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                    "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                    "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                       onlySL = TRUE)$pred
      
      wideSimdata$y4pred0 <- predict.glm(q4,type = "response", newdata = wideSimdata %>% mutate(A_4 = 0, A_3 = 0))
      wideSimdata$y4predmiss0 <- predict.glm(q4miss,type = "response", newdata = wideSimdata %>% mutate(A_4 = 0, A_3 = 0))
      wideSimdata$y4predfuncmiss0 <- predict.glm(q4funcmiss,type = "response", newdata = wideSimdata %>% mutate(A_4 = 0, A_3 = 0))
      wideSimdata$y4predad0 <- 0.0
      wideSimdata$y4predadmiss0 <- 0.0
      wideSimdata$y4predadfuncmiss0 <- 0.0
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predad0 <- predict.glm(q4ad0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadmiss0 <- predict.glm(q4admiss0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadfuncmiss0 <- predict.glm(q4adfuncmiss0, type = "response") 
      
      wideSimdata$y4predSL0 <- predict(q4SL, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                            "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                            "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL0 <- predict(q4funcmissSL, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                    "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                    "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                               onlySL = TRUE)$pred
      wideSimdata <- wideSimdata %>% 
        mutate(y4predad1 = ifelse(y4predad1 == 0.0, NA,y4predad1 ),
               y4predad0 = ifelse(y4predad0 == 0.0, NA,y4predad0 ),
               y4predadmiss1 = ifelse(y4predadmiss1 == 0.0, NA,y4predadmiss1 ),
               y4predadmiss0 = ifelse(y4predadmiss0 == 0.0, NA,y4predadmiss0 ),
               y4predadfuncmiss1 = ifelse(y4predadfuncmiss1 == 0.0, NA,y4predadfuncmiss1 ),
               y4predadfuncmiss0 = ifelse(y4predadfuncmiss0 == 0.0, NA,y4predadfuncmiss0 ))
      
      q3_1 <- glm(y4pred1 ~ X1_3*X2_3 + X1_2*X2_3 + X1_3*X2_2 + X1_2*X2_2 + A_3 + A_2,
                  data = wideSimdata, family = 'gaussian')
      q3miss_1 <- glm(y4predmiss1 ~ X1_3 + X2_3 + X1_2 + X2_2 + A_3 + A_2,
                      data = wideSimdata, family = 'gaussian')
      q3funcmiss_1 <- glm(y4predfuncmiss1 ~ missX1_3*missX2_3 + missX1_2*missX2_3 + missX1_3*missX2_2 + missX1_2*missX2_2 + A_3 + A_2,
                      data = wideSimdata, family = 'gaussian')
      q3ad1 <- glm(y4predad1 ~ X1_3*X2_3 + X1_2*X2_3 + X1_3*X2_2 + X1_2*X2_2,
                   data = wideSimdata[wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      q3admiss1 <- glm(y4predadmiss1 ~ X1_3 + X2_3 + X1_2 + X2_2,
                   data = wideSimdata[wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      q3adfuncmiss1 <- glm(y4predadfuncmiss1 ~ missX1_3*missX2_3 + missX1_2*missX2_3 + missX1_3*missX2_2 + missX1_2*missX2_2,
                   data = wideSimdata[wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      
      q3_0 <- glm(y4pred0 ~ X1_3*X2_3 + X1_2*X2_3 + X1_3*X2_2 + X1_2*X2_2 + A_3 + A_2,
                  data = wideSimdata, family = 'gaussian')
      q3miss_0 <- glm(y4predmiss0 ~ X1_3 + X2_3 + X1_2 + X2_2 + A_3 + A_2,
                      data = wideSimdata, family = 'gaussian')
      q3funcmiss_0 <- glm(y4predfuncmiss0 ~ missX1_3*missX2_3 + missX1_2*missX2_3 + missX1_3*missX2_2 + missX1_2*missX2_2 + A_3 + A_2,
                          data = wideSimdata, family = 'gaussian')
      q3ad0 <- glm(y4predad0 ~ X1_3*X2_3 + X1_2*X2_3 + X1_3*X2_2 + X1_2*X2_2,
                   data = wideSimdata[wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      q3admiss0 <- glm(y4predadmiss0 ~ X1_3 + X2_3 + X1_2 + X2_2,
                       data = wideSimdata[wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      q3adfuncmiss0 <- glm(y4predadfuncmiss0 ~ missX1_3*missX2_3 + missX1_2*missX2_3 + missX1_3*missX2_2 + missX1_2*missX2_2,
                           data = wideSimdata[wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      
      q3SL1 <- SuperLearner(Y =  wideSimdata$y4predSL1,
                           X = wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3",
                                              "X2_0", "X2_1", "X2_2", "X2_3",
                                              "A_0", "A_1", "A_2", "A_3")],
                           family = gaussian(),
                           SL.library = c("SL.glm", "SL.glm.interaction",
                                          "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      q3funcmissSL1 <- SuperLearner(Y =  wideSimdata$y4predfuncmissSL1,
                                   X = wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3",
                                                      "missX2_0", "missX2_1", "missX2_2", "missX2_3",
                                                      "A_0", "A_1", "A_2", "A_3")],
                                   family = gaussian(),
                                   SL.library = c("SL.glm", "SL.glm.interaction",
                                                  "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      
      q3SL0 <- SuperLearner(Y =  wideSimdata$y4predSL0,
                            X = wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3",
                                               "X2_0", "X2_1", "X2_2", "X2_3",
                                               "A_0", "A_1", "A_2", "A_3")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.glm.interaction",
                                           "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      q3funcmissSL1 <- SuperLearner(Y =  wideSimdata$y4predfuncmissSL0,
                                    X = wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3",
                                                       "missX2_0", "missX2_1", "missX2_2", "missX2_3",
                                                       "A_0", "A_1", "A_2", "A_3")],
                                    family = gaussian(),
                                    SL.library = c("SL.glm", "SL.glm.interaction",
                                                   "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      
      
      wideSimdata$y4pred1 <- predict.glm(q3_1,type = "response", newdata = wideSimdata %>% mutate(A_3 = 1, A_2 = 1))
      wideSimdata$y4predmiss1 <- predict.glm(q3miss_1,type = "response", newdata = wideSimdata %>% mutate(A_3 = 1, A_2 = 1))
      wideSimdata$y4predfuncmiss1 <- predict.glm(q3funcmiss_1,type = "response", newdata = wideSimdata %>% mutate(A_3 = 1, A_2 = 1))
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predad1 <- predict.glm(q3ad1, type = "response")
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadmiss1 <- predict.glm(q3admiss1, type = "response")
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadfuncmiss1 <- predict.glm(q3adfuncmiss1, type = "response")
      
      wideSimdata$y4pred0 <- predict.glm(q3_0,type = "response", newdata = wideSimdata %>% mutate(A_3 = 0, A_2 = 0))
      wideSimdata$y4predmiss0 <- predict.glm(q3miss_0,type = "response", newdata = wideSimdata %>% mutate(A_3 = 0, A_2 = 0))
      wideSimdata$y4predfuncmiss0 <- predict.glm(q3funcmiss_0,type = "response", newdata = wideSimdata %>% mutate(A_3 = 0, A_2 = 0))
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predad0 <- predict.glm(q3ad0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadmiss0 <- predict.glm(q3admiss0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadfuncmiss0 <- predict.glm(q3adfuncmiss0, type = "response") 
      
      wideSimdata$y4predSL1 <- predict(q3SL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                            "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                            "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL1 <- predict(q3funcmissSL1, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                            "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                            "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                               onlySL = TRUE)$pred
      
      wideSimdata$y4predSL0 <- predict(q3SL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                            "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                            "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL0 <- predict(q3funcmissSL0, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                            "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                            "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                               onlySL = TRUE)$pred
      q2_1 <- glm(y4pred1 ~  X1_2*X2_2 + X1_1*X2_2 + X1_2*X2_1 + X1_1*X2_1 + A_2 + A_1,
                  data = wideSimdata, family = 'gaussian')
      q2miss_1 <- glm(y4predmiss1 ~ X1_2 + X2_2 + X1_1 + X2_1  + A_2 + A_1,
                      data = wideSimdata, family = 'gaussian')
      q2funcmiss_1 <- glm(y4predfuncmiss1 ~ missX1_2*missX2_2 + missX1_1*missX2_2 + missX1_2*missX2_1 + missX1_1*missX2_1 + A_2 + A_1,
                      data = wideSimdata, family = 'gaussian')
      q2ad1 <- glm(y4predad1 ~ X1_2*X2_2 + X1_1*X2_2 + X1_2*X2_1 + X1_1*X2_1,
                   data = wideSimdata[wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      q2admiss1 <- glm(y4predadmiss1 ~ X1_2 + X2_2 + X1_1 + X2_1,
                   data = wideSimdata[wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      q2adfuncmiss1 <- glm(y4predadfuncmiss1 ~ missX1_2*missX2_2 + missX1_1*missX2_2 + missX1_2*missX2_1 + missX1_1*missX2_1,
                       data = wideSimdata[wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      
      q2_0 <- glm(y4pred0 ~ X1_2*X2_2 + X1_1*X2_2 + X1_2*X2_1 + X1_1*X2_1 + A_2 + A_1,
                  data = wideSimdata, family = 'gaussian')
      q2miss_0 <- glm(y4predmiss0 ~ X1_2 + X2_2 + X1_1 + X2_1  + A_2 + A_1,
                      data = wideSimdata, family = 'gaussian')
      q2funcmiss_0 <- glm(y4predfuncmiss0 ~ missX1_2*missX2_2 + missX1_1*missX2_2 + missX1_2*missX2_1 + missX1_1*missX2_1 + A_2 + A_1,
                          data = wideSimdata, family = 'gaussian')
      q2ad0 <- glm(y4predad0 ~ X1_2*X2_2 + X1_1*X2_2 + X1_2*X2_1 + X1_1*X2_1,
                   data = wideSimdata[wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      q2admiss0 <- glm(y4predadmiss0 ~ X1_2 + X2_2 + X1_1 + X2_1,
                       data = wideSimdata[wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      q2adfuncmiss0 <- glm(y4predadfuncmiss0 ~ missX1_2*missX2_2 + missX1_1*missX2_2 + missX1_2*missX2_1 + missX1_1*missX2_1,
                           data = wideSimdata[wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      
      q2SL1 <- SuperLearner(Y =  wideSimdata$y4predSL1,
                            X = wideSimdata[,c("X1_0", "X1_1", "X1_2",
                                               "X2_0", "X2_1", "X2_2",
                                               "A_0", "A_1", "A_2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.glm.interaction",
                                           "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      q2funcmissSL1 <- SuperLearner(Y =  wideSimdata$y4predfuncmissSL1,
                                    X = wideSimdata[,c("missX1_0", "missX1_1", "missX1_2",
                                                       "missX2_0", "missX2_1", "missX2_2", 
                                                       "A_0", "A_1", "A_2")],
                                    family = gaussian(),
                                    SL.library = c("SL.glm", "SL.glm.interaction",
                                                   "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      
      q2SL0 <- SuperLearner(Y =  wideSimdata$y4predSL0,
                            X = wideSimdata[,c("X1_0", "X1_1", "X1_2",
                                               "X2_0", "X2_1", "X2_2", 
                                               "A_0", "A_1", "A_2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.glm.interaction",
                                           "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      q2funcmissSL1 <- SuperLearner(Y =  wideSimdata$y4predfuncmissSL0,
                                    X = wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", 
                                                       "missX2_0", "missX2_1", "missX2_2",
                                                       "A_0", "A_1", "A_2")],
                                    family = gaussian(),
                                    SL.library = c("SL.glm", "SL.glm.interaction",
                                                   "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      
      wideSimdata$y4pred1 <- predict.glm(q2_1,type = "response", newdata = wideSimdata %>% mutate(A_2 = 1, A_1 = 1))
      wideSimdata$y4predmiss1 <- predict.glm(q2miss_1,type = "response", newdata = wideSimdata %>% mutate(A_2 = 1, A_1 = 1))
      wideSimdata$y4predfuncmiss1 <- predict.glm(q2funcmiss_1,type = "response", newdata = wideSimdata %>% mutate(A_2 = 1, A_1 = 1))
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predad1 <- predict.glm(q2ad1, type = "response")
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadmiss1 <- predict.glm(q2admiss1, type = "response")
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadfuncmiss1 <- predict.glm(q2adfuncmiss1, type = "response")
      
      
      wideSimdata$y4pred0 <- predict.glm(q2_0,type = "response", newdata = wideSimdata %>% mutate(A_2 = 0, A_1 = 0))
      wideSimdata$y4predmiss0 <- predict.glm(q2miss_0,type = "response", newdata = wideSimdata %>% mutate(A_2 = 0, A_1 = 0))
      wideSimdata$y4predfuncmiss0 <- predict.glm(q2funcmiss_0,type = "response", newdata = wideSimdata %>% mutate(A_2 = 0, A_1 = 0))
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predad0 <- predict.glm(q2ad0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadmiss0 <- predict.glm(q2admiss0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadfuncmiss0 <- predict.glm(q2adfuncmiss0, type = "response") 
      
      wideSimdata$y4predSL1 <- predict(q2SL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                             "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL1 <- predict(q2funcmissSL1, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                             "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                               onlySL = TRUE)$pred
      
      wideSimdata$y4predSL0 <- predict(q2SL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                             "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL0 <- predict(q2funcmissSL0, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                             "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                               onlySL = TRUE)$pred
      
      q1_1 <- glm(y4pred1 ~  X1_1*X2_1 + X1_0*X2_1 + X1_1*X2_0 + X1_0*X2_0 + A_1 + A_0,
                  data = wideSimdata, family = 'gaussian')
      q1miss_1 <- glm(y4predmiss1 ~  X1_1+X2_1 + X1_0 + X2_0 + A_1 + A_0,
                      data = wideSimdata, family = 'gaussian')
      q1funcmiss_1 <- glm(y4predfuncmiss1 ~  missX1_1*missX2_1 + missX1_0*missX2_1 + missX1_1*missX2_0 + missX1_0*missX2_0 + A_1 + A_0,
                      data = wideSimdata, family = 'gaussian')
      q1ad1 <- glm(y4predad1 ~ X1_1*X2_1 + X1_0*X2_1 + X1_1*X2_0 + X1_0*X2_0,
                   data = wideSimdata[wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      q1admiss1 <- glm(y4predadmiss1 ~ X1_1+X2_1 + X1_0 + X2_0,
                       data = wideSimdata[wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      q1adfuncmiss1 <- glm(y4predadfuncmiss1 ~ missX1_1*missX2_1 + missX1_0*missX2_1 + missX1_1*missX2_0 + missX1_0*missX2_0,
                           data = wideSimdata[wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ], family = 'gaussian')
      
      q1_0 <- glm(y4pred0 ~  X1_1*X2_1 + X1_0*X2_1 + X1_1*X2_0 + X1_0*X2_0 + A_1 + A_0,
                  data = wideSimdata, family = 'gaussian')
      q1miss_0 <- glm(y4predmiss0 ~  X1_1+X2_1 + X1_0 + X2_0 + A_1 + A_0,
                      data = wideSimdata, family = 'gaussian')
      q1funcmiss_0 <- glm(y4predfuncmiss0 ~ missX1_1*missX2_1 + missX1_0*missX2_1 + missX1_1*missX2_0 + missX1_0*missX2_0 + A_1 + A_0,
                          data = wideSimdata, family = 'gaussian')
      q1ad0 <- glm(y4predad0 ~ X1_1*X2_1 + X1_0*X2_1 + X1_1*X2_0 + X1_0*X2_0,
                   data = wideSimdata[wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      q1admiss0 <- glm(y4predadmiss0 ~ X1_1+X2_1 + X1_0 + X2_0,
                       data = wideSimdata[wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      q1adfuncmiss0 <- glm(y4predadfuncmiss0 ~ missX1_1*missX2_1 + missX1_0*missX2_1 + missX1_1*missX2_0 + missX1_0*missX2_0,
                           data = wideSimdata[wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ], family = 'gaussian')
      
      q1SL1 <- SuperLearner(Y =  wideSimdata$y4predSL1,
                            X = wideSimdata[,c("X1_0", "X1_1",
                                               "X2_0", "X2_1", 
                                               "A_0", "A_1",)],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.glm.interaction",
                                           "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      q1funcmissSL1 <- SuperLearner(Y =  wideSimdata$y4predfuncmissSL1,
                                    X = wideSimdata[,c("missX1_0", "missX1_1",
                                                       "missX2_0", "missX2_1", 
                                                       "A_0", "A_1")],
                                    family = gaussian(),
                                    SL.library = c("SL.glm", "SL.glm.interaction",
                                                   "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      
      q1SL0 <- SuperLearner(Y =  wideSimdata$y4predSL0,
                            X = wideSimdata[,c("X1_0", "X1_1", 
                                               "X2_0", "X2_1",
                                               "A_0", "A_1")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.glm.interaction",
                                           "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      q1funcmissSL1 <- SuperLearner(Y =  wideSimdata$y4predfuncmissSL0,
                                    X = wideSimdata[,c("missX1_0", "missX1_1", 
                                                       "missX2_0", "missX2_1",
                                                       "A_0", "A_1")],
                                    family = gaussian(),
                                    SL.library = c("SL.glm", "SL.glm.interaction",
                                                   "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      
      
      wideSimdata$y4pred1 <- predict.glm(q1_1,type = "response", newdata = wideSimdata %>% mutate(A_1 = 1, A_0 = 1))
      wideSimdata$y4predmiss1 <- predict.glm(q1miss_1,type = "response", newdata = wideSimdata %>% mutate(A_1 = 1, A_0 = 1))
      wideSimdata$y4predfuncmiss1 <- predict.glm(q1funcmiss_1,type = "response", newdata = wideSimdata %>% mutate(A_1 = 1, A_0 = 1))
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predad1 <- predict.glm(q1ad1, type = "response")
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadmiss1 <- predict.glm(q1admiss1, type = "response")
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadfuncmiss1 <- predict.glm(q1adfuncmiss1, type = "response")
      
      wideSimdata$y4pred0 <- predict.glm(q1_0,type = "response", newdata = wideSimdata %>% mutate(A_1 = 0, A_0 = 0))
      wideSimdata$y4predmiss0 <- predict.glm(q1miss_0,type = "response", newdata = wideSimdata %>% mutate(A_1 = 0, A_0 = 0))
      wideSimdata$y4predfuncmiss0 <- predict.glm(q1funcmiss_0,type = "response", newdata = wideSimdata %>% mutate(A_1 = 0, A_0 = 0))
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predad0 <- predict.glm(q1ad0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadmiss0 <- predict.glm(q1admiss0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadfuncmiss0 <- predict.glm(q1adfuncmiss0, type = "response") 
      
      wideSimdata$y4predSL1 <- predict(q1SL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                             "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL1 <- predict(q1funcmissSL1, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                             "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                               onlySL = TRUE)$pred
      
      wideSimdata$y4predSL0 <- predict(q1SL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                             "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL0 <- predict(q1funcmissSL0, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                             "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                               onlySL = TRUE)$pred
      
      
      q0_1 <- glm(y4pred1 ~  X1_0*X2_0 + A_0,
                  data = wideSimdata, family = 'gaussian')
      q0miss_1 <- glm(y4predmiss1 ~  X1_0+X2_0 + A_0,
                      data = wideSimdata, family = 'gaussian')
      q0funcmiss_1 <- glm(y4predfuncmiss1 ~  missX1_0*missX2_0 + A_0,
                      data = wideSimdata, family = 'gaussian')
      q0ad1 <- glm(y4predad1 ~ X1_0*X2_0,
                   data = wideSimdata[wideSimdata$A_0 == 1, ], family = 'gaussian')
      q0admiss1 <- glm(y4predadmiss1 ~ X1_0+X2_0,
                   data = wideSimdata[wideSimdata$A_0 == 1, ], family = 'gaussian')
      q0adfuncmiss1 <- glm(y4predadfuncmiss1 ~ missX1_0*missX2_0,
                   data = wideSimdata[wideSimdata$A_0 == 1, ], family = 'gaussian')
      
      q0_0 <- glm(y4pred0 ~  X1_0*X2_0 + A_0,
                  data = wideSimdata, family = 'gaussian')
      q0miss_0 <- glm(y4predmiss0 ~  X1_0+X2_0 + A_0,
                      data = wideSimdata, family = 'gaussian')
      q0funcmiss_0 <- glm(y4predfuncmiss0 ~  missX1_0*missX2_0 + A_0,
                          data = wideSimdata, family = 'gaussian')
      q0ad0 <- glm(y4predad0 ~ X1_0*X2_0,
                   data = wideSimdata[wideSimdata$A_0 == 0, ], family = 'gaussian')
      q0admiss0 <- glm(y4predadmiss0 ~ X1_0+X2_0,
                       data = wideSimdata[wideSimdata$A_0 == 0, ], family = 'gaussian')
      q0adfuncmiss0 <- glm(y4predadfuncmiss0 ~ missX1_0*missX2_0,
                           data = wideSimdata[wideSimdata$A_0 == 0, ], family = 'gaussian')
      q0SL1 <- SuperLearner(Y =  wideSimdata$y4predSL1,
                            X = wideSimdata[,c("X1_0",
                                               "X2_0",
                                               "A_0")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.glm.interaction",
                                           "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      q0funcmissSL1 <- SuperLearner(Y =  wideSimdata$y4predfuncmissSL1,
                                    X = wideSimdata[,c("missX1_0",
                                                       "missX2_0", 
                                                       "A_0")],
                                    family = gaussian(),
                                    SL.library = c("SL.glm", "SL.glm.interaction",
                                                   "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      
      q0SL0 <- SuperLearner(Y =  wideSimdata$y4predSL0,
                            X = wideSimdata[,c("X1_0",  
                                               "X2_0", 
                                               "A_0")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.glm.interaction",
                                           "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      q0funcmissSL1 <- SuperLearner(Y =  wideSimdata$y4predfuncmissSL0,
                                    X = wideSimdata[,c("missX1_0", 
                                                       "missX2_0", 
                                                       "A_0")],
                                    family = gaussian(),
                                    SL.library = c("SL.glm", "SL.glm.interaction",
                                                   "SL.bayesglm","SL.stepAIC","SL.step.interaction","SL.glmnet","SL.ridge"))
      
      wideSimdata$y4pred1 <- predict.glm(q0_1,type = "response", newdata = wideSimdata %>% mutate(A_0 = 1))
      wideSimdata$y4predmiss1 <- predict.glm(q0miss_1,type = "response", newdata = wideSimdata %>% mutate(A_0 = 1))
      wideSimdata$y4predfuncmiss1 <- predict.glm(q0funcmiss_1,type = "response", newdata = wideSimdata %>% mutate(A_0 = 1))
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predad1 <- predict.glm(q0ad1, type = "response")
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadmiss1 <- predict.glm(q0admiss1, type = "response")
      wideSimdata[wideSimdata$A4 == 1 & wideSimdata$A_3 == 1 & wideSimdata$A_2 == 1& wideSimdata$A_1 == 1& wideSimdata$A_0 == 1, ]$y4predadfuncmiss1 <- predict.glm(q0adfuncmiss1, type = "response")
      
      wideSimdata$y4pred0 <- predict.glm(q0_0,type = "response", newdata = wideSimdata %>% mutate(A_0 = 0))
      wideSimdata$y4predmiss0 <- predict.glm(q0miss_0,type = "response", newdata = wideSimdata %>% mutate(A_0 = 0))
      wideSimdata$y4predfuncmiss0 <- predict.glm(q0funcmiss_0,type = "response", newdata = wideSimdata %>% mutate(A_0 = 0))
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predad0 <- predict.glm(q0ad0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadmiss0 <- predict.glm(q0admiss0, type = "response") 
      wideSimdata[wideSimdata$A_4 == 0 & wideSimdata$A_3 == 0 & wideSimdata$A_2 == 0& wideSimdata$A_1 == 0& wideSimdata$A_0 == 0, ]$y4predadfuncmiss0 <- predict.glm(q0adfuncmiss0, type = "response") 
      
      wideSimdata$y4predSL1 <- predict(q0SL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                             "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL1 <- predict(q0funcmissSL1, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                             "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                               onlySL = TRUE)$pred
      
      wideSimdata$y4predSL0 <- predict(q0SL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                             "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                       onlySL = TRUE)$pred
      wideSimdata$y4predfuncmissSL0 <- predict(q0funcmissSL0, wideSimdata[,c("missX1_0", "missX1_1", "missX1_2", "missX1_3","missX1_4",
                                                                             "missX2_0", "missX2_1", "missX2_2", "missX2_3", "missX2_4",
                                                                             "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                               onlySL = TRUE)$pred
      
      
      
      
      
      wideSimdata$Y1_0 <- predict.glm(q0_1, type = "response", newdata = wideSimdata %>% mutate(A_0 = 1))
      wideSimdata$Y1_1 <- predict.glm(q1_1, type = "response", newdata = wideSimdata %>% mutate(A_1 = 1, A_0 = 1))
      wideSimdata$Y1_2 <- predict.glm(q2_1, type = "response", newdata = wideSimdata%>% mutate(A_2 = 1, A_1 = 1))
      wideSimdata$Y1_3 <- predict.glm(q3_1, type = "response", newdata = wideSimdata%>% mutate(A_3 = 1, A_2 = 1))
      wideSimdata$Y1_4 <- predict.glm(q4, type = "response", newdata = wideSimdata %>% mutate(A_4 = 1, A_3 = 1))
      
      GcompY1 <- mean(q0ad1$fitted.values)
      
      wideSimdata$missY1_0 <- predict.glm(q0miss_1, type = "response", newdata = wideSimdata %>% mutate(A_0 = 1))
      wideSimdata$missY1_1 <- predict.glm(q1miss_1, type = "response", newdata = wideSimdata %>% mutate(A_1 = 1, A_0 = 1))
      wideSimdata$missY1_2 <- predict.glm(q2miss_1, type = "response", newdata = wideSimdata%>% mutate(A_2 = 1, A_1 = 1))
      wideSimdata$missY1_3 <- predict.glm(q3miss_1, type = "response", newdata = wideSimdata%>% mutate(A_3 = 1, A_2 = 1))
      wideSimdata$missY1_4 <- predict.glm(q4miss, type = "response", newdata = wideSimdata %>% mutate(A_4 = 1, A_3 = 1))
      
      GcompmissY1 <- mean(q0admiss1$fitted.values)
      
      wideSimdata$funcmissY1_0 <- predict.glm(q0funcmiss_1, type = "response", newdata = wideSimdata %>% mutate(A_0 = 1))
      wideSimdata$funcmissY1_1 <- predict.glm(q1funcmiss_1, type = "response", newdata = wideSimdata %>% mutate(A_1 = 1, A_0 = 1))
      wideSimdata$funcmissY1_2 <- predict.glm(q2funcmiss_1, type = "response", newdata = wideSimdata%>% mutate(A_2 = 1, A_1 = 1))
      wideSimdata$funcmissY1_3 <- predict.glm(q3funcmiss_1, type = "response", newdata = wideSimdata%>% mutate(A_3 = 1, A_2 = 1))
      wideSimdata$funcmissY1_4 <- predict.glm(q4funcmiss, type = "response", newdata = wideSimdata %>% mutate(A_4 = 1, A_3 = 1))
      
      GcompfuncmissY1 <- mean(q0adfuncmiss1$fitted.values)
      
      
      wideSimdata$SLY1_0 <- predict(q0SL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$SLY1_1 <- predict(q1SL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$SLY1_2 <- predict(q2SL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$SLY1_3 <- predict(q3SL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$SLY1_4 <- predict(q4SL, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$funcmissSLY1_0 <- predict(q0funcmissSL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$funcmissSLY1_1 <- predict(q1funcmissSL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$funcmissSLY1_2 <- predict(q2funcmissSL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$funcmissSLY1_3 <- predict(q3funcmissSL1, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$funcmissSLY1_4 <- predict(q4funcmissSL, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 1, A_3 = 1,A_2 = 1, A_1 = 1,A_0 = 1),
                                    onlySL = TRUE)$pred
      wideSimdata$Y0_0 <- predict.glm(q0_0, type = "response", newdata = wideSimdata%>% mutate(A_0 = 0))
      wideSimdata$Y0_1 <- predict.glm(q1_0, type = "response", newdata = wideSimdata%>% mutate(A_1 = 0, A_0 = 0))
      wideSimdata$Y0_2 <- predict.glm(q2_0, type = "response", newdata = wideSimdata%>% mutate(A_2 = 0, A_1 = 0))
      wideSimdata$Y0_3 <- predict.glm(q3_0, type = "response", newdata = wideSimdata%>% mutate(A_3 = 0, A_2 = 0))
      wideSimdata$Y0_4 <- predict.glm(q4, type = "response", newdata = wideSimdata %>% mutate(A_4 = 0, A_3 = 0))
      
      GcompY0 <- mean(q0ad0$fitted.values)
      
      wideSimdata$missY0_0 <- predict.glm(q0miss_0, type = "response", newdata = wideSimdata%>% mutate(A_0 = 0))
      wideSimdata$missY0_1 <- predict.glm(q1miss_0, type = "response", newdata = wideSimdata%>% mutate(A_1 = 0, A_0 = 0))
      wideSimdata$missY0_2 <- predict.glm(q2miss_0, type = "response", newdata = wideSimdata%>% mutate(A_2 = 0, A_1 = 0))
      wideSimdata$missY0_3 <- predict.glm(q3miss_0, type = "response", newdata = wideSimdata%>% mutate(A_3 = 0, A_2 = 0))
      wideSimdata$missY0_4 <- predict.glm(q4miss, type = "response", newdata = wideSimdata %>% mutate(A_4 = 0, A_3 = 0))
      
      GcompmissY0 <- mean(q0admiss0$fitted.values)
      
      wideSimdata$funcmissY0_0 <- predict.glm(q0funcmiss_0, type = "response", newdata = wideSimdata%>% mutate(A_0 = 0))
      wideSimdata$funcmissY0_1 <- predict.glm(q1funcmiss_0, type = "response", newdata = wideSimdata%>% mutate(A_1 = 0, A_0 = 0))
      wideSimdata$funcmissY0_2 <- predict.glm(q2funcmiss_0, type = "response", newdata = wideSimdata%>% mutate(A_2 = 0, A_1 = 0))
      wideSimdata$funcmissY0_3 <- predict.glm(q3funcmiss_0, type = "response", newdata = wideSimdata%>% mutate(A_3 = 0, A_2 = 0))
      wideSimdata$funcmissY0_4 <- predict.glm(q4funcmiss, type = "response", newdata = wideSimdata %>% mutate(A_4 = 0, A_3 = 0))
      
      GcompfuncmissY0 <- mean(q0adfuncmiss0$fitted.values)
      
      wideSimdata$SLY0_0 <- predict(q0SL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 1,A_0 = 0),
                                    onlySL = TRUE)$pred
      wideSimdata$SLY0_1 <- predict(q1SL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                    onlySL = TRUE)$pred
      wideSimdata$SLY0_2 <- predict(q2SL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                    onlySL = TRUE)$pred
      wideSimdata$SLY0_3 <- predict(q3SL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                    onlySL = TRUE)$pred
      wideSimdata$SLY0_4 <- predict(q4SL, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                         "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                         "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                    onlySL = TRUE)$pred
      wideSimdata$funcmissSLY0_0 <- predict(q0funcmissSL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                            onlySL = TRUE)$pred
      wideSimdata$funcmissSLY0_1 <- predict(q1funcmissSL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                            onlySL = TRUE)$pred
      wideSimdata$funcmissSLY0_2 <- predict(q2funcmissSL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                            onlySL = TRUE)$pred
      wideSimdata$funcmissSLY0_3 <- predict(q3funcmissSL0, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                                          "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                                          "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                            onlySL = TRUE)$pred
      wideSimdata$funcmissSLY0_4 <- predict(q4funcmissSL, wideSimdata[,c("X1_0", "X1_1", "X1_2", "X1_3","X1_4",
                                                                         "X2_0", "X2_1", "X2_2", "X2_3", "X2_4",
                                                                         "A_0", "A_1", "A_2", "A_3", "A_4")] %>% mutate(A_4 = 0, A_3 = 0,A_2 = 0, A_1 = 0,A_0 = 0),
                                            onlySL = TRUE)$pred
      simdata_imputed <- pivot_longer(wideSimdata %>% dplyr::select(ID, Y1_0, Y1_1, Y1_2, Y1_3, Y1_4),
                                      cols = Y1_0:Y1_4, names_to = 't', names_prefix = "Y1_", values_to = "Y1hat") %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, Y0_0, Y0_1, Y0_2, Y0_3, Y0_4),
                           cols = Y0_0:Y0_4, names_to = 't', names_prefix = "Y0_", values_to = "Y0hat"),by = c("ID", "t")) %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, missY0_0, missY0_1, missY0_2, missY0_3, missY0_4),
                           cols = missY0_0:missY0_4, names_to = 't', names_prefix = "missY0_", values_to = "missY0hat"),by = c("ID", "t")) %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, missY1_0, missY1_1, missY1_2, missY1_3, missY1_4),
                           cols = missY1_0:missY1_4, names_to = 't', names_prefix = "missY1_", values_to = "missY1hat"),by = c("ID", "t")) %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, funcmissY0_0, funcmissY0_1, funcmissY0_2, funcmissY0_3, funcmissY0_4),
                           cols = funcmissY0_0:funcmissY0_4, names_to = 't', names_prefix = "funcmissY0_", values_to = "funcmissY0hat"),by = c("ID", "t")) %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, funcmissY1_0, funcmissY1_1, funcmissY1_2, funcmissY1_3, funcmissY1_4),
                           cols = funcmissY1_0:funcmissY1_4, names_to = 't', names_prefix = "funcmissY1_", values_to = "funcmissY1hat"),by = c("ID", "t")) %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, funcmissSLY0_0, funcmissSLY0_1, funcmissSLY0_2, funcmissSLY0_3, funcmissSLY0_4),
                           cols = funcmissSLY0_0:funcmissSLY0_4, names_to = 't', names_prefix = "funcmissSLY0_", values_to = "funcmissSLY0hat"),by = c("ID", "t")) %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, funcmissSLY1_0, funcmissSLY1_1, funcmissSLY1_2, funcmissSLY1_3, funcmissSLY1_4),
                           cols = funcmissSLY1_0:funcmissSLY1_4, names_to = 't', names_prefix = "funcmissSLY1_", values_to = "funcmissSLY1hat"),by = c("ID", "t")) %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, SLY0_0, SLY0_1, SLY0_2, SLY0_3, SLY0_4),
                           cols = SLY0_0:SLY0_4, names_to = 't', names_prefix = "SLY0_", values_to = "SLY0hat"),by = c("ID", "t")) %>% 
        merge(pivot_longer(wideSimdata %>% dplyr::select(ID, SLY1_0, SLY1_1, SLY1_2, SLY1_3, SLY1_4),
                           cols = SLY1_0:SLY1_4, names_to = 't', names_prefix = "SLY1_", values_to = "SLY1hat"),by = c("ID", "t")) %>% 
        dplyr::mutate(t = as.numeric(t))
      
      
      ########### Data prep for calibration ###############
      
      data_restric <- simdata %>% 
        merge(simdata_imputed, by = c("ID", "t")) %>% 
        dplyr::group_by(ID) %>% 
        dplyr::mutate(A_0 = first(A), RA = ifelse(t !=0, ifelse(A == first(A) & A==Ap, 1, 0),1)) %>% 
        dplyr::mutate(CRA = cumsum(RA),
                      RA = ifelse(CRA == t+1,1,0),
                      RC = ifelse(lag(C) == 0,1,0),
                      A1X2 = A_0*X2,
                      A0X2 = (1-A_0)*X2,
                      A1X1 = A_0*X1,
                      A0X1 = (1-A_0)*X1,
                      A1X1X2 = A_0*X1*X2,
                      A0X1X2 = (1-A_0)*X1*X2,
                      A1missX2 = A_0*missX2,
                      A0missX2 = (1-A_0)*missX2,
                      A1missX1 = A_0*missX1,
                      A0missX1 = (1-A_0)*missX1,
                      A1missX1X2 = A_0*missX1*missX2,
                      A0missX1X2 = (1-A_0)*missX1*missX2,
                      A1 = A_0,
                      A0 = 1-A_0,
                      A1Y1hat = A_0*Y1hat,
                      A0Y0hat = (1-A_0)*Y0hat,
                      A1missY1hat = A_0*missY1hat,
                      A0missY0hat = (1-A_0)*missY0hat,
                      A1funcmissY1hat = A_0*funcmissY1hat,
                      A0funcmissY0hat = (1-A_0)*funcmissY0hat,
                      A1SLY1hat = A_0*SLY1hat,
                      A0SLY0hat = (1-A_0)*SLY0hat,
                      A1funcmissSLY1hat = A_0*funcmissSLY1hat,
                      A0funcmissSLY0hat = (1-A_0)*funcmissSLY0hat,
                      sub = ID,
                      tall = t,
                      One = 1.0) %>% 
        merge(dplyr::select(switch_data,id, followup_time, weight), 
              by.x = c('ID', 't'), by.y = c('id', 'followup_time'), all.x = T) %>% 
        dplyr::mutate(weights = ifelse(!is.na(weight), weight, 0)) %>% 
        dplyr::arrange(ID, t) 
      
      ############# IPW estimation ################
      weight_model <- glm(data = data_restric[data_restric$t != 0,], 
                          formula = A ~ Ap + X1*X2 + X1p + X2p , family = 'binomial')
      weight_model_miss <- glm(data = data_restric[data_restric$t != 0,], 
                               formula = A ~ Ap + X1 + X2 + X1p + X2p , family = 'binomial')
      weight_model_funcmiss <- glm(data = data_restric[data_restric$t != 0,], 
                               formula = A ~ Ap + missX1*missX2 + missX1p + missX2p , family = 'binomial')
      
      data_restric$p_1 <- 1.0
      data_restric$p_1miss <- 1.0
      data_restric$p_1funcmiss <- 1.0
      data_restric[data_restric$t != 0,]$p_1 <- predict.glm(weight_model, data_restric[data_restric$t != 0,], type = 'response')
      data_restric[data_restric$t != 0,]$p_1miss <- predict.glm(weight_model_miss, data_restric[data_restric$t != 0,], type = 'response')
      data_restric[data_restric$t != 0,]$p_1funcmiss <- predict.glm(weight_model_funcmiss, data_restric[data_restric$t != 0,], type = 'response')
      
      data_restric <- data_restric %>% 
        arrange(ID, t) %>% 
        group_by(ID) %>% 
        dplyr::mutate(
          wt = ifelse( t == 0, 1.0,ifelse(A == 1, 1/p_1, 1/(1-p_1))),
          wtprod = cumprod(wt),
          weights = ifelse(weights !=0.0,wtprod,0.0),
          wtmiss = ifelse( t == 0, 1.0,ifelse(A == 1, 1/p_1miss, 1/(1-p_1miss))),
          wtprodmiss = cumprod(wtmiss),
          weights_miss = ifelse(weights !=0.0,wtprodmiss,0.0),
          wtfuncmiss = ifelse( t == 0, 1.0,ifelse(A == 1, 1/p_1funcmiss, 1/(1-p_1funcmiss))),
          wtprodfuncmiss = cumprod(wtfuncmiss),
          weights_funcmiss = ifelse(weights !=0.0,wtprodfuncmiss,0.0))
      
      ################### Calibration by time Lk #######################
      simdatafinal1 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1X1', 'A1X2','A1X1X2',
                                                   'A0','A0X1','A0X2','A0X1X2'))
      
      
      
      ################### Calibration by time g(Lk)  #######################
      simdatafinal2 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1Y1hat',
                                                   'A0','A0Y0hat'))
      
      
      ################## Calibration by time Lk, g(Lk) ###########################
      simdatafinal3 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1X1', 'A1X2','A1X1X2','A1Y1hat',
                                                   'A0','A0X1','A0X2','A0X1X2','A0Y0hat'))
      
      
      
      ################## Calibration by time g(Lk) misspecified, IPW correct #######################
      simdatafinal4 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1missY1hat',
                                                   'A0','A0missY0hat'))
      
      ################## Calibration by time Lk correct, g(Lk) misspecified ###########################
      simdatafinal5 <- calibration_by_time(simdatafinal = data_restric, 
                                           var = c('A1','A1X1', 'A1X2','A1X1X2','A1missY1hat',
                                                   'A0','A0X1','A0X2','A0X1X2','A0missY0hat'))
      
      ################### Calibration by time Lk miss #######################
      simdatafinal6 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1X1', 'A1X2',
                                                   'A0','A0X1', 'A0X2'))
      ################### Calibration by time Lk miss, g(Lk) only but correct #######################
      simdatafinal7 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1Y1hat',
                                                   'A0','A0Y0hat'))
      ################### Calibration by time Lk miss with g(Lk) correct #######################
      simdatafinal8 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1X1', 'A1X2', 'A1Y1hat',
                                                   'A0','A0X1','A0X2', 'A0Y0hat'))
      
      ################### Calibration by time Lk miss, g(Lk) miss only #######################
      simdatafinal9 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                           var = c('A1','A1missY1hat',
                                                   'A0','A0missY0hat'))
      
      
      ################### Calibration by time all miss #######################
      simdatafinal10 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss), 
                                            var = c('A1','A1X1', 'A1X2', 'A1missY1hat',
                                                    'A0','A0X1', 'A0X2','A0missY0hat'))
      
      
      ################## Calibration by time func miss Lk only ###########################
      simdatafinal11 <- calibration_by_time(simdatafinal = data_restric %>%  mutate(weights = weights_funcmiss), 
                                            var = c('A1','A1missX1', 'A1missX2','A1missX1X2',
                                                    'A0','A0missX1', 'A0missX2','A0missX1X2'))
      
      
      
      ################## Calibration by time func miss g(Lk) only #######################
      simdatafinal12 <- calibration_by_time(simdatafinal = data_restric %>%  mutate(weights = weights_funcmiss), 
                                            var = c('A1','A1funcmissY1hat',
                                                    'A0','A0funcmissY0hat'))
      
      ################## Calibration by time func miss Lk, g(Lk) ###########################
      simdatafinal13 <- calibration_by_time(simdatafinal = data_restric %>%  mutate(weights = weights_funcmiss), 
                                            var = c('A1','A1missX1', 'A1missX2','A1missX1X2','A1funcmissY1hat',
                                                    'A0','A0missX1', 'A0missX2','A0missX1X2','A0funcmissY0hat'))
      
      ################## Calibration by time g(Lk) only SL###########################
      simdatafinal14 <- calibration_by_time(simdatafinal = data_restric,
                                            var = c('A1','A1SLY1hat',
                                                    'A0','A0SLY0hat'))
      
      ################## Calibration by time Lk, g(Lk)  SL###########################
      simdatafinal15 <- calibration_by_time(simdatafinal = data_restric,
                                            var = c('A1','A1X1', 'A1X2', 'A1X1X2', 'A1SLY1hat',
                                                    'A0','A0X1', 'A0X2', 'A0X1X2','A0SLY0hat'))
      
      ################## Calibration by time miss, g(Lk) only SL###########################
      simdatafinal16 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss),
                                            var = c('A1','A1SLY1hat',
                                                    'A0','A0SLY0hat'))
      
      ################## Calibration by time miss Lk, g(Lk)  SL###########################
      simdatafinal17 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss),
                                            var = c('A1','A1X1', 'A1X2',  'A1SLY1hat',
                                                    'A0','A0X1', 'A0X2', 'A0SLY0hat'))
      
      ################## Calibration by time func miss, g(Lk) only SL###########################
      simdatafinal18 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_funcmiss),
                                            var = c('A1','A1funcmissSLY1hat',
                                                    'A0','A0funcmissSLY0hat'))
      
      ################## Calibration by time funcmiss Lk, g(Lk)  SL###########################
      simdatafinal19 <- calibration_by_time(simdatafinal = data_restric %>% mutate(weights = weights_miss),
                                            var = c('A1','A1missX1', 'A1missX2','A1missX1X2','A1funcmissSLY1hat',
                                                    'A0','A0missX1', 'A0missX2','A0missX1X2','A0funcmissSLY0hat'))
      
      switch_data$IPWcorrect <- simdatafinal1$data[simdatafinal1$data$RA == 1,]$weights
      switch_data$CaliLkcorrect <- simdatafinal1$data[simdatafinal1$data$RA == 1,]$Cweights
      switch_data$CaliGLkcorrect <- simdatafinal2$data[simdatafinal2$data$RA == 1,]$Cweights
      switch_data$CaliAllCorrect <- simdatafinal3$data[simdatafinal3$data$RA == 1,]$Cweights
      switch_data$CaliGLkmiss <- simdatafinal4$data[simdatafinal4$data$RA == 1,]$Cweights
      switch_data$CaliAllGLkmiss <- simdatafinal5$data[simdatafinal5$data$RA == 1,]$Cweights
      switch_data$CaliLkmiss <- simdatafinal6$data[simdatafinal6$data$RA == 1,]$Cweights
      switch_data$CaliMissGLkcorrect <- simdatafinal7$data[simdatafinal7$data$RA == 1,]$Cweights
      switch_data$CaliMissAllGLkcorrect <- simdatafinal8$data[simdatafinal8$data$RA == 1,]$Cweights
      switch_data$CaliMissGLKmiss <- simdatafinal9$data[simdatafinal9$data$RA == 1,]$Cweights
      switch_data$IPWmiss <- simdatafinal9$data[simdatafinal9$data$RA == 1,]$weights
      switch_data$CaliAllmiss <- simdatafinal10$data[simdatafinal10$data$RA == 1,]$Cweights
      switch_data$CalifuncmissLk <- simdatafinal11$data[simdatafinal11$data$RA == 1,]$Cweights
      switch_data$CalifuncmissGLk <- simdatafinal12$data[simdatafinal12$data$RA == 1,]$Cweights
      switch_data$CalifuncmissAll <- simdatafinal13$data[simdatafinal13$data$RA == 1,]$Cweights
      switch_data$IPWfuncmiss <- simdatafinal11$data[simdatafinal11$data$RA == 1,]$weights
      switch_data$CaliCorrectGLkSL <- simdatafinal14$data[simdatafinal14$data$RA == 1,]$Cweights
      switch_data$CaliAllCorrectSL <- simdatafinal15$data[simdatafinal15$data$RA == 1,]$Cweights
      switch_data$CaliMissGLkSL <- simdatafinal16$data[simdatafinal16$data$RA == 1,]$Cweights
      switch_data$CaliMissAllSL <- simdatafinal17$data[simdatafinal17$data$RA == 1,]$Cweights
      switch_data$CalifuncmissGLkSL <- simdatafinal18$data[simdatafinal18$data$RA == 1,]$Cweights
      switch_data$CalifuncmissAllSL <- simdatafinal19$data[simdatafinal19$data$RA == 1,]$Cweights
      
      MSM_data  <- wideSimdata <- data.table::dcast(setDT(switch_data), id ~ t, value.var = c("assigned_treatment", "Y",'IPWcorrect', 
                                                                                              'CaliLkcorrect', 'CaliGLkcorrect', 
                                                                                              'CaliAllCorrect', 'CaliGLkmiss', 
                                                                                              'CaliAllGLkmiss', 'CaliLkmiss', 
                                                                                              'CaliMissGLkcorrect', 'CaliMissAllGLkcorrect',
                                                                                              'CaliMissGLKmiss', 'IPWmiss', 
                                                                                              'CaliAllmiss', 'CalifuncmissLk', 
                                                                                              'CalifuncmissGLk', 
                                                                                              'CalifuncmissAll', 'IPWfuncmiss',
                                                                                              'CaliCorrectGLkSL','CaliAllCorrectSL',
                                                                                              'CaliMissGLkSL', 'CaliMissAllSL',
                                                                                              'CalifuncmissGLkSL','CalifuncmissAllSL'
 ))
      
      PP_naive <- glm(data = MSM_data,
                           formula = Y_4 ~ assigned_treatment_0 ,
                           weights = NULL, family = 'gaussian')
      summary(PP_naive)
      
      PP_IPWcorrect <- glm(data = MSM_data,
                           formula = Y_4 ~ assigned_treatment_0 ,
                           weights = IPWcorrect_4, family = 'gaussian')
      summary(PP_IPWcorrect)
      
      PP_CaliLkcorrect <- glm(data = MSM_data,
                              formula = Y_4 ~ assigned_treatment_0 ,
                              weights = CaliLkcorrect_4, family = 'gaussian')
      summary(PP_CaliLkcorrect)
      
      PP_CaliGLkcorrect <- glm(data = MSM_data,
                               formula = Y_4 ~ assigned_treatment_0 ,
                               weights = CaliGLkcorrect_4, family = 'gaussian')
      summary(PP_CaliGLkcorrect)
      
      PP_CaliAllCorrect <- glm(data = MSM_data,
                               formula = Y_4 ~ assigned_treatment_0 ,
                               weights = CaliAllCorrect_4, family = 'gaussian')
      summary(PP_CaliAllCorrect)
      
      PP_CaliGLkmiss <- glm(data = MSM_data,
                            formula = Y_4 ~ assigned_treatment_0 ,
                            weights = CaliGLkmiss_4, family = 'gaussian')
      summary(PP_CaliGLkmiss)
      
      PP_CaliAllGLkmiss <- glm(data = MSM_data,
                               formula = Y_4 ~ assigned_treatment_0 ,
                               weights = CaliAllGLkmiss_4, family = 'gaussian')
      summary(PP_CaliAllGLkmiss)
      
      PP_CaliLkmiss <- glm(data = MSM_data,
                           formula = Y_4 ~ assigned_treatment_0 ,
                           weights = CaliLkmiss_4, family = 'gaussian')
      summary(PP_CaliLkmiss)
      
      PP_CaliMissGLkcorrect <- glm(data = MSM_data,
                                   formula = Y_4 ~ assigned_treatment_0 ,
                                   weights = CaliMissGLkcorrect_4, family = 'gaussian')
      summary(PP_CaliMissGLkcorrect)
      
      PP_CaliMissAllGLkcorrect <- glm(data = MSM_data,
                                      formula = Y_4 ~ assigned_treatment_0 ,
                                      weights = CaliMissAllGLkcorrect_4, family = 'gaussian')
      summary(PP_CaliMissAllGLkcorrect)
      
      PP_CaliMissGLKmiss <- glm(data = MSM_data,
                                formula = Y_4 ~ assigned_treatment_0 ,
                                weights = CaliMissGLKmiss_4, family = 'gaussian')
      summary(PP_CaliMissGLKmiss)
      
      PP_IPWmiss <- glm(data = MSM_data,
                        formula = Y_4 ~ assigned_treatment_0 ,
                        weights = IPWmiss_4, family = 'gaussian')
      summary(PP_IPWmiss)
      
      PP_CaliAllmiss <- glm(data = MSM_data,
                            formula = Y_4 ~ assigned_treatment_0 ,
                            weights = CaliAllmiss_4, family = 'gaussian')
      summary(PP_CaliAllmiss)
      
      PP_CalifuncmissLk <- glm(data = MSM_data,
                               formula = Y_4 ~ assigned_treatment_0 ,
                               weights = CalifuncmissLk_4, family = 'gaussian')
      summary(PP_CalifuncmissLk)
      
      PP_CalifuncmissGLk <- glm(data = MSM_data,
                                formula = Y_4 ~ assigned_treatment_0 ,
                                weights = CalifuncmissGLk_4, family = 'gaussian')
      summary(PP_CalifuncmissGLk)
      
      PP_CalifuncmissAll <- glm(data = MSM_data,
                                formula = Y_4 ~ assigned_treatment_0 ,
                                weights = CalifuncmissAll_4, family = 'gaussian')
      summary(PP_CalifuncmissAll)
      
      PP_IPWfuncmiss <- glm(data = MSM_data,
                        formula = Y_4 ~ assigned_treatment_0 ,
                        weights = IPWfuncmiss_4, family = 'gaussian')
      summary(PP_IPWfuncmiss)
      
      PP_CaliCorrectGLkSL <- glm(data = MSM_data,
                            formula = Y_4 ~ assigned_treatment_0 ,
                            weights = CaliCorrectGLkSL_4, family = 'gaussian')
      summary(PP_CaliCorrectGLkSL)
      
      PP_CaliAllCorrectSL <- glm(data = MSM_data,
                                 formula = Y_4 ~ assigned_treatment_0 ,
                                 weights = CaliAllCorrectSL_4, family = 'gaussian')
      summary(PP_CaliAllCorrectSL)
      
      PP_CaliMissGLkSL <- glm(data = MSM_data,
                                 formula = Y_4 ~ assigned_treatment_0 ,
                                 weights = CaliMissGLkSL_4, family = 'gaussian')
      summary(PP_CaliMissGLkSL)
      
      PP_CaliMissAllSL <- glm(data = MSM_data,
                             formula = Y_4 ~ assigned_treatment_0 ,
                             weights = CaliMissAllSL_4, family = 'gaussian')
      summary(PP_CaliMissAllSL)
      
      PP_CalifuncmissGLkSL <- glm(data = MSM_data,
                              formula = Y_4 ~ assigned_treatment_0 ,
                              weights = CalifuncmissGLkSL_4, family = 'gaussian')
      summary(PP_CalifuncmissGLkSL)
      
      PP_CalifuncmissAllSL <- glm(data = MSM_data,
                                  formula = Y_4 ~ assigned_treatment_0 ,
                                  weights = CalifuncmissAllSL_4, family = 'gaussian')
      summary(PP_CalifuncmissAllSL)
      

      ATE <- data.frame(k = c(4)) %>% 
        mutate(ATE_naive = predict.glm(PP_naive, data.frame(assigned_treatment_0 = 1))- predict.glm(PP_naive, data.frame(assigned_treatment_0 = 0)),
               ATE_GcompCorrect = GcompY1 - GcompY0,
               ATE_IPWcorrect = predict.glm(PP_IPWcorrect, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_IPWcorrect, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliLkcorrect = predict.glm(PP_CaliLkcorrect, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliLkcorrect, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliGLkcorrect = predict.glm(PP_CaliGLkcorrect, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliGLkcorrect, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliAllCorrect = predict.glm(PP_CaliAllCorrect, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliAllCorrect, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliGLkmiss = predict.glm(PP_CaliGLkmiss, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliGLkmiss, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliAllGLkmiss = predict.glm(PP_CaliAllGLkmiss, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliAllGLkmiss, data.frame(assigned_treatment_0 = 0)),
               ATE_Gcompmiss = GcompmissY1 - GcompmissY0,
               ATE_IPWmiss = predict.glm(PP_IPWmiss, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_IPWmiss, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliLkmiss = predict.glm(PP_CaliLkmiss, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliLkmiss, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliMissGLkcorrect = predict.glm(PP_CaliMissGLkcorrect, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliMissGLkcorrect, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliMissAllGLkcorrect = predict.glm(PP_CaliMissAllGLkcorrect, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliMissAllGLkcorrect, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliMissGLKmiss = predict.glm(PP_CaliMissGLKmiss, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliMissGLKmiss, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliAllmiss = predict.glm(PP_CaliAllmiss, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliAllmiss, data.frame(assigned_treatment_0 = 0)),
               ATE_Gcompfuncmiss = GcompfuncmissY1 - GcompfuncmissY0,
               ATE_IPWfuncmiss = predict.glm(PP_IPWfuncmiss, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_IPWfuncmiss, data.frame(assigned_treatment_0 = 0)),
               ATE_CalifuncmissLk = predict.glm(PP_CalifuncmissLk, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CalifuncmissLk, data.frame(assigned_treatment_0 = 0)),
               ATE_CalifuncmissGLk = predict.glm(PP_CalifuncmissGLk, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CalifuncmissGLk, data.frame(assigned_treatment_0 = 0)),
               ATE_CalifuncmissAll = predict.glm(PP_CalifuncmissAll, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CalifuncmissAll, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliCorrectGLkSL = predict.glm(PP_CaliCorrectGLkSL, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliCorrectGLkSL, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliAllCorrectSL = predict.glm(PP_CaliAllCorrectSL, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliAllCorrectSL, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliMissGLkSL = predict.glm(PP_CaliMissGLkSL, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliMissGLkSL, data.frame(assigned_treatment_0 = 0)),
               ATE_CaliMissAllSL = predict.glm(PP_CaliMissAllSL, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CaliMissAllSL, data.frame(assigned_treatment_0 = 0)),
               ATE_CalifuncmissGLkSL = predict.glm(PP_CalifuncmissGLkSL, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CalifuncmissGLkSL, data.frame(assigned_treatment_0 = 0)),
               ATE_CalifuncmissAllSL = predict.glm(PP_CalifuncmissAllSL, data.frame(assigned_treatment_0 = 1)) - predict.glm(PP_CalifuncmissAllSL, data.frame(assigned_treatment_0 = 0))
        )
      
      Y1 <- data.frame(k = c(4)) %>% 
        mutate(EY1_naive = predict.glm(PP_naive, data.frame(assigned_treatment_0 = 1)),
               EY1_GcompCorrect = GcompY1,
               EY1_IPWcorrect = predict.glm(PP_IPWcorrect, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliLkcorrect = predict.glm(PP_CaliLkcorrect, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliGLkcorrect = predict.glm(PP_CaliGLkcorrect, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliAllCorrect = predict.glm(PP_CaliAllCorrect, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliGLkmiss = predict.glm(PP_CaliGLkmiss, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliAllGLkmiss = predict.glm(PP_CaliAllGLkmiss, data.frame(assigned_treatment_0 = 1)),
               EY1_Gcompmiss = GcompmissY1,
               EY1_IPWmiss = predict.glm(PP_IPWmiss, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliLkmiss = predict.glm(PP_CaliLkmiss, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliMissGLkcorrect = predict.glm(PP_CaliMissGLkcorrect, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliMissAllGLkcorrect = predict.glm(PP_CaliMissAllGLkcorrect, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliMissGLKmiss = predict.glm(PP_CaliMissGLKmiss, data.frame(assigned_treatment_0 = 1)),
               EY1_CaliAllmiss = predict.glm(PP_CaliAllmiss, data.frame(assigned_treatment_0 = 1)),
               EY1_Gcompfuncmiss = GcompfuncmissY1,
               EY1_IPWfuncmiss = predict.glm(PP_IPWfuncmiss, data.frame(assigned_treatment_0 = 1)),
               EY1_CalifuncmissLk = predict.glm(PP_CalifuncmissLk, data.frame(assigned_treatment_0 = 1)),
               EY1_CalifuncmissGLk = predict.glm(PP_CalifuncmissGLk, data.frame(assigned_treatment_0 = 1)),
               EY1_CalifuncmissAll = predict.glm(PP_CalifuncmissAll, data.frame(assigned_treatment_0 = 1)), 
               EY1_CaliCorrectGLkSL = predict.glm(PP_CaliCorrectGLkSL, data.frame(assigned_treatment_0 = 1)) ,
               EY1_CaliAllCorrectSL = predict.glm(PP_CaliAllCorrectSL, data.frame(assigned_treatment_0 = 1)) ,
               EY1_CaliMissGLkSL = predict.glm(PP_CaliMissGLkSL, data.frame(assigned_treatment_0 = 1)) ,
               EY1_CaliMissAllSL = predict.glm(PP_CaliMissAllSL, data.frame(assigned_treatment_0 = 1)) ,
               EY1_CalifuncmissGLkSL = predict.glm(PP_CalifuncmissGLkSL, data.frame(assigned_treatment_0 = 1)) ,
               EY1_CalifuncmissAllSL = predict.glm(PP_CalifuncmissAllSL, data.frame(assigned_treatment_0 = 1)) 
        )
      Y0 <- data.frame(k = c(4)) %>% 
        mutate(EY0_naive = predict.glm(PP_naive, data.frame(assigned_treatment_0 = 0)),
               EY0_GcompCorrect = GcompY0,
               EY0_IPWcorrect = predict.glm(PP_IPWcorrect, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliLkcorrect = predict.glm(PP_CaliLkcorrect, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliGLkcorrect = predict.glm(PP_CaliGLkcorrect, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliAllCorrect = predict.glm(PP_CaliAllCorrect, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliGLkmiss = predict.glm(PP_CaliGLkmiss, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliAllGLkmiss = predict.glm(PP_CaliAllGLkmiss, data.frame(assigned_treatment_0 = 0)),
               EY0_Gcompmiss = GcompmissY0,
               EY0_IPWmiss = predict.glm(PP_IPWmiss, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliLkmiss = predict.glm(PP_CaliLkmiss, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliMissGLkcorrect = predict.glm(PP_CaliMissGLkcorrect, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliMissAllGLkcorrect = predict.glm(PP_CaliMissAllGLkcorrect, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliMissGLKmiss = predict.glm(PP_CaliMissGLKmiss, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliAllmiss = predict.glm(PP_CaliAllmiss, data.frame(assigned_treatment_0 = 0)),
               EY0_Gcompfuncmiss = GcompfuncmissY0,
               EY0_IPWfuncmiss = predict.glm(PP_IPWfuncmiss, data.frame(assigned_treatment_0 = 0)),
               EY0_CalifuncmissLk = predict.glm(PP_CalifuncmissLk, data.frame(assigned_treatment_0 = 0)),
               EY0_CalifuncmissGLk = predict.glm(PP_CalifuncmissGLk, data.frame(assigned_treatment_0 = 0)),
               EY0_CalifuncmissAll = predict.glm(PP_CalifuncmissAll, data.frame(assigned_treatment_0 = 0)),
               EY0_CaliCorrectGLkSL = predict.glm(PP_CaliCorrectGLkSL, data.frame(assigned_treatment_0 = 0)) ,
               EY0_CaliAllCorrectSL = predict.glm(PP_CaliAllCorrectSL, data.frame(assigned_treatment_0 = 0)) ,
               EY0_CaliMissGLkSL = predict.glm(PP_CaliMissGLkSL, data.frame(assigned_treatment_0 = 0)) ,
               EY0_CaliMissAllSL = predict.glm(PP_CaliMissAllSL, data.frame(assigned_treatment_0 = 0)) ,
               EY0_CalifuncmissGLkSL = predict.glm(PP_CalifuncmissGLkSL, data.frame(assigned_treatment_0 = 0)) ,
               EY0_CalifuncmissAllSL = predict.glm(PP_CalifuncmissAllSL, data.frame(assigned_treatment_0 = 0)) 
        )
      result$predict_estimates <- ATE
      result$EY1 <- Y1
      result$EY0 <- Y0
      
      return(result)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  save(oper, file = paste("Simulation results/result_simu_linear_outcome_cali_test_",as.character(l),".rda", sep = ""))
}

scenarios <- tidyr::crossing(size, treat,conf)

bias_ate <- array(,dim = c(26,7))
sd_ate <- array(,dim = c(26,7))
mae_ate <- array(,dim = c(26,7))
medae_ate <- array(,dim = c(26,7))
rootmse_ate <- array(,dim = c(26,7))

bias_EY1 <- array(,dim = c(26,7))
sd_EY1 <- array(,dim = c(26,7))
mae_EY1 <- array(,dim = c(26,7))
medae_EY1 <- array(,dim = c(26,7))
rootmse_EY1 <- array(,dim = c(26,7))

bias_EY0 <- array(,dim = c(26,7))
sd_EY0 <- array(,dim = c(26,7))
mae_EY0 <- array(,dim = c(26,7))
medae_EY0 <- array(,dim = c(26,7))
rootmse_EY0 <- array(,dim = c(26,7))

library(rlist)
library(matrixStats)
for (l in 1:6){
  load(paste0("Simulation results/result_simu_linear_outcome_cali_test_",as.character(l),".rda"))
  simu.t <- as.data.frame(1:iters)
  ate_all <- as.data.frame(list.rbind(oper[3,]))
  ate_all$simu <- simu.t[,1]
  ate_all <- ate_all[,-1]
  
  EY1_all <- as.data.frame(list.rbind(oper[1,]))
  EY1_all$simu <- simu.t[,1]
  EY1_all <- EY1_all[,-1]
  
  EY0_all <- as.data.frame(list.rbind(oper[2,]))
  EY0_all$simu <- simu.t[,1]
  EY0_all <- EY0_all[,-1]
  
  bias_ate[,1] <- colnames(ate_all[,-27])
  sd_ate[,1] <- colnames(ate_all[,-27])
  mae_ate[,1] <- colnames(ate_all[,-27])
  medae_ate[,1] <- colnames(ate_all[,-27])
  rootmse_ate[,1] <- colnames(ate_all[,-27])
  
  bias_ate[,l+1] <- colMeans(ate_all[,-27]) - true_ATE
  sd_ate[,l+1] <- colSds(as.matrix(ate_all[,-27]))
  mae_ate[,l+1] <- colMeans(abs(ate_all[,-27] - true_ATE))
  medae_ate[,l+1] <- colMedians(as.matrix(abs(ate_all[,-27] - true_ATE)))
  rootmse_ate[,l+1] <- sqrt(as.numeric(bias_ate[,l+1])^2 +as.numeric(sd_ate[,l+1])^2)
  
  bias_EY1[,1] <- colnames(EY1_all[,-27])
  sd_EY1[,1] <- colnames(EY1_all[,-27])
  mae_EY1[,1] <- colnames(EY1_all[,-27])
  medae_EY1[,1] <- colnames(EY1_all[,-27])
  rootmse_EY1[,1] <- colnames(EY1_all[,-27])
  
  bias_EY1[,l+1] <- colMeans(EY1_all[,-27]) - true_EY1
  sd_EY1[,l+1] <- colSds(as.matrix(EY1_all[,-27]))
  mae_EY1[,l+1] <- colMeans(abs(EY1_all[,-27] - true_EY1))
  medae_EY1[,l+1] <- colMedians(as.matrix(abs(EY1_all[,-27] - true_EY1)))
  rootmse_EY1[,l+1] <- sqrt(as.numeric(bias_EY1[,l+1])^2 +as.numeric(sd_EY1[,l+1])^2)
  
  bias_EY0[,1] <- colnames(EY0_all[,-27])
  sd_EY0[,1] <- colnames(EY0_all[,-27])
  mae_EY0[,1] <- colnames(EY0_all[,-27])
  medae_EY0[,1] <- colnames(EY0_all[,-27])
  rootmse_EY0[,1] <- colnames(EY0_all[,-27])
  
  bias_EY0[,l+1] <- colMeans(EY0_all[,-27]) - true_EY0
  sd_EY0[,l+1] <- colSds(as.matrix(EY0_all[,-27]))
  mae_EY0[,l+1] <- colMeans(abs(EY0_all[,-27] - true_EY0))
  medae_EY0[,l+1] <- colMedians(as.matrix(abs(EY0_all[,-27] - true_EY0)))
  rootmse_EY0[,l+1] <- sqrt(as.numeric(bias_EY0[,l+1])^2 +as.numeric(sd_EY0[,l+1])^2)
  
  
}

colnames(bias_ate) <- c('Method', 'N = 200', 'N = 500', 'N = 1000', 'N = 2500', 'N = 5000', 'N = 10000')
colnames(sd_ate) <- colnames(mae_ate) <- colnames(medae_ate) <-colnames(rootmse_ate) <-
  colnames(bias_EY1) <- colnames(sd_EY1) <- colnames(mae_EY1) <- colnames(medae_EY1) <-colnames(rootmse_EY1) <-
  colnames(bias_EY0) <- colnames(sd_EY0) <- colnames(mae_EY0) <- colnames(medae_EY0) <-colnames(rootmse_EY0) <-colnames(bias_ate)
  
bias_EY0 <- as.data.frame(bias_EY0)
bias_EY0[,-1] <- lapply(bias_EY0[,-1],as.numeric)
bias_EY0 <- bias_EY0 %>% 
  mutate_if(is.numeric, round,digits = 3)


bias_EY1 <- as.data.frame(bias_EY1)
bias_EY1[,-1] <- lapply(bias_EY1[,-1],as.numeric)
bias_EY1 <- bias_EY1 %>% 
  mutate_if(is.numeric, round,digits = 3)

rootmse_EY0 <- as.data.frame(rootmse_EY0)
rootmse_EY0[,-1] <- lapply(rootmse_EY0[,-1],as.numeric)
rootmse_EY0 <- rootmse_EY0 %>% 
  mutate_if(is.numeric, round,digits = 3)


rootmse_EY1 <- as.data.frame(rootmse_EY1)
rootmse_EY1[,-1] <- lapply(rootmse_EY1[,-1],as.numeric)
rootmse_EY1 <- rootmse_EY1 %>% 
  mutate_if(is.numeric, round,digits = 3)

mae_EY0 <- as.data.frame(mae_EY0)
mae_EY0[,-1] <- lapply(mae_EY0[,-1],as.numeric)
mae_EY0 <- mae_EY0 %>% 
  mutate_if(is.numeric, round,digits = 3)

mae_EY1 <- as.data.frame(mae_EY1)
mae_EY1[,-1] <- lapply(mae_EY1[,-1],as.numeric)
mae_EY1 <- mae_EY1 %>% 
  mutate_if(is.numeric, round,digits = 3)

medae_EY0 <- as.data.frame(medae_EY0)
medae_EY0[,-1] <- lapply(medae_EY0[,-1],as.numeric)
medae_EY0 <- medae_EY0 %>% 
  mutate_if(is.numeric, round,digits = 3)

medae_EY1 <- as.data.frame(medae_EY1)
medae_EY1[,-1] <- lapply(medae_EY1[,-1],as.numeric)
medae_EY1 <- medae_EY1 %>% 
  mutate_if(is.numeric, round,digits = 3)

