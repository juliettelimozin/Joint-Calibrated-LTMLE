
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

simdata_alltreat <- DATA_GEN_cont_outcome_treatment_switch_history(ns = as.numeric(size[l]),nv = 5,treat_prev =  0.5,
                                                                   conf =  1.5,
                                                                   censor = F, all_treat = TRUE)
simdata_allcontrol <- DATA_GEN_cont_outcome_treatment_switch_history(ns = as.numeric(size[l]),nv = 5,treat_prev =  0.5,
                                                                     conf =  1.5,
                                                                     censor = F, all_control = TRUE)
mean(simdata_alltreat[simdata_alltreat$t == 4,]$Y)
mean(simdata_allcontrol[simdata_allcontrol$t == 4,]$Y)
min(simdata$Y)

############################### SIMULATION ####################################

#Simulation setp up:
# Correct specification: treatment or outcome model includes an interaction term
# Misspecification: No interaction term
# Functional misspecification: transformed covariates U1 = log(|X1|)/3, U2 = sqrt(|X2|)/4






