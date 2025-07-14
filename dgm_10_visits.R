DATA_GEN_TEN<-function(ns, #ns=number of patients
                   conf = 0.2, 
                   treat_prev_0 = 0, 
                   treat_prev_d1 = c(1.85,1.65,1.45,1.25,1.05,0.85,0.65,0.45,0.25), 
                   treat_prev_d0 = c(-2.15, -2.15, -2.15, -2.15, -2.15, -2.15, -2.15, -2.15, -2.15)){
  Z1<-rnorm(10*ns,0,1)
  Z2<-rnorm(10*ns,0,1)
  Z3<-rnorm(10*ns,0,1)
  Z4<-rnorm(10*ns,0,1)
  
  X1<-rep(0,10*ns)
  X2<-rep(0,10*ns)
  X3<-rep(0,10*ns)
  X4<-rep(0,10*ns)
  
  #Baseline
  seq1<-seq(1,10*ns - (10-1) ,10)
  
  Ap<-rep(0,10*ns)
  CAp<-rep(0,10*ns)
  A <- rep(0,10*ns)
  X1[seq1]<-Z1[seq1]
  X2[seq1]<-Z2[seq1]
  X3[seq1]<-Z3[seq1]
  X4[seq1]<-Z4[seq1]
  
  P01<-1/(1+exp(-treat_prev_0-0.5*X1[seq1]-0.5*X2[seq1]-as.numeric(conf)*X3[seq1]-as.numeric(conf)*X4[seq1]))
  
  A[seq1]<-rbinom(ns,1,P01)
  
  Y<-rep(0,10*ns)
  
  ### mean of the longitudinal outcome at baseline given treatment history and covariate history
  lp1<-2*A[seq1]+X1[seq1]+X2[seq1]+X3[seq1]+X4[seq1]
  
  Y[seq1]<-200+5*lp1+rnorm(ns,0,20)
  
  ## update covariates at t = 1
  seq2<-seq1+1
  
  Ap[seq2]<-A[seq1]
  CAp[seq2]<- CAp[seq1]+ A[seq1]
  
  U2 <- 1-0.3*Ap[seq2]
  
  X1[seq2]<-U2*Z1[seq2]
  X2[seq2]<-U2*Z2[seq2]
  X3[seq2]<-Z3[seq2]+0.5*CAp[seq2]
  X4[seq2]<-Z4[seq2]+0.5*CAp[seq2]
  
  P1<-1/(1+exp(-as.numeric(treat_prev_d1[1])*Ap[seq2] -as.numeric(treat_prev_d0[1])*(1-Ap[seq2])-0.5*X1[seq2]-0.5*X2[seq2]-as.numeric(conf)*X3[seq2]-as.numeric(conf)*X4[seq2]))
  
  A[seq2]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq2]+Ap[seq2]+X1[seq2]+X2[seq2]+X3[seq2]+X4[seq2] + X1[seq1]+X2[seq1]+X3[seq1]+X4[seq1]
  
  Y[seq2]<-200+5*lp+rnorm(ns,0,20)
  
  ## update covariates at t = 2
  seq3<-seq2+1
  
  Ap[seq3]<-A[seq2]
  CAp[seq3]<- CAp[seq2]+ A[seq2]
  
  U2 <- 1-0.3*Ap[seq3]
  
  X1[seq3]<-U2*Z1[seq3]
  X2[seq3]<-U2*Z2[seq3]
  X3[seq3]<-Z3[seq3]+0.5*CAp[seq3]
  X4[seq3]<-Z4[seq3]+0.5*CAp[seq3]
  
  P1<-1/(1+exp(-as.numeric(treat_prev_d1[2])*Ap[seq3] -as.numeric(treat_prev_d0[2])*(1-Ap[seq3])-0.5*X1[seq3]-0.5*X2[seq3]-as.numeric(conf)*X3[seq3]-as.numeric(conf)*X4[seq3]))
  
  A[seq3]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq3]+Ap[seq3]+X1[seq3]+X2[seq3]+X3[seq3]+X4[seq3] + X1[seq2]+X2[seq2]+X3[seq2]+X4[seq2]
  
  Y[seq3]<-200+5*lp+rnorm(ns,0,20)
  
  
  ## update covariates at t = 3
  seq4<-seq3+1
  
  Ap[seq4]<-A[seq3]
  CAp[seq4]<- CAp[seq3]+ A[seq3]
  
  U2 <- 1-0.3*Ap[seq4]
  
  X1[seq4]<-U2*Z1[seq4]
  X2[seq4]<-U2*Z2[seq4]
  X3[seq4]<-Z3[seq4]+0.5*CAp[seq4]
  X4[seq4]<-Z4[seq4]+0.5*CAp[seq4]
  
  P1<-1/(1+exp(-as.numeric(treat_prev_d1[3])*Ap[seq4] -as.numeric(treat_prev_d0[3])*(1-Ap[seq4])-0.5*X1[seq4]-0.5*X2[seq4]-as.numeric(conf)*X3[seq4]-as.numeric(conf)*X4[seq4]))
  
  A[seq4]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq4]+Ap[seq4]+X1[seq4]+X2[seq4]+X3[seq4]+X4[seq4] + X1[seq3]+X2[seq3]+X3[seq3]+X4[seq3]
  
  Y[seq4]<-200+5*lp+rnorm(ns,0,20)
  
  
  ## update covariates at t = 4
  seq5<-seq4+1
  
  Ap[seq5]<-A[seq4]
  CAp[seq5]<- CAp[seq4]+ A[seq4]
  
  U2 <- 1-0.3*Ap[seq5]
  
  X1[seq5]<-U2*Z1[seq5]
  X2[seq5]<-U2*Z2[seq5]
  X3[seq5]<-Z3[seq5]+0.5*CAp[seq5]
  X4[seq5]<-Z4[seq5]+0.5*CAp[seq5]
  
  P1<-1/(1+exp(-as.numeric(treat_prev_d1[4])*Ap[seq5] -as.numeric(treat_prev_d0[4])*(1-Ap[seq5])-0.5*X1[seq5]-0.5*X2[seq5]-as.numeric(conf)*X3[seq5]-as.numeric(conf)*X4[seq5]))
  
  A[seq5]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq5]+Ap[seq5]+X1[seq5]+X2[seq5]+X3[seq5]+X4[seq5] + X1[seq4]+X2[seq4]+X3[seq4]+X4[seq4]
  
  Y[seq5]<-200+5*lp+rnorm(ns,0,20)
  
  ## update covariates at t = 5
  seq6<-seq5+1
  
  Ap[seq6]<-A[seq5]
  CAp[seq6]<- CAp[seq5]+ A[seq5]
  
  U2 <- 1-0.3*Ap[seq6]
  
  X1[seq6]<-U2*Z1[seq6]
  X2[seq6]<-U2*Z2[seq6]
  X3[seq6]<-Z3[seq6]+0.5*CAp[seq6]
  X4[seq6]<-Z4[seq6]+0.5*CAp[seq6]
  
  P1<-1/(1+exp(-as.numeric(treat_prev_d1[5])*Ap[seq6] -as.numeric(treat_prev_d0[5])*(1-Ap[seq6])-0.5*X1[seq6]-0.5*X2[seq6]-as.numeric(conf)*X3[seq6]-as.numeric(conf)*X4[seq6]))
  
  A[seq6]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq6]+Ap[seq6]+X1[seq6]+X2[seq6]+X3[seq6]+X4[seq6] + X1[seq5]+X2[seq5]+X3[seq5]+X4[seq5]
  
  Y[seq6]<-200+5*lp+rnorm(ns,0,20)
  
############### update covariates at t = 6 ##############
  seq7<-seq6+1
  
  Ap[seq7]<-A[seq6]
  CAp[seq7]<- CAp[seq6]+ A[seq6]
  
  U2 <- 1-0.3*Ap[seq7]
  
  X1[seq7]<-U2*Z1[seq7]
  X2[seq7]<-U2*Z2[seq7]
  X3[seq7]<-Z3[seq7]+0.5*CAp[seq7]
  X4[seq7]<-Z4[seq7]+0.5*CAp[seq7]
  
  P1<-1/(1+exp(-as.numeric(treat_prev_d1[6])*Ap[seq7] -as.numeric(treat_prev_d0[6])*(1-Ap[seq7])-0.5*X1[seq7]-0.5*X2[seq7]-as.numeric(conf)*X3[seq7]-as.numeric(conf)*X4[seq7]))
  
  A[seq7]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq7]+Ap[seq7]+X1[seq7]+X2[seq7]+X3[seq7]+X4[seq7] + X1[seq6]+X2[seq6]+X3[seq6]+X4[seq6]
  
  Y[seq7]<-200+5*lp+rnorm(ns,0,20)

  ############### update covariates at t = 7 ##############
  seq8<-seq7+1
  
  Ap[seq8]<-A[seq7]
  CAp[seq8]<- CAp[seq7]+ A[seq7]
  
  U2 <- 1-0.3*Ap[seq8]
  
  X1[seq8]<-U2*Z1[seq8]
  X2[seq8]<-U2*Z2[seq8]
  X3[seq8]<-Z3[seq8]+0.5*CAp[seq8]
  X4[seq8]<-Z4[seq8]+0.5*CAp[seq8]
  
  P1<-1/(1+exp(-as.numeric(treat_prev_d1[7])*Ap[seq8] -as.numeric(treat_prev_d0[7])*(1-Ap[seq8])-0.5*X1[seq8]-0.5*X2[seq8]-as.numeric(conf)*X3[seq8]-as.numeric(conf)*X4[seq8]))
  
  A[seq8]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq8]+Ap[seq8]+X1[seq8]+X2[seq8]+X3[seq8]+X4[seq8] + X1[seq7]+X2[seq7]+X3[seq7]+X4[seq7]
  
  Y[seq8]<-200+5*lp+rnorm(ns,0,20) 
  
  ############### update covariates at t = 8 ##############
  seq9<-seq8+1
  
  Ap[seq9]<-A[seq8]
  CAp[seq9]<- CAp[seq8]+ A[seq8]
  
  U2 <- 1-0.3*Ap[seq9]
  
  X1[seq9]<-U2*Z1[seq9]
  X2[seq9]<-U2*Z2[seq9]
  X3[seq9]<-Z3[seq9]+0.5*CAp[seq9]
  X4[seq9]<-Z4[seq9]+0.5*CAp[seq9]
  
  P1<-1/(1+exp(-as.numeric(treat_prev_d1[8])*Ap[seq9] -as.numeric(treat_prev_d0[8])*(1-Ap[seq9])-0.5*X1[seq9]-0.5*X2[seq9]-as.numeric(conf)*X3[seq9]-as.numeric(conf)*X4[seq9]))
  
  A[seq9]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq9]+Ap[seq9]+X1[seq9]+X2[seq9]+X3[seq9]+X4[seq9] + X1[seq8]+X2[seq8]+X3[seq8]+X4[seq8]
  
  Y[seq9]<-200+5*lp+rnorm(ns,0,20) 
  
  ############### update covariates at t = 9 ##############
  seq10<-seq9+1
  
  Ap[seq10]<-A[seq9]
  CAp[seq10]<- CAp[seq9]+ A[seq9]
  
  U2 <- 1-0.3*Ap[seq10]
  
  X1[seq10]<-U2*Z1[seq10]
  X2[seq10]<-U2*Z2[seq10]
  X3[seq10]<-Z3[seq10]+0.5*CAp[seq10]
  X4[seq10]<-Z4[seq10]+0.5*CAp[seq10]
  
  P1<-1/(1+exp(-as.numeric(treat_prev_d1[9])*Ap[seq10] -as.numeric(treat_prev_d0[9])*(1-Ap[seq10])-0.5*X1[seq10]-0.5*X2[seq10]-as.numeric(conf)*X3[seq10]-as.numeric(conf)*X4[seq10]))
  
  A[seq10]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq10]+Ap[seq10]+X1[seq10]+X2[seq10]+X3[seq10]+X4[seq10] + X1[seq9]+X2[seq9]+X3[seq9]+X4[seq9]
  
  Y[seq10]<-200+5*lp+rnorm(ns,0,20) 

  ID<-rep(1:ns,each=10)
  
  CA<-ave(A,ID,FUN=cumsum)
  
  ### transformed covariates
  TX1<-X1^3/9
  
  TX2<-X1*X2
  
  TX3<-log(abs(X3))+4
  
  TX4<-1/(1+exp(X4))
  
  DATA<-data.frame(ID,t=rep(0:9,ns),A,Ap,CA,X1,X2,X3,X4,TX1,TX2,TX3,TX4,Y)
  
  DATA
}
