DATA_GEN_censored<-function(ns, #ns=number of patients
                   conf = 0.2, 
                   treat_prev_0 = 0, 
                   treat_prev_d1_1 = 1, 
                   treat_prev_d0_1 = -1.25, 
                   treat_prev_d1_2 =0.8, 
                   treat_prev_d0_2 = -1.25){
  Z1<-rnorm(3*ns,0,1)
  Z2<-rnorm(3*ns,0,1)
  Z3<-rnorm(3*ns,0,1)
  Z4<-rnorm(3*ns,0,1)
  
  X1<-rep(0,3*ns)
  X2<-rep(0,3*ns)
  X3<-rep(0,3*ns)
  X4<-rep(0,3*ns)
  
  #Baseline
  seq1<-seq(1,3*ns - (3-1) ,3)
  
  Ap<-rep(0,3*ns)
  CAp<-rep(0,3*ns)
  A <- rep(0,3*ns)
  
  X1[seq1]<-Z1[seq1]
  X2[seq1]<-Z2[seq1]
  X3[seq1]<-Z3[seq1]
  X4[seq1]<-Z4[seq1]
  
  P01<-1/(1+exp(-treat_prev_0-0.5*X1[seq1]-0.5*X2[seq1]-as.numeric(conf)*X3[seq1]-as.numeric(conf)*X4[seq1]))
  
  A[seq1]<-rbinom(ns,1,P01)
  
  Y<-rep(0,3*ns)

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
  
  P1<-1/(1+exp(-treat_prev_d1_1*Ap[seq2] -treat_prev_d0_1*(1-Ap[seq2])-0.5*X1[seq2]-0.5*X2[seq2]-as.numeric(conf)*X3[seq2]-as.numeric(conf)*X4[seq2]))
  
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
  
  P1<-1/(1+exp(-treat_prev_d1_2*Ap[seq3] -treat_prev_d0_2*(1-Ap[seq3])-0.5*X1[seq3]-0.5*X2[seq3]-as.numeric(conf)*X3[seq3]-as.numeric(conf)*X4[seq3]))
  
  A[seq3]<-rbinom(ns,1,P1)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp<-2*A[seq3]+Ap[seq3]+X1[seq3]+X2[seq3]+X3[seq3]+X4[seq3] + X1[seq2]+X2[seq2]+X3[seq2]+X4[seq2]
  
  Y[seq3]<-200+5*lp+rnorm(ns,0,20)
  
  ID<-rep(1:ns,each=3)
  
  CA<-ave(A,ID,FUN=cumsum)
  
  ### Dropout
  C<- rep(0,3*ns)
  Dprob<-1/(1+exp(2.5 -0.5*X1[seq1] +0.5*X2[seq1] -0.2*X3[seq1] +0.2*X4[seq1])) ##Probability of dropout
  
  C[seq1]<-rbinom(ns,1,Dprob)
  
  Dprob<-1/(1+exp(2.5 -0.5*X1[seq2] +0.5*X2[seq2] -0.2*X3[seq2] +0.2*X4[seq2])) ##Probability of dropout
  
  C[seq2]<-rbinom(ns,1,Dprob)
  
  indfun<-function(n){
    if (sum(n)==0) {rep(0,3)}
    else{k<-min(which(n==1))
    c(rep(0,k),rep(1,3-k))}}
  
  RL<-ave(C,ID,FUN=indfun)
  ### transformed covariates
  TX1<-X1^3/9
  
  TX2<-X1*X2
  
  TX3<-log(abs(X3))+4
  
  TX4<-1/(1+exp(X4))
  
  DATA<-data.frame(ID,t=rep(0:2,ns),A,Ap,CA,X1,X2,X3,X4,TX1,TX2,TX3,TX4,Y,C)
  DATA[RL ==0,]
}
