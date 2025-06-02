DATA_GEN<-function(ns, treat_prev, conf){  #ns=number of patients
  Z1<-rnorm(5*ns,0,1)
  Z2<-rnorm(5*ns,0,1)
  Z3<-rnorm(5*ns,0,1)
  Z4<-rnorm(5*ns,0,1)
  
  X1<-rep(0,5*ns)
  X2<-rep(0,5*ns)
  X3<-rep(0,5*ns)
  X4<-rep(0,5*ns)
  
  seq1<-seq(1,5*ns-4,5)
  
  Ap<-rep(0,5*ns)
  CAp<-rep(0,5*ns)
  
  X1[seq1]<-Z1[seq1]
  X2[seq1]<-Z2[seq1]
  X3[seq1]<-Z3[seq1]
  X4[seq1]<-Z4[seq1]
  
  P01<-1/(1+exp(-as.numeric(treat_prev)-0.5*X1[seq1]-0.5*X2[seq1]+0.2*X3[seq1]+0.2*X4[seq1]))
  
  A0<-rbinom(ns,1,P01)
  
  Y<-rep(0,5*ns)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp1<-A0+X1[seq1]+X2[seq1]+X3[seq1]+X4[seq1]
  
  Y[seq1]<-200+5*lp1+rnorm(ns,0,20)
  
  ### treatment assignment probability at 1st follow-up
  seq2<-seq(2,5*ns-3,5)
  
  Ap[seq2]<-A0
  CAp[seq2]<-A0
  
  U2 <- 1-0.3*A0
  
  X1[seq2]<-U2*Z1[seq2]
  X2[seq2]<-U2*Z2[seq2]
  X3[seq2]<-Z3[seq2]+as.numeric(conf)*CAp[seq2]
  X4[seq2]<-Z4[seq2]+as.numeric(conf)*CAp[seq2]
  
  P01<-1/(1+exp(-2.5*A0 + 2.5*(1-A0)-0.5*X1[seq2]-0.5*X2[seq2]+0.2*X3[seq2]+0.2*X4[seq2]))
  
  A1<-rbinom(ns,1,P01)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp2<-A1+X1[seq2]+X2[seq2]+X3[seq2]+X4[seq2]
  
  Y[seq2]<-200+5*lp2+rnorm(ns,0,20)
  
  ### treatment assignment probability at 2nd follow-up
  seq3<-seq(3,5*ns-2,5)
  Ap[seq3]<-A1
  CAp[seq3]<-CAp[seq2]+A1
  
  U3<-1-0.3*A1
  
  X1[seq3]<-U3*Z1[seq3]
  X2[seq3]<-U3*Z2[seq3]
  X3[seq3]<-Z3[seq3]+as.numeric(conf)*CAp[seq3]
  X4[seq3]<-Z4[seq3]+as.numeric(conf)*CAp[seq3]
  
  P02<-1/(1+exp(-2.5*A1 + 2.5*(1-A1)-0.5*X1[seq3]-0.5*X2[seq3]+0.2*X3[seq3]+0.2*X4[seq3]))
  
  A2<-rbinom(ns,1,P02)
  
  ### mean of the longitudinal outcome at 2nd visit given treatment history and covariate history
  
  lp3<-A2+X1[seq3]+X2[seq3]+X3[seq3]+X4[seq3]
  
  Y[seq3]<-200+5*lp3+rnorm(ns,0,20)
  
  ### treatment assignment probability at 3rd follow-up
  seq4<-seq(4,5*ns-1,5)
  Ap[seq4]<-A2
  CAp[seq4]<-CAp[seq3]+A2
  
  U4<-1-0.3*A2
  
  X1[seq4]<-U4*Z1[seq4]
  X2[seq4]<-U4*Z2[seq4]
  X3[seq4]<-Z3[seq4]+as.numeric(conf)*CAp[seq4]
  X4[seq4]<-Z4[seq4]+as.numeric(conf)*CAp[seq4]
  
  P03<-1/(1+exp(-2.5*A2 + 2.5*(1-A2)-0.5*X1[seq4]-0.5*X2[seq4]+0.2*X3[seq4]+0.2*X4[seq4]))
  
  A3<-rbinom(ns,1,P03)
  
 
  ### mean of the longitudinal outcome at 3rd visit given treatment history and covariate history
  
  lp4<-A3+X1[seq4]+X2[seq4]+X3[seq4]+X4[seq4]
  
  Y[seq4]<-200+5*lp4+rnorm(ns,0,20)
  
  ### treatment assignment probability at 4th follow-up
  seq5<-seq(5,5*ns,5)
  Ap[seq5]<-A3
  CAp[seq5]<-CAp[seq4]+A3
  
  U5<-1-0.3*A3
  
  X1[seq5]<-U5*Z1[seq5]
  X2[seq5]<-U5*Z2[seq5]
  X3[seq5]<-Z3[seq5]+as.numeric(conf)*CAp[seq5]
  X4[seq5]<-Z4[seq5]+as.numeric(conf)*CAp[seq5]
  
  P04<-1/(1+exp(-2.5*A3 + 2.5*(1-A3)-0.5*X1[seq5]-0.5*X2[seq5]+0.2*X3[seq5]+0.2*X4[seq5]))
  
  A4<-rbinom(ns,1,P04)
  
  
  
  ### mean of the longitudinal outcome at 4th visit given treatment history and covariate history
  
  lp5<-A4+X1[seq5]+X2[seq5]+X3[seq5]+X4[seq5]
  
  Y[seq5]<-200+5*lp5+rnorm(ns,0,20)
  
  
  
  ID<-rep(1:ns,each=5)
  
  A<-rep(0,5*ns)
  A[seq1]<-A0
  A[seq2]<-A1
  A[seq3]<-A2
  A[seq4]<-A3
  A[seq5]<-A4
  
  CA<-ave(A,ID,FUN=cumsum)
  
  ### transformed covariates
  TX1<-X1^3/9
  
  TX2<-X1*X2
  
  TX3<-log(abs(X3))+4
  
  TX4<-1/(1+exp(X4))
  
  DATA<-data.frame(ID,t=rep(c(0:4),ns),A,Ap,CA,X1,X2,X3,X4,TX1,TX2,TX3,TX4,Y)
  
  DATA
}
