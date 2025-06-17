DATA_GEN<-function(ns, nv, conf = 0.5){  #ns=number of patients
  Z1<-rnorm(nv*ns,0,1)
  Z2<-rnorm(nv*ns,0,1)
  Z3<-rnorm(nv*ns,0,1)
  Z4<-rnorm(nv*ns,0,1)
  
  X1<-rep(0,nv*ns)
  X2<-rep(0,nv*ns)
  X3<-rep(0,nv*ns)
  X4<-rep(0,nv*ns)
  
  seq1<-seq(1,nv*ns - (nv-1) ,nv)
  
  Ap<-rep(0,nv*ns)
  CAp<-rep(0,nv*ns)
  A <- rep(0,nv*ns)
  X1[seq1]<-Z1[seq1]
  X2[seq1]<-Z2[seq1]
  X3[seq1]<-Z3[seq1]
  X4[seq1]<-Z4[seq1]
  
  P01<-1/(1+exp(-0.5*X1[seq1]-0.5*X2[seq1]+as.numeric(conf)*X3[seq1]+as.numeric(conf)*X4[seq1]))
  
  A[seq1]<-rbinom(ns,1,P01)
  
  Y<-rep(0,nv*ns)
  
  ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
  lp1<-2*A[seq1]+X1[seq1]+X2[seq1]+X3[seq1]+X4[seq1]
  
  Y[seq1]<-200+5*lp1+rnorm(ns,0,20)
  
  seqlist<-list()                              
  seqlist[[1]]<-seq1
  
  for (k in 2:nv){  
    ## update covariates
    seqlist[[k]]<-seqlist[[k-1]]+1
  
    Ap[seqlist[[k]]]<-A[seqlist[[k-1]]]
    CAp[seqlist[[k]]]<- CAp[seqlist[[k-1]]]+ A[seqlist[[k-1]]]
    
    U2 <- 1-0.3*Ap[seqlist[[k]]]
    
    X1[seqlist[[k]]]<-U2*Z1[seqlist[[k]]]
    X2[seqlist[[k]]]<-U2*Z2[seqlist[[k]]]
    X3[seqlist[[k]]]<-Z3[seqlist[[k]]]+0.5*CAp[seqlist[[k]]]
    X4[seqlist[[k]]]<-Z4[seqlist[[k]]]+0.5*CAp[seqlist[[k]]]
    
    P1<-1/(1+exp(-2.5*Ap[seqlist[[k]]] + 2.5*(1-Ap[seqlist[[k]]])-0.5*X1[seqlist[[k]]]-0.5*X2[seqlist[[k]]]+as.numeric(conf)*X3[seqlist[[k]]]+as.numeric(conf)*X4[seqlist[[k]]]))
    
    A[seqlist[[k]]]<-rbinom(ns,1,P1)
    
    ### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
    lp<-2*A[seqlist[[k]]]+Ap[seqlist[[k]]]+X1[seqlist[[k]]]+X2[seqlist[[k]]]+X3[seqlist[[k]]]+X4[seqlist[[k]]] + X1[seqlist[[k-1]]]+X2[seqlist[[k-1]]]+X3[seqlist[[k-1]]]+X4[seqlist[[k-1]]]
    
    Y[seqlist[[k]]]<-200+5*lp+rnorm(ns,0,20)
  }
  
  ID<-rep(1:ns,each=nv)
  
  CA<-ave(A,ID,FUN=cumsum)
  
  ### transformed covariates
  TX1<-X1^3/9
  
  TX2<-X1*X2
  
  TX3<-log(abs(X3))+4
  
  TX4<-1/(1+exp(X4))
  
  DATA<-data.frame(ID,t=rep(c(0:(nv-1)),ns),A,Ap,CA,X1,X2,X3,X4,TX1,TX2,TX3,TX4,Y)
  
  DATA
}
