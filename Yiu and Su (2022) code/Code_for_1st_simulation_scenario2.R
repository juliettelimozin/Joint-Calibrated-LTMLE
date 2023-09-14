##Data generating function

DATA_GEN<-function(ns){  ##Generate Z
Z1<-rnorm(11*ns,0,1)
Z2<-rnorm(11*ns,0,1)
Z3<-rnorm(11*ns,0,1)
Z4<-rnorm(11*ns,0,1)

A00<-rbinom(ns,1,0.5)     ##Generate initial treatment at random
A10<-rbinom(ns,1,0.5)*A00  

X1<-rep(0,11*ns)          ##Covariate vectors
X2<-rep(0,11*ns)
X3<-rep(0,11*ns)
X4<-rep(0,11*ns)

A0P<-rep(0,11*ns)        ##Previous and current treatment
A1P<-rep(0,11*ns)
A0<-rep(0,11*ns)
A1<-rep(0,11*ns)
CA0P<-rep(0,11*ns)   ##vector for sum of previous treatment A0, these are used to affect current covariates
CA1P<-rep(0,11*ns)   ##vector for sum of previous treatment A1, these are used to affect current covariates

Y<-rep(0,11*ns)     ##Outcome vector

                    ##Fill in initial values
seq1<-seq(1,11*ns-10,11)
A0P[seq1]<-A00
A1P[seq1]<-A10
CA0P[seq1]<-A00
CA1P[seq1]<-A10

U<-rep(0,11*ns)     ##Define U in the simulation set-up
U[seq1]<-1-0.3*A00+0.3*A10

X1[seq1]<-U[seq1]*Z1[seq1]
X2[seq1]<-U[seq1]*Z2[seq1]
X3[seq1]<-Z3[seq1]+0.5*A00+0.5*A10
X4[seq1]<-Z4[seq1]+0.5*A00+0.5*A10

P0<-list()                                   ##list of treatment probabilities 
P1<-list()
seqlist<-list()                              
seqlist[[1]]<-seq1

for (k in 2:11){                             ##Generate Data

P0[[k-1]]<-1/(1+exp(-A0P[seqlist[[k-1]]]-A1P[seqlist[[k-1]]]-0.5*X1[seqlist[[k-1]]]-0.5*X2[seqlist[[k-1]]]+0.2*X3[seqlist[[k-1]]]+
0.2*X4[seqlist[[k-1]]]))

A0[seqlist[[k-1]]]<-rbinom(ns,1,P0[[k-1]]) ##Generate treatment at current visit based on previous covariates

P1[[k-1]]<-1/(1+exp(-A0P[seqlist[[k-1]]]-A1P[seqlist[[k-1]]]-0.5*X1[seqlist[[k-1]]]-0.5*X2[seqlist[[k-1]]]+0.2*X3[seqlist[[k-1]]]+
0.2*X4[seqlist[[k-1]]]))

A1[seqlist[[k-1]]]<-rbinom(ns,1,P1[[k-1]])*A0[seqlist[[k-1]]]

##Update values

seqlist[[k]]<-seqlist[[k-1]]+1

A0P[seqlist[[k]]]<-A0[seqlist[[k-1]]]
A1P[seqlist[[k]]]<-A1[seqlist[[k-1]]]
CA0P[seqlist[[k]]]<-CA0P[seqlist[[k-1]]]+A0P[seqlist[[k]]]
CA1P[seqlist[[k]]]<-CA1P[seqlist[[k-1]]]+A1P[seqlist[[k]]]

U[seqlist[[k]]]<-1-0.3*A0P[seqlist[[k]]]-0.3*A1P[seqlist[[k]]] 

X1[seqlist[[k]]]<-U[seqlist[[k]]]*Z1[seqlist[[k]]]
X2[seqlist[[k]]]<-U[seqlist[[k]]]*Z2[seqlist[[k]]]
X3[seqlist[[k]]]<-Z3[seqlist[[k]]]+0.5*CA0P[seqlist[[k]]]+0.5*CA1P[seqlist[[k]]]
X4[seqlist[[k]]]<-Z4[seqlist[[k]]]+0.5*CA0P[seqlist[[k]]]+0.5*CA1P[seqlist[[k]]]

lp<-A0[seqlist[[k-1]]]+A1[seqlist[[k-1]]]+X1[seqlist[[k]]]+X2[seqlist[[k]]]+X3[seqlist[[k]]]+X4[seqlist[[k]]]+
X1[seqlist[[k-1]]]+X2[seqlist[[k-1]]]+X3[seqlist[[k-1]]]+X4[seqlist[[k-1]]]

##Genrate outcome

Y[seqlist[[k]]]<-200+5*lp+rnorm(ns,0,20)}

##Make data frame

ID<-rep(1:ns,each=10)

##Align data by removing values 
NSEQ<-seq(11,11*ns,11)

X1<-X1[-NSEQ]
X2<-X2[-NSEQ]
X3<-X3[-NSEQ]
X4<-X4[-NSEQ]
A0P<-A0P[-NSEQ]
A1P<-A1P[-NSEQ]
A0<-A0[-NSEQ]
A1<-A1[-NSEQ]

Y<-Y[-seq(1,11*ns-10,11)]

##Cumulative treatment effect

CA0<-ave(A0,ID,FUN=cumsum)+rep(A00,each=10)

CA1<-ave(A1,ID,FUN=cumsum)+rep(A10,each=10)

##Define transform covariates

TX1<-X1^3/9

TX2<-X1*X2

TX3<-log(abs(X3))+4

TX4<-1/(1+exp(X4))

##Dropout

Dprob<-1/(1+exp(1+A0P+A1P+0.5*X1+0.5*X2+0.2*X3+0.2*X4)) ##Probability of dropout

R<-rbinom(10*ns,1,Dprob) ##R=0 is remain in the study

##Create Data frame

DATA<-data.frame(ID,T=rep(c(1:10),ns),A0,A1,A0P,A1P,CA0,CA1,X1,X2,X3,X4,TX1,TX2,TX3,TX4,Y,R)

##Remove observations after dropout

indfun<-function(n){
if (sum(n)==0) {rep(0,10)}
else{k<-min(which(n==1))
c(rep(0,k),rep(1,10-k))}}

RL<-ave(R,DATA$ID,FUN=indfun)

indR<-which(RL==0)

DATA1<-DATA[indR,]
DATA1}

######

library(compiler)
enableJIT(3)

library(doParallel)
registerDoParallel(2) 
library(foreach)


##################################
####################

RE<-foreach (i=1:2500,.packages=c('compiler','nleqslv')) %dopar% {
    enableJIT(3)
    
    DATADROP<-DATA_GEN(1000) ##Generate data including dropout visit
    
    indR0<-which(DATADROP$R==0)
    DATA<-DATADROP[indR0,]   ##Complete data not including dropout visit
    
    ##MLE method
    
    ##Treatment weights
    
    Proj_A0<-glm(A0~A0P+A1P,data=DATA,binomial)                          ##fit numerator model for A0   
    
    Prob_A0<-glm(A0~A0P+A1P+X1+X2+X3+X4,data=DATA,binomial)              ##fit denominator model for A0 with correct covariates
    Prob_A0T<-glm(A0~A0P+A1P+TX1+TX2+TX3+TX4,data=DATA,binomial)         ##fit denominator model for A0 with transformed covariates
    
    weProj_A0<-fitted(Proj_A0)*DATA$A0+(1-DATA$A0)*(1-fitted(Proj_A0))      ##numerator of weight for A0
    
    weProb_A0<-fitted(Prob_A0)*DATA$A0+(1-DATA$A0)*(1-fitted(Prob_A0))      ##denominator of weight for A0 with correct covariates
    weProb_A0T<-fitted(Prob_A0T)*DATA$A0+(1-DATA$A0)*(1-fitted(Prob_A0T))   ##denominator of weight for A0 with transformed covariates
    
    indA1<-which(DATA$A0==1)
    
    Proj_A1<-glm(A1~A0P+A1P,data=DATA[indA1,],binomial)                    ##fit numerator model for A1
    
    Prob_A1<-glm(A1~A0P+A1P+X1+X2+X3+X4,data=DATA[indA1,],binomial)           ##fit denominator model for A1 with correct covariates
    Prob_A1T<-glm(A1~A0P+A1P+TX1+TX2+TX3+TX4,data=DATA[indA1,],binomial)      ##fit denominator model for A1 with transformed covariates
    
    weProj_A1<-rep(1,nrow(DATA))
    
    weProj_A1[indA1]<-fitted(Proj_A1)*DATA$A1[indA1]+(1-DATA$A1[indA1])*(1-fitted(Proj_A1))  ##numerator of weight for A1
    
    weProb_A1<-rep(1,nrow(DATA))
    weProb_A1T<-rep(1,nrow(DATA))
    
    weProb_A1[indA1]<-fitted(Prob_A1)*DATA$A1[indA1]+(1-DATA$A1[indA1])*(1-fitted(Prob_A1))     ##denominator of weight for A1 with correct covariates
    weProb_A1T[indA1]<-fitted(Prob_A1T)*DATA$A1[indA1]+(1-DATA$A1[indA1])*(1-fitted(Prob_A1T))  ##denominator of weight for A1 with transformed covariates
    
    we_Prob_AN<-weProj_A0*weProj_A1/(weProb_A0*weProb_A1)        ##MLE weights with correct covariates                         
    we_Prob_ANT<-weProj_A0*weProj_A1/(weProb_A0T*weProb_A1T)     ##MLE weights with transformed covariates
    
    weightsMLEA<-ave(we_Prob_AN,DATA$ID,FUN=cumprod)             ##MLE of cumulative weights with correct covariates
    weightsMLEAT<-ave(we_Prob_ANT,DATA$ID,FUN=cumprod)           ##MLE of cumulative weights with transformed covariates
    
    WeightsMLEA<-weightsMLEA/mean(weightsMLEA)                   ## Scaled weights with correct covariates
    WeightsMLEAT<-weightsMLEAT/mean(weightsMLEAT)                ## Scaled weights with transformed covariates
    
    ##Dropout weights
    
    Prob_R<-glm((1-R)~A0P+A1P+X1+X2+X3+X4,data=DATADROP,binomial)       ##fit denominator model for 1-R (remain in the study) with correct covariates
    Prob_RT<-glm((1-R)~A0P+A1P+TX1+TX2+TX3+TX4,data=DATADROP,binomial)  ##fit denominator model for 1-R (remain in the study) with transformed covariates
    
    we_R<-1/fitted(Prob_R)       ##MLE weights for dropout with correct covariates
    we_RT<-1/fitted(Prob_RT)     ##MLE weights for dropout with transformed covariates
    
    we_Prob_RN<-we_R[indR0]
    we_Prob_RNT<-we_RT[indR0]
    
    weightsMLER<-ave(we_Prob_RN,DATA$ID,FUN=cumprod)   ##cumulative treatment and dropout weights with correct covariates
    weightsMLERT<-ave(we_Prob_RNT,DATA$ID,FUN=cumprod) ##cumulative treatment and dropout weights with transformed covariates
    
    Npop<-length(which(DATADROP$T==1)) ##Number of subjects in the study at visit 1
    
    WeightsMLEAR<-10*Npop*weightsMLER*weightsMLEA/sum(weightsMLER*weightsMLEA)      ##Scaled weights with correct covariates
    WeightsMLEART<-10*Npop*weightsMLERT*weightsMLEAT/sum(weightsMLERT*weightsMLEAT) ##Scaled weights with transformed covariates
    
    COEF_MLE<-lm(Y~CA0+CA1,data=DATA,weights=WeightsMLEAR)$coef[2:3]       ##Treatment effect with correct covariates
    COEF_MLET<-lm(Y~CA0+CA1,data=DATA,weights=WeightsMLEART)$coef[2:3]     ##Treatment effect with transformed covariates
     
    
    ###############################Calibration
    
    library(nleqslv)
    
    ##construct restrictions (7) in the main text for the correct covariates
    ##Construct restrictions for A0
    
    A0pi<-fitted(Proj_A0)
    
    A0kernal<-DATA$A0-A0pi
    
    DMAT<-cbind(rep(1,nrow(DATA)),DATA$X1,DATA$X2,DATA$X3,DATA$X4,DATA$A0P,DATA$A1P)
    
    DMATA0<-A0kernal*DMAT
    
    avecumsum<-function(n){ave(n,DATA$ID,FUN=cumsum)}
    
    DMATA0A<-apply(DMATA0,2,avecumsum) ##DMATA0A is a matrix containing the coefficients of the weights in the restrictions for A0
    
    
    ##Construct restrictions for A1
    
    A1pi<-1/(1+exp(-Proj_A1$coef[1]-Proj_A1$coef[2]*DATA$A0P-Proj_A1$coef[3]*DATA$A1P))
    
    A1kernal<-DATA$A0*(DATA$A1-A1pi)
    
    DMATA1<-A1kernal*DMAT
    
    DMATA1A<-apply(DMATA1,2,avecumsum) ##DMATA1A is a matrix containing the coefficients of the weights in the restrictions for A1
    
    ##Overall constraints
    
    DMATAO<-cbind(DMATA0A,DMATA1A)
    
    ##Construct restrictions (7) in the main text for the transformed covariates
    
    DMATT<-cbind(rep(1,nrow(DATA)),DATA$TX1,DATA$TX2,DATA$TX3,DATA$TX4,DATA$A0P,DATA$A1P)
    
    DMATA0T<-A0kernal*DMATT
    
    DMATA0AT<-apply(DMATA0T,2,avecumsum)
    
    ##
    
    DMATA1T<-A1kernal*DMATT
    
    DMATA1AT<-apply(DMATA1T,2,avecumsum)
    
    DMATAOT<-cbind(DMATA0AT,DMATA1AT)   
    
    
    
    ##Construct restrictions (5) in the main text
    
    visit1<-as.numeric(DATA$T==1)
    visit2<-as.numeric(DATA$T==2)
    visit3<-as.numeric(DATA$T==3)
    visit4<-as.numeric(DATA$T==4)
    visit5<-as.numeric(DATA$T==5)
    visit6<-as.numeric(DATA$T==6)
    visit7<-as.numeric(DATA$T==7)
    visit8<-as.numeric(DATA$T==8)
    visit9<-as.numeric(DATA$T==9)
    visit10<-as.numeric(DATA$T==10)
    
    TEMAT<-cbind(visit1,visit2,visit3,visit4,visit5,visit6,visit7,visit8,visit9,visit10)
    EMAT<-t(TEMAT)
    
    ##Dropout restrictions (15) without intercept for correct covariates
    
    lagfun1<-function(n){
    n1<-n[-1]
    c(n1,0)}
    
    lagfun2<-function(n){
    ave(n,DATADROP$ID,FUN=lagfun1)}
    
    lenID<-ave(DATADROP$ID,DATADROP$ID,FUN=length)
    lenID1<-ave(DATA$ID,DATA$ID,FUN=length)
    maxlen<-max(lenID)
    
    clen<-function(n){
    maxlen:(maxlen-n[1]+1)}
    
    lenvec<-ave(lenID,DATADROP$ID,FUN=clen)
    
    DMATR<-cbind(DATADROP$X1,DATADROP$X2,DATADROP$X3,
    DATADROP$X4,DATADROP$A0P,DATADROP$A1P)
    
    DMATRvec1<-lenvec*DMATR
    
    DMATRvec2<-apply(lenvec*DMATR,2,lagfun2)
    	
    DMATR12<-DMATRvec1-DMATRvec2
    
    DMATRO<-DMATR12[indR0,]              ##DMATRO is a matrix for the coefficients of the weights for the right hand side of the second equation in (15)
    
    indT0<-which(DATADROP$T==1)
    
    Baselinecond<-maxlen*colSums(DMATR[indT0,]) ##Baselinecond is a vector for the left hand side of the second equation in (15) 
    
    
    ##Weight estimation with correct covariates
    
    TDMATAO<-t(DMATAO)
    TDMATRO<-t(DMATRO)
    
    gfunAR<-function(w){            ##Objective function
    lp1<-colSums(TDMATAO*w[1:14])           
    lp2<-colSums(EMAT*w[15:24])
    lp3<-colSums(TDMATRO*w[25:30])
    we<-weightsMLEA*weightsMLER*exp(lp1+lp2+lp3)
    m1<-colMeans(DMATAO*we)                        ##Restrictions (7)
    m2<-(colSums(TEMAT*we)-Npop)/nrow(TEMAT)      ##Restrictions (5)
    m3<-(colSums(DMATRO*we)-Baselinecond)/nrow(DMATRO) ##Restrictions (15)
    c(m1,m2,m3)}
    
    ##Hessian matrix
    
    N1MATAR<-rbind(TDMATAO,EMAT,TDMATRO)
    N2MATAR<-cbind(DMATAO,TEMAT,DMATRO)
    
    DgAR<-function(w){
    lp1<-colSums(TDMATAO*w[1:14])           
    lp2<-colSums(EMAT*w[15:24])
    lp3<-colSums(TDMATRO*w[25:30])
    we<-weightsMLEA*weightsMLER*exp(lp1+lp2+lp3)
    MAT<-N2MATAR*we
    MAT1<-diag(30)
    for (k in 1:30){
    MAT1[,k]<-colMeans(N1MATAR[k,]*MAT)}
    MAT1}
    
    weioptAR<-nleqslv(rep(0,30),gfunAR,DgAR,method="Broyden",
    control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16)))
    
    ##Dropout restrictions (15) without intercept for transformed covariates
    
    DMATRT<-cbind(DATADROP$TX1,DATADROP$TX2,DATADROP$TX3,
    DATADROP$TX4,DATADROP$A0P,DATADROP$A1P)
    
    DMATRvec1T<-lenvec*DMATRT
    
    DMATRvec2T<-apply(lenvec*DMATRT,2,lagfun2)
    	
    DMATR12T<-DMATRvec1T-DMATRvec2T
    
    DMATROT<-DMATR12T[indR0,]  
    
    BaselinecondT<-maxlen*colSums(DMATRT[indT0,])
    
    ##Weight estimation with transformed covariates
    
    TDMATAOT<-t(DMATAOT)
    TDMATROT<-t(DMATROT)
    
    gfunART<-function(w){            ##Objective function
    lp1<-colSums(TDMATAOT*w[1:14])           
    lp2<-colSums(EMAT*w[15:24])
    lp3<-colSums(TDMATROT*w[25:30])
    we<-weightsMLEAT*weightsMLERT*exp(lp1+lp2+lp3)
    m1<-colMeans(DMATAOT*we)                        ##Restrictions (7)
    m2<-(colSums(TEMAT*we)-Npop)/nrow(TEMAT)      ##Restrictions (5)
    m3<-(colSums(DMATROT*we)-Baselinecond)/nrow(DMATROT) ##Restrictions (15)
    c(m1,m2,m3)}
    
    ##Hessian matrix
    
    N1MATART<-rbind(TDMATAOT,EMAT,TDMATROT)
    N2MATART<-cbind(DMATAOT,TEMAT,DMATROT)
    
    DgART<-function(w){
    lp1<-colSums(TDMATAOT*w[1:14])           
    lp2<-colSums(EMAT*w[15:24])
    lp3<-colSums(TDMATROT*w[25:30])
    we<-weightsMLEAT*weightsMLERT*exp(lp1+lp2+lp3)
    MAT<-N2MATART*we
    MAT1<-diag(30)
    for (k in 1:30){
    MAT1[,k]<-colMeans(N1MATART[k,]*MAT)}
    MAT1}
    
    weioptART<-nleqslv(rep(0,30),gfunART,DgART,method="Broyden",
    control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16)))
    
    ##Weight construction
    
    weconsAR<-function(w){
    lp1<-colSums(TDMATAO*w[1:14])           
    lp2<-colSums(EMAT*w[15:24])
    lp3<-colSums(TDMATRO*w[25:30])
    we<-weightsMLEA*weightsMLER*exp(lp1+lp2+lp3)}
    
    weconsART<-function(w){
    lp1<-colSums(TDMATAOT*w[1:14])           
    lp2<-colSums(EMAT*w[15:24])
    lp3<-colSums(TDMATROT*w[25:30])
    we<-weightsMLEAT*weightsMLERT*exp(lp1+lp2+lp3)}
    
    CALW<-weconsAR(weioptAR$x)
    CALWT<-weconsART(weioptART$x)
    
    COEF_CALW<-lm(Y~CA0+CA1,data=DATA,weights=CALW)$coef[2:3]
    COEF_CALWT<-lm(Y~CA0+CA1,data=DATA,weights=CALWT)$coef[2:3]
    
    
    c(cumsum(COEF_MLE),cumsum(COEF_CALW),
    cumsum(COEF_MLET),cumsum(COEF_CALWT))
  }


REmat<-do.call(cbind,RE)

bias<-rowMeans(REmat)-rep(c(10,20),4)
bias

variance<-apply(REmat,1,var)
sqrt(variance)

rmse<-sqrt(bias^2+variance)
rmse

       
