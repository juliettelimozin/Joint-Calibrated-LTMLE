###### R code  for the 2nd simulation study

###### Function to generate data #######

DATA_GEN<-function(ns){  #ns=number of patients
Z1<-rnorm(5*ns,0,1)
Z2<-rnorm(5*ns,0,1)
Z3<-rnorm(5*ns,0,1)
Z4<-rnorm(5*ns,0,1)

A00<-rbinom(ns,1,0.5)# baseline treatment

X1<-rep(0,5*ns)
X2<-rep(0,5*ns)
X3<-rep(0,5*ns)
X4<-rep(0,5*ns)

A0P<-rep(0,5*ns)
CA0P<-rep(0,5*ns)

seq1<-seq(1,5*ns-4,5)
A0P[seq1]<-A00
CA0P[seq1]<-A00

U1<-1-0.3*A00

X1[seq1]<-U1*Z1[seq1]
X2[seq1]<-U1*Z2[seq1]
X3[seq1]<-Z3[seq1]+0.5*A00
X4[seq1]<-Z4[seq1]+0.5*A00

### treatment assignment probability at 1st follow-up
P01<-1/(1+exp(-A00-0.5*X1[seq1]-0.5*X2[seq1]+0.2*X3[seq1]+0.2*X4[seq1]))

A01<-rbinom(ns,1,P01)

seq2<-seq(2,5*ns-3,5)
A0P[seq2]<-A01
CA0P[seq2]<-CA0P[seq1]+A01

U2<-1-0.3*A01

X1[seq2]<-U2*Z1[seq2]
X2[seq2]<-U2*Z2[seq2]
X3[seq2]<-Z3[seq2]+0.5*CA0P[seq2]
X4[seq2]<-Z4[seq2]+0.5*CA0P[seq2]

Y<-rep(0,5*ns)

### mean of the longitudinal outcome at 1st visit given treatment history and covariate history
lp2<-A0P[seq2]+X1[seq2]+X2[seq2]+X3[seq2]+X4[seq2]+X1[seq1]+X2[seq1]+X3[seq1]+X4[seq1]

Y[seq1]<-200+5*lp2+rnorm(ns,0,20)

### treatment assignment probability at 2nd follow-up

P02<-1/(1+exp(-A01-0.5*X1[seq2]-0.5*X2[seq2]+0.2*X3[seq2]+0.2*X4[seq2]))

A02<-rbinom(ns,1,P02)

seq3<-seq(3,5*ns-2,5)
A0P[seq3]<-A02
CA0P[seq3]<-CA0P[seq2]+A02

U3<-1-0.3*A02

X1[seq3]<-U3*Z1[seq3]
X2[seq3]<-U3*Z2[seq3]
X3[seq3]<-Z3[seq3]+0.5*CA0P[seq3]
X4[seq3]<-Z4[seq3]+0.5*CA0P[seq3]

### mean of the longitudinal outcome at 2nd visit given treatment history and covariate history

lp3<-A0P[seq3]+X1[seq3]+X2[seq3]+X3[seq3]+X4[seq3]+X1[seq2]+X2[seq2]+X3[seq2]+X4[seq2]

Y[seq2]<-200+5*lp3+rnorm(ns,0,20)

### treatment assignment probability at 3rd follow-up

P03<-1/(1+exp(-A02-0.5*X1[seq3]-0.5*X2[seq3]+0.2*X3[seq3]+0.2*X4[seq3]))

A03<-rbinom(ns,1,P03)

seq4<-seq(4,5*ns-1,5)
A0P[seq4]<-A03
CA0P[seq4]<-CA0P[seq3]+A03

U4<-1-0.3*A03

X1[seq4]<-U4*Z1[seq4]
X2[seq4]<-U4*Z2[seq4]
X3[seq4]<-Z3[seq4]+0.5*CA0P[seq4]
X4[seq4]<-Z4[seq4]+0.5*CA0P[seq4]

### mean of the longitudinal outcome at 3rd visit given treatment history and covariate history

lp4<-A0P[seq4]+X1[seq4]+X2[seq4]+X3[seq4]+X4[seq4]+X1[seq3]+X2[seq3]+X3[seq3]+X4[seq3]

Y[seq3]<-200+5*lp4+rnorm(ns,0,20)

### treatment assignment probability at 4th follow-up

P04<-1/(1+exp(-A03-0.5*X1[seq4]-0.5*X2[seq4]+0.2*X3[seq4]+0.2*X4[seq4]))

A04<-rbinom(ns,1,P04)

seq5<-seq(5,5*ns,5)
A0P[seq5]<-A04
CA0P[seq5]<-CA0P[seq4]+A04

U5<-1-0.3*A04

X1[seq5]<-U5*Z1[seq5]
X2[seq5]<-U5*Z2[seq5]
X3[seq5]<-Z3[seq5]+0.5*CA0P[seq5]
X4[seq5]<-Z4[seq5]+0.5*CA0P[seq5]

### mean of the longitudinal outcome at 4th visit given treatment history and covariate history

lp5<-A0P[seq5]+X1[seq5]+X2[seq5]+X3[seq5]+X4[seq5]+X1[seq4]+X2[seq4]+X3[seq4]+X4[seq4]

Y[seq4]<-200+5*lp5+rnorm(ns,0,20)

### treatment assignment probability at 5th follow-up

P05<-1/(1+exp(-A04-0.5*X1[seq5]-0.5*X2[seq5]+0.2*X3[seq5]+0.2*X4[seq5]))

A05<-rbinom(ns,1,P05)

Z16<-rnorm(ns,0,1)
Z26<-rnorm(ns,0,1)
Z36<-rnorm(ns,0,1)
Z46<-rnorm(ns,0,1)

U6<-1-0.3*A05

X16<-U6*Z16
X26<-U6*Z26
X36<-Z36+0.5*A05+0.5*CA0P[seq5]
X46<-Z46+0.5*A05+0.5*CA0P[seq5]

### the longitudinal outcome at 5th visit given treatment history and covariate history

Y[seq5]<-200+5*(A05+X16+X26+X36+X46+X1[seq5]+X2[seq5]+X3[seq5]+X4[seq5])+rnorm(ns,0,20)


ID<-rep(1:ns,each=5)

A0<-rep(0,5*ns)
A0[seq1]<-A01
A0[seq2]<-A02
A0[seq3]<-A03
A0[seq4]<-A04
A0[seq5]<-A05

CA0<-ave(A0,ID,FUN=cumsum)+rep(A00,each=5)

### transformed covariates
TX1<-X1^3/9

TX2<-X1*X2

TX3<-log(abs(X3))+4

TX4<-1/(1+exp(X4))

DATA<-data.frame(ID,T=rep(c(1:5),ns),A0,A0P,CA0,X1,X2,X3,X4,TX1,TX2,TX3,TX4,Y)

DATA
}

##########################simulations ########

library(compiler)
enableJIT(3)

library(nleqslv)
library(compiler)
library(CBPS)



simCL<-function(j,Np=2500) #Np=number of patients
{        

enableJIT(3)


DATA<-DATA_GEN(Np)


#########################MLE method for weight estimation#################

###Correct variates

starttime.mle<-Sys.time()
Proj_A0<-glm(A0~A0P,data=DATA,binomial)

Prob_A0<-glm(A0~A0P+X1+X2+X3+X4,data=DATA,binomial)
Prob_A0T<-glm(A0~A0P+TX1+TX2+TX3+TX4,data=DATA,binomial)

weProj_A0<-fitted(Proj_A0)*DATA$A0+(1-DATA$A0)*(1-fitted(Proj_A0))

weProb_A0<-fitted(Prob_A0)*DATA$A0+(1-DATA$A0)*(1-fitted(Prob_A0))
weProb_A0T<-fitted(Prob_A0T)*DATA$A0+(1-DATA$A0)*(1-fitted(Prob_A0T))

weightsMLEN<-weProj_A0/weProb_A0
weightsMLENT<-weProj_A0/weProb_A0T

weightsMLE<-ave(weightsMLEN,DATA$ID,FUN=cumprod)
weightsMLET<-ave(weightsMLENT,DATA$ID,FUN=cumprod)

WeightsMLE<-weightsMLE/mean(weightsMLE)  ### scaling MLE weights
WeightsMLET<-weightsMLET/mean(weightsMLET)

COEF_MLE<-lm(Y~CA0,data=DATA,weights=WeightsMLE)$coef[2]
COEF_MLET<-lm(Y~CA0,data=DATA,weights=WeightsMLET)$coef[2]

endtime.mle<-Sys.time()

###################### calibraion by Type (2) method##############
library(nleqslv)

starttime.MB=Sys.time()

Proj_A0A<-Proj_A0

A0pi3A<-fitted(Proj_A0A)

A0kern3A<-DATA$A0-A0pi3A

A0SQ1T3A<-A0kern3A*DATA$X1
A0SQ2T3A<-A0kern3A*DATA$X2
A0SQ3T3A<-A0kern3A*DATA$X3
A0SQ4T3A<-A0kern3A*DATA$X4
A0SQ5T3A<-A0kern3A*DATA$A0P

A0NSQ1<-ave(A0kern3A,DATA$ID,FUN=cumsum)
A0NSQ2<-ave(A0SQ1T3A,DATA$ID,FUN=cumsum)
A0NSQ3<-ave(A0SQ2T3A,DATA$ID,FUN=cumsum)
A0NSQ4<-ave(A0SQ3T3A,DATA$ID,FUN=cumsum)
A0NSQ5<-ave(A0SQ4T3A,DATA$ID,FUN=cumsum)
A0NSQ6<-ave(A0SQ5T3A,DATA$ID,FUN=cumsum)

A0AMAT3A<-cbind(A0NSQ1,A0NSQ2,A0NSQ3,A0NSQ4,A0NSQ5,A0NSQ6)

ProjEA0A<-fitted(Proj_A0A)*DATA$A0+(1-fitted(Proj_A0A))*(1-DATA$A0)

ProjEA<-ave(ProjEA0A,DATA$ID,FUN=cumprod)

gfunA<-function(w){
Prob<-1/(1+exp(-(2*DATA$A0-1)*(w[1]+w[2]*DATA$A0P+w[3]*DATA$X1+w[4]*DATA$X2+w[5]*DATA$X3+w[6]*DATA$X4)))
bv<-ProjEA/ave(Prob,DATA$ID,FUN=cumprod)
colMeans(A0AMAT3A*bv)}

weioptA<-nleqslv(c(Prob_A0$coef),gfunA,method="Broyden",
control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16)))


##########Transformed covariates

A0SQ1T3AT<-A0kern3A*DATA$TX1
A0SQ2T3AT<-A0kern3A*DATA$TX2
A0SQ3T3AT<-A0kern3A*DATA$TX3
A0SQ4T3AT<-A0kern3A*DATA$TX4

A0AMAT3AT<-cbind(A0NSQ1,ave(A0SQ1T3AT,DATA$ID,FUN=cumsum),
ave(A0SQ2T3AT,DATA$ID,FUN=cumsum),ave(A0SQ3T3AT,DATA$ID,FUN=cumsum),
ave(A0SQ4T3AT,DATA$ID,FUN=cumsum),A0NSQ6)

gfunAT<-function(w){
Prob<-1/(1+exp(-(2*DATA$A0-1)*(w[1]+w[2]*DATA$A0P+w[3]*DATA$TX1+w[4]*DATA$TX2+w[5]*DATA$TX3+w[6]*DATA$TX4)))
bv<-ProjEA/ave(Prob,DATA$ID,FUN=cumprod)
colMeans(A0AMAT3AT*bv)}


weioptAT<-nleqslv(c(Prob_A0T$coef),gfunAT,method="Broyden",
control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16)))


weicalcA<-function(w){
Prob<-1/(1+exp(-(2*DATA$A0-1)*(w[1]+w[2]*DATA$A0P+w[3]*DATA$X1+w[4]*DATA$X2+w[5]*DATA$X3+w[6]*DATA$X4)))
ProjEA/ave(Prob,DATA$ID,FUN=cumprod)}

weicalc1A<-function(w){
Prob<-1/(1+exp(-(2*DATA$A0-1)*(w[1]+w[2]*DATA$A0P+w[3]*DATA$TX1+w[4]*DATA$TX2+w[5]*DATA$TX3+w[6]*DATA$TX4)))
ProjEA/ave(Prob,DATA$ID,FUN=cumprod)}

ASETWN<-weicalcA(weioptA$x)
ASETWTN<-weicalc1A(weioptAT$x)

ASETW<-ASETWN/mean(ASETWN)
ASETWT<-ASETWTN/mean(ASETWTN)

COEF_SETW<-lm(Y~CA0,data=DATA,weights=ASETW)$coef[2]
COEF_SETWT<-lm(Y~CA0,data=DATA,weights=ASETWT)$coef[2]

endtime.MB<-Sys.time()

###################### calibraion by Type (1) method#########################

starttime.cmle<-Sys.time()
visit1<-as.numeric(DATA$T==1)
visit2<-as.numeric(DATA$T==2)
visit3<-as.numeric(DATA$T==3)
visit4<-as.numeric(DATA$T==4)
visit5<-as.numeric(DATA$T==5)

TCCMAT<-cbind(visit1,visit2,visit3,visit4,visit5)
CCMAT<-t(TCCMAT)

Baselinecond<-c(sum(visit1),sum(visit2),sum(visit3),sum(visit4),sum(visit5))

NAMATA<-t(A0AMAT3A)

gfunAD<-function(w){
lp1<-colSums(NAMATA*w[1:6])
lp2<-colSums(CCMAT*w[7:11])
we<-weightsMLE*exp(lp1+lp2)
m1<-colMeans(A0AMAT3A*we)
m2<-(colSums(TCCMAT*we)-Baselinecond)/nrow(TCCMAT)
c(m1,m2)}

weioptAD<-nleqslv(rep(0,11),gfunAD,method="Broyden",
control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16)))

###########################Transformed covariates

NAMATAT<-t(A0AMAT3AT)

gfunADT<-function(w){
lp1<-colSums(NAMATAT*w[1:6])
lp2<-colSums(CCMAT*w[7:11])
we<-weightsMLET*exp(lp1+lp2)
m1<-colMeans(A0AMAT3AT*we)
m2<-(colSums(TCCMAT*we)-Baselinecond)/nrow(TCCMAT)
c(m1,m2)}

weioptADT<-nleqslv(rep(0,11),gfunADT,method="Broyden",
control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16),allowSingular=TRUE))


weconsAD<-function(w){
lp1<-colSums(NAMATA*w[1:6])
lp2<-colSums(CCMAT*w[7:11])
weightsMLE*exp(lp1+lp2)}

weconsADT<-function(w){
lp1<-colSums(NAMATAT*w[1:6])
lp2<-colSums(CCMAT*w[7:11])
weightsMLET*exp(lp1+lp2)}

ADSETW<-weconsAD(weioptAD$x)
ADSETWT<-weconsADT(weioptADT$x)

COEF_ADSETW<-lm(Y~CA0,data=DATA,weights=ADSETW)$coef[2]
COEF_ADSETWT<-lm(Y~CA0,data=DATA,weights=ADSETWT)$coef[2]

endtime.cmle<-Sys.time()

########################CBPS##################

starttime.CBPS<-Sys.time()
library(CBPS)

ind_cbps1<-seq(1,5*Np-4,5)
ind_cbps2<-rep(seq(0,5*Np-5,5),each=2)+1:2
ind_cbps3<-rep(seq(0,5*Np-5,5),each=3)+1:3
ind_cbps4<-rep(seq(0,5*Np-5,5),each=4)+1:4

CBPS_MOD1<-CBMSM(A0~X1+X2+X3+X4,id=DATA$ID[ind_cbps1],time=DATA$T[ind_cbps1],data=DATA[ind_cbps1,],type="MSM",twostep = TRUE,msm.variance = "approx",
time.vary = TRUE)
CBPS_MOD2<-CBMSM(A0~X1+X2+X3+X4,id=DATA$ID[ind_cbps2],time=DATA$T[ind_cbps2],data=DATA[ind_cbps2,],type="MSM",twostep = TRUE,msm.variance = "approx",
time.vary = TRUE)
CBPS_MOD3<-CBMSM(A0~X1+X2+X3+X4,id=DATA$ID[ind_cbps3],time=DATA$T[ind_cbps3],data=DATA[ind_cbps3,],type="MSM",twostep = TRUE,msm.variance = "approx",
time.vary = TRUE)
CBPS_MOD4<-CBMSM(A0~X1+X2+X3+X4,id=DATA$ID[ind_cbps4],time=DATA$T[ind_cbps4],data=DATA[ind_cbps4,],type="MSM",twostep = TRUE,msm.variance = "approx",
time.vary = TRUE)
CBPS_MOD5<-CBMSM(A0~X1+X2+X3+X4,id=DATA$ID,time=DATA$T,data=DATA,type="MSM",twostep = TRUE,msm.variance = "approx",
time.vary = TRUE)

CBPS_weights<-rep(0,Np*5)
CBPS_weights[seq(1,Np*5-4,5)]<-CBPS_MOD1$weights
CBPS_weights[seq(2,Np*5-3,5)]<-CBPS_MOD2$weights
CBPS_weights[seq(3,Np*5-2,5)]<-CBPS_MOD3$weights
CBPS_weights[seq(4,Np*5-1,5)]<-CBPS_MOD4$weights
CBPS_weights[seq(5,Np*5,5)]<-CBPS_MOD5$weights

##Transformed covariates

CBPS_MOD1T<-CBMSM(A0~TX1+TX2+TX3+TX4,id=DATA$ID[ind_cbps1],time=DATA$T[ind_cbps1],data=DATA[ind_cbps1,],type="MSM",twostep=TRUE,
msm.variance = "approx",time.vary = TRUE)
CBPS_MOD2T<-CBMSM(A0~TX1+TX2+TX3+TX4,id=DATA$ID[ind_cbps2],time=DATA$T[ind_cbps2],data=DATA[ind_cbps2,],type="MSM",twostep = TRUE,
msm.variance = "approx",time.vary = TRUE)
CBPS_MOD3T<-CBMSM(A0~TX1+TX2+TX3+TX4,id=DATA$ID[ind_cbps3],time=DATA$T[ind_cbps3],data=DATA[ind_cbps3,],type="MSM",twostep = TRUE,
msm.variance = "approx",time.vary = TRUE)
CBPS_MOD4T<-CBMSM(A0~TX1+TX2+TX3+TX4,id=DATA$ID[ind_cbps4],time=DATA$T[ind_cbps4],data=DATA[ind_cbps4,],type="MSM",twostep = TRUE,
msm.variance = "approx",time.vary = TRUE)
CBPS_MOD5T<-CBMSM(A0~TX1+TX2+TX3+TX4,id=DATA$ID,time=DATA$T,data=DATA,type="MSM",twostep = TRUE,
msm.variance = "approx",time.vary = TRUE)

CBPS_weightsT<-rep(0,Np*5)
CBPS_weightsT[seq(1,Np*5-4,5)]<-CBPS_MOD1T$weights
CBPS_weightsT[seq(2,Np*5-3,5)]<-CBPS_MOD2T$weights
CBPS_weightsT[seq(3,Np*5-2,5)]<-CBPS_MOD3T$weights
CBPS_weightsT[seq(4,Np*5-1,5)]<-CBPS_MOD4T$weights
CBPS_weightsT[seq(5,Np*5,5)]<-CBPS_MOD5T$weights

COEF_CBPS<-lm(Y~CA0,data=DATA,weights=CBPS_weights)$coef[2]
COEF_CBPST<-lm(Y~CA0,data=DATA,weights=CBPS_weightsT)$coef[2]

endtime.CBPS<-Sys.time()

result<-c(COEF_MLE,COEF_SETW,COEF_ADSETW,COEF_CBPS,
COEF_MLET,COEF_SETWT,COEF_ADSETWT,COEF_CBPST,
endtime.mle-starttime.mle,
endtime.MB-starttime.MB,
endtime.cmle-starttime.cmle,
endtime.CBPS-starttime.CBPS)

save(result, file=paste0('/simulation/result_',Np,'_',j, '.Rda'))

}





   
simCLseq<-function(k)
{
    seqin<-(k-1)*40+1:40
    
    for (j in seqin)
    {
         simCL(j, Np=500)
         simCL(j, Np=1000)
         simCL(j, Np=2500)
    }
}


### parallel simulations #######
library(parallel)
no_core<-64
mclapply(1:64, simCLseq, mc.cores=no_core) 

#### summarize simulation results################

CL<-NULL
for(Np in c(500, 1000, 2500))
{
    
    
    for( k in 1:64)
    {
        seqin<-(k-1)*40+1:40
        for (j in seqin)
        {
            load(file=paste0('/simulation/result_',Np,'_',j, '.Rda'))
            CL<-rbind(CL,result)
        }
    }
    
}


corejob=40

###N=500
k=1

index=(k-1)*corejob*64+(1:(64*corejob))
REmat=CL[index[1:2500],]

bias<-colMeans(REmat)-rep(10,8)
bias

mase<-apply(abs(REmat-rep(10,8)),2,median)
mase

variance<-apply(REmat,2,var)
sqrt(variance)

rmse<-sqrt(bias^2+variance)
rmse


result500<-rbind(bias,sqrt(variance),mase,rmse)

###N=1000
k=2

index=(k-1)*corejob*64+(1:(64*corejob))
REmat=CL[index[1:2500],]

bias<-colMeans(REmat)-rep(10,8)
bias

mase<-apply(abs(REmat-rep(10,8)),2,median)
mase

variance<-apply(REmat,2,var)
sqrt(variance)

rmse<-sqrt(bias^2+variance)
rmse

result1000<-rbind(bias,sqrt(variance),mase,rmse)

##2500
k=3

index=(k-1)*corejob*64+(1:(64*corejob))
REmat=CL[index[1:2500],]

bias<-colMeans(REmat)-rep(10,8)
bias

mase<-apply(abs(REmat-rep(10,8)),2,median)
mase

variance<-apply(REmat,2,var)
sqrt(variance)

rmse<-sqrt(bias^2+variance)
rmse

result2500<-rbind(bias,sqrt(variance),mase,rmse)



library(xtable)
xtable(result500[, c(1,3,2,4, 5,7,6,8)])
xtable(result1000[, c(1,3,2,4, 5,7,6,8)])
xtable(result2500[, c(1,3,2,4, 5,7,6,8)])


