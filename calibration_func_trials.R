
##### expit function
expit<-function(x){exp(x)/(1+exp(x))}


##Remove observations after dropout
indfun<-function(n)
{
  if (sum(n)==0) {rep(0,9)}
  else{k<-min(which(n==1))
  c(rep(0,k),rep(1,9-k))}
}


#### recalculate weights, with correct covariates, selection function

newweight<-function(alpha=c(2, -1, 0.5, -0.25, -0.1,0.2), gamma=1,x=dropoutdata[,c('X1', 'X2','X3', 'X4', 'y_c', 'ycomp_c')])
{
  x<-cbind(1,x)
  w=expit(as.matrix(x)%*%as.matrix(c(alpha,gamma)))
  1/w
} 


#### estimating equations in (4), using correct covariates and previous outcome

score<-function(alpha=c(2, -1, 0.5, -0.25, -0.1, 0.2), gamma=1, x=dropoutdata[,c('X1', 'X2','X3', 'X4', 'y_c', 'ycomp_c')], ind=dropoutdata$drop)
{
  ### assume selection function known
  
  x<-cbind(1,x)
  expx<-exp(-as.matrix(x)%*%as.matrix(c(alpha,gamma)))
  
  obsx<-x[ind==0,-(length(alpha)+1)]
  expxx<-(expx%*%t(as.matrix(rep(1, dim(x)[2]-1))))
  obsxx<-obsx*expxx[ind==0,]
  
  missx<--x[ind==1,-(length(alpha)+1)]
  
  allscore<-rbind(obsxx,missx)
  
  finalscore<-colSums(allscore)
  as.matrix(finalscore)
}



####first derivative of  the estimating equations in (4), using correct covariates and previous outcome



scoreder<-function(alpha=c(2, -1, 0.5, -0.25, -0.1, 0.2), gamma=1, x=dropoutdata[,c('X1', 'X2','X3', 'X4', 'y_c', 'ycomp_c')], ind=dropoutdata$drop)
{
  
  
  x<-cbind(1,x)
  expx<-exp(-as.matrix(x)%*%as.matrix(c(alpha,gamma)))
  obsx<-x[ind==0,-c((length(alpha)+1))]
  expxx<-expx[ind==0]
  
  mat=0
  for (i in 1:length(expxx))
  {
    xvec<-obsx[i,]
    
    mat=mat-(t(as.matrix(xvec))%*%as.matrix(xvec))*expxx[i]
  }
  mat
}



# Function to calculate the difference-quotient approx gradient
# (matrix) of an arbitrary input (vector) function infcn
# Now recoded to use central differences !

Gradmat<-function(parvec, infcn, eps = 1e-06)
{
  
  dd = length(parvec)
  aa = length(infcn(parvec))
  epsmat = (diag(dd) * eps)/2
  gmat = array(0, dim = c(aa, dd))
  for(i in 1:dd)
    gmat[, i] = (infcn(parvec + epsmat[, i]) -
                   infcn(parvec - epsmat[, i]))/eps
  if(aa > 1) gmat else c(gmat)
}


###### Newton-Raphson 
NRroot<-function(inipar, infcn, nmax = 25, stoptol = 1e-05,
                 eps = 1e-06, gradfunc = NULL)
{
  if(is.null(gradfunc))
    gradfunc = function(x) Gradmat(x, infcn, eps)
  ctr = 0
  newpar = inipar
  oldpar = inipar - 1
  while(ctr < nmax & sqrt(sum((newpar - oldpar)^2)) > stoptol) {
    oldpar = newpar
    newpar = oldpar - solve(gradfunc(oldpar), infcn(oldpar))
    ctr = ctr + 1
  }
  result<-list(nstep = ctr, initial = inipar, final = newpar,
               funcval = infcn(newpar))
  newpar
}



calibration<-function(simdatafinal, var=c('tall', 'X1', 'X2', 'X3', 'X4'))  
{
  
  ##restrictions (7) #########
  lenID<-ave(simdatafinal$sub,simdatafinal$sub,FUN=length)
  maxlen<-max(lenID)
  
  clen<-function(n){
    maxlen:(maxlen-n[1]+1)}
  
  lenvec<-ave(lenID,simdatafinal$sub,FUN=clen)
  lenvec2<-lenvec-1   ### remove the record corresponding to the last obs 
  
  DMATR<-cbind(1, simdatafinal$sub,  simdatafinal[, var] )
  DMATR2<-DMATR[,-2] ### remove the sub index
  
  DMATRvec1<-lenvec2*DMATR2 #(T - t+1)X_t-1 
  
  lagfun1<-function(n)
  {
    n1<-n[-1]
    c(n1,0)
  }
  
  lagfun2<-function(n)
  {
    ave(n,DMATR[,2],FUN=lagfun1)
  }
  
  DMATRvec2<-apply(DMATRvec1,2,lagfun2) #
  
  DMATR12<-simdatafinal$RA*(DMATRvec1-DMATRvec2)
  
  DMATRO<-DMATR12[simdatafinal$tall+1<lenID,]  ##DMATRO is a matrix for the coefficients of the weights for the left hand side of the equation in (7)
  DMATRO <- DMATRO[,-1] ### NEW 
  indT0<-which(simdatafinal$tall==0)  ### indicator for baseline measurements
  
  Baselinecond<-(maxlen-1)*colSums(DMATR[indT0,-1:-2])## NEW ##Baselinecond is a vector for the right hand side of the  equation in (7)  
  
  
  ##Weight estimation with correct covariates##########
  
  TDMATRO<-t(DMATRO)
  
  ##Objective function
  gfunAR<-function(w)
  {            
    lp3<-colSums(TDMATRO*w)
    we<-simdatafinal$weights[simdatafinal$tall>0]*exp(lp3)
    m3<-(colSums(DMATRO*we)-Baselinecond)/nrow(DMATRO) ##Restrictions (6)
    c(m3)
  }
  
  
  ##Hessian matrix
  
  N1MATAR<-TDMATRO
  N2MATAR<-DMATRO
  
  DgAR<-function(w){
    
    lp3<-colSums(TDMATRO*w)
    we<-simdatafinal$weights[simdatafinal$tall>0]*exp(lp3)
    MAT<-N2MATAR*we
    MAT1<-diag(length(w))
    for (k in 1:length(w)){
      MAT1[,k]<-colMeans(N1MATAR[k,]*MAT)}
    MAT1
  }
  
  
  weioptAR<-nleqslv::nleqslv(rep(0,dim(DMATR2)[2]-1),gfunAR,DgAR,method="Broyden",
                             control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16)), jacobian=T)
  
  
  
  
  weconsAR<-function(w){
    lp3<-colSums(TDMATRO*w)
    we<-simdatafinal$weights[simdatafinal$tall>0]*exp(lp3)
    return(we)
  }
  
  CALW<-weconsAR(weioptAR$x)
  
  simdatafinal$Cweights<-simdatafinal$weights
  simdatafinal$Cweights[simdatafinal$tall>0]<-CALW
  
  list(data = simdatafinal, 
       objective.IPW =  gfunAR(rep(0,dim(DMATR2)[2]-1)), 
       objective.Cali = weioptAR$fvec)
}


calibration_by_time<-function(simdatafinal, var=c('A1', 'A1X1')){
  T <- max(simdatafinal$tall)
  data1 <- simdatafinal[simdatafinal$tall == 1,]
  data1 <- data1$RA*data1[, var]
  Tdata1<- t(data1)
  data0 <- simdatafinal[simdatafinal$tall == 0,]
  data0 <- data0$RA*simdatafinal[simdatafinal$tall == 1, var]
  RHS <- colSums(data0)
  print(RHS)
  restrictions_weight1 <- function(w){
    ##Objective function
    lp3<-colSums(Tdata1*w)
    we<-simdatafinal$weights[simdatafinal$tall == 1]*exp(lp3)
    m3<-(colSums(data1*we) - RHS)/nrow(data1) ##Restrictions (6)
    c(m3)
  }
  
  ##Hessian matrix
  
  N1MATAR<-Tdata1
  N2MATAR<-data1
  
  DgAR1<-function(w){
    
    lp3<-colSums(Tdata1*w)
    we<-simdatafinal$weights[simdatafinal$tall == 1]*exp(lp3)
    MAT<-N2MATAR*we
    MAT1<-diag(length(w))
    for (k in 1:length(w)){
      MAT1[,k]<-colMeans(N1MATAR[k,]*MAT)}
    MAT1
  }
  
  wei1optAR<-nleqslv::nleqslv(rep(0,dim(data1)[2]),restrictions_weight1,DgAR1,method="Broyden",
                              control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16)), jacobian=T)
  
  weconsAR<-function(w){
    lp3<-colSums(Tdata1*w)
    we<-simdatafinal$weights[simdatafinal$tall==1]*exp(lp3)
    return(we)
  }
  
  CALW1<-weconsAR(wei1optAR$x)
  
  simdatafinal$Cweights<-simdatafinal$weights
  simdatafinal$Cweights[simdatafinal$tall == 1]<-CALW1
  
  
  for (k in 2:T){
    data1 <- simdatafinal[simdatafinal$tall == k,]
    data1 <- data1$RA*data1[, var]
    Tdata1<- t(data1)
    
    data0 <- simdatafinal[simdatafinal$tall == k-1,]
    data0 <- data0$RA*data0$Cweights*simdatafinal[simdatafinal$tall == k, var]
    RHS <- colSums(data0)
    print(RHS)
    restrictions_weightk <- function(w){
      ##Objective function
      lp3<-colSums(Tdata1*w)
      we<-simdatafinal$weights[simdatafinal$tall == k]*exp(lp3)
      m3<-(colSums(data1*we) - RHS)/nrow(data1) ##Restrictions (6)
      c(m3)
    }
    
    ##Hessian matrix
    
    N1MATAR<-Tdata1
    N2MATAR<-data1
    
    DgAR1<-function(w){
      
      lp3<-colSums(Tdata1*w)
      we<-simdatafinal$weights[simdatafinal$tall == k]*exp(lp3)
      MAT<-N2MATAR*we
      MAT1<-diag(length(w))
      for (j in 1:length(w)){
        MAT1[,j]<-colMeans(N1MATAR[j,]*MAT)}
      MAT1
    }
    
    weikoptAR<-nleqslv::nleqslv(rep(0,dim(data1)[2]),restrictions_weightk,DgAR1,method="Broyden",
                                control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16)), jacobian=T)
    
    weconsAR<-function(w){
      lp3<-colSums(Tdata1*w)
      we<-simdatafinal$weights[simdatafinal$tall==k]*exp(lp3)
      return(we)
    }
    
    CALW1<-weconsAR(weikoptAR$x)
    
    simdatafinal$Cweights[simdatafinal$tall == k]<-CALW1
    
    
  }
  list(data = simdatafinal, 
       objective.IPW =  restrictions_weight1(rep(0,dim(data1)[2])), 
       objective.Cali = wei1optAR$fvec)
}
