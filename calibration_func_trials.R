
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
  
  DMATRvec1<-lenvec2*DMATR2
  View(DMATRvec1)
  lagfun1<-function(n)
  {
    n1<-n[-1]
    c(n1,0)
  }
  
  lagfun2<-function(n)
  {
    ave(n,DMATR[,2],FUN=lagfun1)
  }
  
  DMATRvec2<-apply(DMATRvec1,2,lagfun2)
  View(DMATRvec2)
  
  DMATR12<-DMATRvec1-DMATRvec2
  
  DMATRO<-DMATR12[simdatafinal$tall+1<lenID,]         ##DMATRO is a matrix for the coefficients of the weights for the left hand side of the equation in (7)
  
  DMATRO <- DMATRO[,-1] ### NEW 
  View(cbind(DMATRO, simdatafinal[simdatafinal$tall>0, c('ID', 't')]))
  
  indT0<-which(simdatafinal$tall==0)  ### indicator for baseline measurements
  
  Baselinecond<-(maxlen-1)*colSums(DMATR[indT0,-1:-2])##Baselinecond is a vector for the right hand side of the  equation in (7)  
  
  View(Baselinecond)
  
  
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
  
  
  weioptAR<-nleqslv(rep(0,dim(DMATR2)[2]-1),gfunAR,DgAR,method="Broyden",
                    control=list(maxit=10000,ftol=10^(-16),xtol=10^(-16)), jacobian=T)
  
  
  
  
  weconsAR<-function(w){
    lp3<-colSums(TDMATRO*w)
    we<-simdatafinal$weights[simdatafinal$tall>0]*exp(lp3)
  }
  
  CALW<-weconsAR(weioptAR$x)
  
  simdatafinal$Cweights<-simdatafinal$weights
  simdatafinal$Cweights[simdatafinal$tall>0]<-CALW
  
  simdatafinal
}

