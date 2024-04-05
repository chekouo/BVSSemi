
library(truncnorm)

library(gear)

YLatent2 = function(Yobs, X,Xcov,betaR) {
  n=nrow(X)
  ll = rep(0, n)
  ll[Yobs!=0] = -Inf
  uu = rep(0, n)
  uu[Yobs == 0] = Inf
  Xb=cbind(rep(1,n),Xcov,X)%*%betaR
  Ys=sapply(1:n, function(x) rtruncnorm(1,a=ll[x],b=uu[x],mean=Xb[x],sd=1))
  return(Ys)
}

SampleGammaCombProb=function(N2,Gamma,U,Y, X,Xcov,tau2,Bigtau2,nu,sigma2){
  p=length(Gamma)
  if (is.null(Xcov)){
    pc=0
  } else {
    pc=ncol(Xcov)
  }
  beta=rep(0,p+pc+1)
  logEX=nu
  GammaNew=proposalGam(Gamma)
  #  GammaNew=Gamma2
  sumX=sum(log(1+exp(logEX)))
  logprior_old=sum(Gamma*logEX)-sumX
  logprior_new=sum(GammaNew*logEX)-sumX
  loglikOldF=loglik(U, X,Xcov, Gamma, 1,tau2,Bigtau2)
  loglikOld=loglikOldF$logl
  cholMat=loglikOldF$cholMat
  betaMean=loglikOldF$betaMean
  loglikOldFR=loglik(Y[N2], X[N2,],Xcov[N2,], Gamma, sigma2,tau2,Bigtau2)
  loglikOld=loglikOld+loglikOldFR$logl
  uSu=loglikOldFR$uSu;
  loglikNewF=loglik(U, X,Xcov, GammaNew, 1,tau2,Bigtau2)
  #loglikNewF=loglikLogit(Y1, X, Xcov,GammaNew, sigma2r,tau2,Bigtau2,R)
  loglikNew=loglikNewF$logl
  loglikNewFR=loglik(Y[N2], X[N2,], Xcov[N2,], GammaNew, sigma2,tau2,Bigtau2)
  loglikNew=loglikNew+loglikNewFR$logl
  logratio=loglikNew+logprior_new-(loglikOld+logprior_old)
  u1=runif(1,0,1)
  if (log(u1)<logratio){
    Gamma=GammaNew
    uSu=loglikNewFR$uSu
    betaMean=loglikNewF$betaMean
    cholMat=loglikNewF$cholMat
  }
  beta=rep(0,p+pc+1)
  wh=which(Gamma==1)
  pp=sum(Gamma==1)
  UU=rnorm(pp+pc+1)
  Bet=betaMean +backsolve(cholMat,UU)
  beta[1:(1+pc)]=Bet[1:(1+pc)]
  if (pp>=1){
    beta[wh+1+pc]=Bet[(2+pc):(pp+1+pc)]
  }
  return (list(Gamma=Gamma,uSu=uSu,beta=beta))
}

#library(gear)
proposalGam=function(gamma){
  prop=gamma;
  p=length(gamma)
  u=runif(1,0,1)
  id1=which(gamma==1)
  #print(id1)
  L=length(id1)
  if ((u<0.5) || ((L==0) || (L==p))){ ## Add/Delete
    l=sample.int(p,1)
    prop[l]=1-gamma[l]
  } else { ## Swap
    id2=setdiff(1:p,id1)
    l1=sample.int(length(id1),1)
    l2=sample.int(length(id2),1)
    prop[l1]=0; prop[l2]=1;
  }
  return(prop)
}


loglik=function(y, X,Xcov, gamma, sigma2,tau2,Bigtau2){
  n=length(y)
  if (is.null(Xcov)){
    pc=0
  } else if (is.vector(Xcov)){
    pc=1;nn=length(Xcov)
    Xcov=matrix(Xcov,nn,1)
  } else {
    pc=ncol(Xcov)
  }
  XX=cbind(rep(1,n),Xcov,X[,gamma==1])
  p=sum(gamma==1)+pc+1
  u=t(y)%*%XX
  if (p==1){
    XX=matrix(XX,n,1);
    Mat=1/Bigtau2+t(XX)%*%XX
  } else if (p==pc+1){
    Mat=diag(c(rep(1/Bigtau2,1+pc)))+t(XX)%*%XX
  } else {
    Mat=diag(c(rep(1/Bigtau2,1+pc),rep(1/tau2,p-pc-1)))+t(XX)%*%XX
  }
  
  #print(paste("PC=",pc))
  #print(dim(XX))
  #print(p)
  #print(t(XX)%*%XX+diag(c(rep(1/Bigtau2,1+pc),rep(1/tau2,p-pc-1))))
  
  #######
  CholMat=chol(Mat)
  logdet=sum(log(diag(CholMat)^2))
  # betaMean=solve(Mat,t(u))
  betaMean = solve_chol(CholMat, t(u))
  #print(mean(betaMean-betaMean1))
  #betaMean=chol2inv(CholMat)%*%t(u)#solve(Mat,t(u))
  uSu=sum(y^2)-u%*%betaMean
  
  
  #######
  
  #print("BetaReg=")
  #   print(solve(Mat,t(u)))
  loglik1=-(0.5/(sigma2))*(uSu)-0.5*logdet-0.5*n*log(sigma2)-0.5*(pc+1)*log(Bigtau2)-0.5*(p-pc-1)*log(tau2)-0.5*n*log(2*pi)
  return(list(logl=loglik1,uSu=uSu,betaMean=betaMean,cholMat=CholMat))
}

### sample from sigma2

Sigma2=function(n,uSu,aa,ba){
  return(1/rgamma(1,shape=.5*n+aa,rate=0.5*uSu+ba))
}

### Sample Gamma
## prior proba
SampleGamma=function(Gamma,y, X,Xcov,pc, sigma2,tau2,Bigtau2,theta,Gamma2,nu){
  
  if (!is.null(Xcov)){
    Xcov=as.matrix(Xcov,nrow=length(y),ncol=pc)
  }
  
  p=length(Gamma)
  logEX=nu+theta*Gamma2
  GammaNew=proposalGam(Gamma)
  #  GammaNew=Gamma
  logprior_old=sum(Gamma*logEX)#-sum(log(1+exp(logEX)))
  logprior_new=sum(GammaNew*logEX)#-sum(log(1+exp(logEX)))
  loglikOldF=loglik(y, X,Xcov, Gamma, sigma2,tau2,Bigtau2)
  loglikOld=loglikOldF$logl
  #XX=cbind(rep(1,length(y)),Xcov,X[,Gamma==1])
  #loggg=dmvnorm(y,mean=rep(0,length(y)),sigma=sigma2*XX%*%diag(c(rep(Bigtau2,1+pc),rep(tau2,sum(Gamma==1))))%*%t(XX)+sigma2*diag(length(y)),log=T,checkSymmetry=TRUE)
  
  #print(paste("N=",length(y)))
  #print(paste("difflogll=",loggg-loglikOld))
  
  uSu=loglikOldF$uSu;
  betaMean=loglikOldF$betaMean
  cholMat=loglikOldF$cholMat;
  loglikNewF=loglik(y, X, Xcov, GammaNew, sigma2,tau2,Bigtau2)
  loglikNew=loglikNewF$logl
  logratio=loglikNew+logprior_new-(loglikOld+logprior_old)
  u2=runif(1,0,1)
  if (log(u2)<logratio){
    Gamma=GammaNew
    uSu=loglikNewF$uSu
    betaMean=loglikNewF$betaMean
    cholMat=loglikNewF$cholMat
  }
  beta=rep(0,p+pc+1)
  wh=which(Gamma==1)
  pp=sum(Gamma==1)
  UU=rnorm(pp+pc+1)
  Bet=betaMean +sqrt(sigma2)*backsolve(cholMat,UU)
  beta[1:(1+pc)]=Bet[1:(1+pc)]
  if (pp>=1){
    beta[wh+1+pc]=Bet[(2+pc):(pp+1+pc)]
  }
  return (list(Gamma=Gamma,uSu=uSu,beta=beta))
}

## sample theta

SampleTheta=function(theta, nu1,nu2,Gamma1,Gamma2,alpha1,beta1,varpropo=1){
  p=length(Gamma1)
  NormaCost=1+exp(nu1+nu2+theta)+exp(nu1)+exp(nu2)
  logEX=theta*sum(Gamma1*Gamma2)-p*log(NormaCost)
  alphanew=theta^2/varpropo;betanew=theta/varpropo
  thetaProp=rgamma(1,alphanew,rate=betanew)
  NormaCostnew=1+exp(nu1+nu2+thetaProp)+exp(nu1)+exp(nu2)
  logEXnew=thetaProp*sum(Gamma1*Gamma2)-p*log(NormaCostnew)
  logratio=logEXnew+dgamma(thetaProp,shape=alpha1,rate=beta1,log = T)+dgamma(theta,shape=thetaProp^2/varpropo,rate=thetaProp/varpropo,log = T)-logEX-dgamma(theta,shape=alpha1,rate=beta1,log = T)-dgamma(thetaProp,alphanew,rate=betanew,log = T)
  u3=runif(1,0,1)
  acceptTheta=0;
  #print(paste("Alphanew=",alphanew));print(paste("Betanew=",betanew))
  if (log(u3)<logratio){
    theta=thetaProp
    #    print("acceptTheta")
    acceptTheta=1;
  }
  return (list(theta=theta,accept=acceptTheta))
}
















