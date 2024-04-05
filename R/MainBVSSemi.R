#methods=BVSSemiMRF, BVSSemiIndep, BVSSemiComb
MainBVSSemi<-function(Method="BVSSemiMRF",Y=Y,X=X,Xcov=NULL,seed=1,atheta=1,btheta=1,tau2cont=1,
                      tau2bin=0.5,nu1cont=-3,nu2bin=-3,varpropTheta=.25,Bigtau2=100,asigma=.1,bsigma=.1,
                      mcmcsample=10000,burnin=5000){
  tau21=tau2cont;tau22=tau2bin
  nu1=nu1cont;nu2=nu2bin
  set.seed(seed)
  n=length(Y)
  N2=which(Y!=0)
  y2=Y[N2]
  if (is.null(Xcov)){
    pc=0
  } else if (is.vector(Xcov)){
    Xcov=matrix(Xcov,n,1)
    pc=ncol(Xcov)
  } else {
    Xcov=as.matrix(Xcov)
    pc=ncol(Xcov)
  }
  X=as.matrix(X)
  p=ncol(X)

  aa=asigma
  ba=bsigma
  #aa=ba=0.001

  #####
  alpha1=atheta;beta1=btheta;
  
  ### Initialization
  Gamma1=rbinom(p,size=1,prob = 1/(1+exp(-nu1))) ### Regression model
  Gamma2=rbinom(p,size=1,prob = 1/(1+exp(-nu2)))# Binary model
  sigma2=.1;
   if (Method=="BVSSemiMRF"){
      theta=.5;
  } else {
  theta=0
     }
  beta=rep(0,p+pc+1)
    U=rep(0,n)
    U[Y!=0]=-abs(rnorm(sum(Y!=0)))
    U[Y==0]=abs(rnorm(sum(Y==0)))

  #### Output values
  NN=mcmcsample
  acceptTheta=0;Gam1Mean=Gam2Mean=rep(0,p);thetaMean=0;
  #betaSample=matrix(0,NN-burnin,p+pc+1);
  thetasample=rep(0,NN-burnin);sigma2Sample=rep(0,NN-burnin);
  
 
  
  NR=length(N2)
  
  for (s in 1:NN){
    if ((Method=="BVSSemiIndep")||(Method=="BVSSemiMRF")){
      ### Sample Gamma1
      Gamma1F=SampleGamma(Gamma1,y2, X[N2,],Xcov[N2,],pc, sigma2,tau2=tau21,Bigtau2,theta,Gamma2,nu1)
      Gamma1=Gamma1F$Gamma
      ### Sample sigma2
      sigma2=Sigma2(NR,Gamma1F$uSu,aa,ba)
        ### Sample Gamma2 binary model
      Gamma2F=SampleGamma(Gamma2,U, X,Xcov,pc,sigma2=1,tau2=tau22,Bigtau2,theta,Gamma1,nu2)
      Gamma2=Gamma2F$Gamma
      beta=Gamma2F$beta
      
    } else if (Method=="BVSSemiComb") {
      GammaF=SampleGammaCombProb(N2,Gamma1,U,Y, X,Xcov,tau2=tau21,Bigtau2,nu1,sigma2)
      beta=GammaF$beta
      Gamma1=GammaF$Gamma
      ### Sample sigma2
      sigma2=Sigma2(NR,GammaF$uSu,aa,ba)
      Gamma2=Gamma1
    }
    if (s>burnin){
      Gam1Mean=Gam1Mean+Gamma1/(NN-burnin)
      Gam2Mean=Gam2Mean+Gamma2/(NN-burnin)
      #betaSample[s-(burnin),]=beta
      sigma2Sample[s-(burnin)]=sigma2
    }
    
      U=YLatent2 (Y, X, Xcov,beta);
    
    if (s%%(NN/5)==1){
      print(paste("Number of mcmc Sample is =",s))
    #  print(beta[1:30])
    }
    if (Method=="BVSSemiMRF"){
      ### Sample theta
      thetaF=SampleTheta(theta, nu1,nu2,Gamma1,Gamma2,alpha1,beta1,varpropo=varpropTheta)
      theta=thetaF$theta
      #thetaMean=thetaMean+theta/NN
      acceptTheta=acceptTheta+thetaF$accept
      if (s>burnin) thetasample[s-(burnin)]= theta;
    } else {
      theta=0
    }
    
  }
  return(list(prob.Z.Cont=Gam1Mean,prob.Z.Bin=Gam2Mean,thetasample=thetasample,AcceptanceRateTheta=acceptTheta/NN,sigma2Sample=sigma2Sample))
}
  
  
  
  
