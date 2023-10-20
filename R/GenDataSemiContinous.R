GenDataSemiContinous<-function(Xreal=FALSE,n=n,p=p,X=NULL,sd=1,impf=20,beta=1,
                               percentOverlap="Full",seed=1){
  set.seed(seed)
  # impf is the number of important covariates. if impf=10 then the 10 first features are important
  # percentOverlap is the percentage of overlap important features between model 1 and model 2: "Full", "Medium" and "NoOverlap"
  Gamma1=rep(0,p) ### continuous model
  Gamma2=rep(0,p)
  Gamma2[1:impf]=1 ## binary model
  if (Xreal==T){ if (is.null(X)){
    stop("Provide the matrice of the set of features")
  } } else {
    X=matrix(rnorm(n*p),n,p)
  }
  X=scale(X)
    Ystar=X[,1:impf]%*%c(rep(beta,impf))+rnorm(n)
    Y=1-(Ystar>0)
  N1=which(Y==0)
  N2=setdiff(1:n,N1)
  error=sd*rnorm(length(N2));
  if (percentOverlap=="Full"){
    Y[N2]=X[N2,1:impf]%*%c(rep(beta,impf))+error
    Gamma1[1:impf]=1
  }
  if (percentOverlap=="Medium"){
    Y[N2]=X[N2,1:(impf/2)]%*%c(rep(beta,impf/2))+X[N2,(1+impf):(impf+impf/2)]%*%c(rep(beta,impf/2))+error
    Gamma1[1:(impf/2)]=1;Gamma1[(1+impf):(impf+impf/2)]=1
  }
  if (percentOverlap=="NoOverlap"){
    Y[N2]=X[N2,(1+impf):(2*impf)]%*%c(rep(beta,impf))+error
    Gamma1[(1+impf):(2*impf)]=1
  }
  return(list(Y=Y,X=X,Z.cont=Gamma1,Z.bin=Gamma2))
}













  
