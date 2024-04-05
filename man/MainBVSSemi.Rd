\name{MainBVSSemi}
\alias{MainBVSSemi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	An MCMC algorithm to perform a Bayesian variable selection for a semicontinuous response via a two-part model: continous model (i.e. linear model with continuous response) and binary model (Probit model with binary response). 
}
\description{
The algorithm is implemented on the three methods: i) BVSSemiMRF: this method encourages the common selection of important features between the two models (continuous and binary model).  ii) BVSSemiComb: this method selects the same set of features associated for both the continuous and binary models such that each selected feature is related to our semicontinuous response (both zeros and positive values) and iii) BVSSemiIndep: it  assumes that the set of selected features for the log-linear model is not necessarily the same with the logistic model, and the two models are fitted independently.  The algorithm computes mainly the marginal posterior probabilities of inclusion of each feature  for each model.  
}

\usage{
MainBVSSemi(Method="BVSSemiMRF",Y=Y,X=X,Xcov=NULL,seed=1,atheta=1,btheta=1,tau2cont=1,
                      tau2bin=0.5,nu1cont=-3,nu2bin=-3,varpropTheta=.25,Bigtau2=100,asigma=.1,bsigma=.1,
                      mcmcsample=10000,burnin=5000)
}

\arguments{
	\item{Method}{ It's one of the three methods: "BVSSemiMRF","BVSSemiComb" or "BVSSemiIndep". The default value is "BVSSemiMRF".
	}
		\item{Y}{ A semicontinuous response.   
	}
	\item{X}{ A matrix of features  of dimension \eqn{n\times p}}
	\item{Xcov}{ A set of covariates that are not subject of variable selection (e.g. clinical/demographic variables such as sex, age, race, etc..)}
	\item{seed}{ Set a seed number to generate distributions in the MCMC algorithm.}
	\item{atheta}{ It's the shape (hyper)parameter of the gamma prior distribution of theta from the BVSSemiMRF method. The parameter theta measures the borrows strength between the two model and encourages the common selection between the two models.}
		\item{btheta}{ It's the scale (hyper)parameter of the gamma prior distribution of theta from the BVSSemiMRF method.}
		\item{tau2cont}{ It's the variance (hyper)parameter of the normal prior distribution of regression effects on the continuous model}
		\item{tau2bin}{ It's the variance (hyper)parameter of the normal prior distribution of regression effects on the binary model}
		\item{nu1cont}{ log-odds of prior prob. of feature inclusion in the continuous model}
			\item{nu2bin}{ log-odds of prior prob. of feature inclusion in the binary model}
		\item{varpropTheta}{ Variance of the proposal distribution (in the Metro. Hasting step) of theta in our MCMC algorithm. It should be chosen to have an acceptance rate of theta around 20-60 percent. The default value is .25.}
			\item{Bigtau2}{ It's the variance (hyper)parameter of the normal prior distribution of regression effects of features that are not subject of variable selection (e.g., intercept and other clinical covariates)}
\item{asigma}{ Shape parameter for the inverse gamma prior distribution of sigma2}
\item{bsigma}{ Scale parameter for the inverse gamma prior distribution of sigma2}
\item{mcmcsample}{ Total number of MCMC draws. It must be larger than burnin.}
\item{burnin}{ 	Number of draws to discard for burn-in.} 

}
\details{
	The function will return several R objects, which can be assigned to a variable. To see the results, use the ``$" operator. 
}

%  return(list(Z.Cont=Gam1Mean,Z.Bin=Gam2Mean,thetasample=thetasample,AcceptanceRateTheta=acceptTheta/NN,sigma2Sample=sigma2Sample,betaSampleLogit=betaSample))

\value{
\item{prob.Z.Cont}{ Marginal posterior probabilities of each feauture in the continuous model.}
\item{prob.Z.Bin}{ Marginal posterior probabilities of each feauture in the binary model.}
\item{thetasample}{ Posterior samples of theta from the BVSSemiMRF method. }
\item{AcceptanceRateTheta}{ MCMC Acceptance rate of theta from the BVSSemiMRF.}
\item{sigma2Sample}{ Posterior samples of the variance sigma2 from the continous model. }

}

\references{
Samuel Babatunde, Tolulope Sajobi and Thierry Chekouo (2023), \emph{A Bayesian Variable Selection  for Semicontinuous Response data: Application to cardiovascular disease}, \emph{submitted}.
}

\seealso{
\code{\link{GenDataSemiContinous}}
}

\examples{

library(BVSSemi);
Dat=GenDataSemiContinous(n=500,p=200,sd=1,impf=20,beta=0.3,percentOverlap="Full",seed=1)
str(Dat)
result<-MainBVSSemi(Method="BVSSemiMRF",Y=Dat$Y,X=Dat$X,Xcov=NULL,seed=1,atheta=1,btheta=1,tau2cont=1,tau2bin=1,nu1cont=-4,nu2bin=-4,varpropTheta=.25,Bigtau2=100,mcmcsample=50000,burnin=10000);
str(result)
library(pROC)
AUC1=as.numeric(auc(roc(as.factor(Dat$Z.cont),result$prob.Z.Cont)))
AUC2=as.numeric(auc(roc(as.factor(Dat$Z.bin),result$prob.Z.Bin)))
print(AUC1);print(AUC2)

}

