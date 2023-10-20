\name{GenDataSemiContinous}
\alias{GenDataSemiContinous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation of simulated data as explained in the reference manuscript.	      
}
\description{ 
	This function generates data described in the manuscript.
}

\usage{
GenDataSemiContinous(Xreal=FALSE,n=n,p=p,X=NULL,sd=1,impf=20,beta=1,
                               percentOverlap="Full",seed=1)
}

\arguments{
	        \item{Xreal}{ If TRUE, then provide the matrix of covariates X to simulate the semicontinuous response. If FALSE, then X is simulated as well.}
		\item{n}{ Number of individuals}
				\item{p}{ Number of features in X}
		\item{X}{ matrix of the set of features}
			\item{sd}{ Standard deviation of the error of the continuous response}
			\item{impf}{Number of known important features}
			\item{beta}{ Beta regression effect}
				\item{percentOverlap}{Data type. If percentOverlap="Full", then the set of important features is the same between the continous and binary models.  If percentOverlap="Medium", then 50 percent of important features are the same between the two models. If percentOverlap="NoOverlap", then there is no commeon important features between the two models.}
						\item{seed}{ Seed to generate random numbers}
}

\details{
The function will generate data as explained in the manuscript. To see the results, use the ``$" operator. 
} 
\value{
\item{Y}{ A semicontinuous response variable  of dimension \eqn{n}}
\item{X}{ A matrix of features  of dimension \eqn{n\times p}}
\item{Z.cont}{ A binary vector that indicates whether a feature is important in the continous model (vector of dimension \eqn{p})}
\item{Z.bin}{ A binary vector that indicates whether a feature is important in the binary model (vector of dimension \eqn{p})}
}

\references{
Samuel Babatunde, Tolulope Sajobi and Thierry Chekouo (2023), \emph{A Bayesian Variable Selection  for Semicontinuous Response data: Application to cardiovascular disease}, \emph{submitted}.
}
\seealso{
\code{\link{GenDataSemiContinous}}
}

\examples{
library(BVSSemi);
Dat=GenDataSemiContinous(n=500,p=500,sd=1,impf=20,beta=0.3,percentOverlap="Full",seed=1)
str(Dat)
}
