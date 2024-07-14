# dyn.load("NonLocal.so");
##tau0 = 0.08

LogLikNonLocal <- function(sel = sel, X = X, Xcov=NULL, y=y, alpha = 0.01, lambda = 0.01, tau0 = tau0, tau1 = tau1, r = 1, maxiter = 40) {
    
     if (is.null(Xcov)) {
    K1 <- 0
     } else {
   K1 <- ncol(Xcov)
     }
    k1 <- as.integer(1 + K1 + length(sel))
    N <- nrow(X)
    ng1 <- as.integer(length(sel))
    result <- .C("LogLikNonLocal",
        k1 = k1, K1 = as.integer(K1), ng1 = ng1, N1 = as.integer(N),
        alpha1 = as.double(alpha), psi1 = as.double(lambda), y = as.double(y),
        Ps1 = as.double(as.vector(as.matrix(t(cbind(Xcov, X[, sel]))))),
        betaMax = as.double(as.vector(rep(0, k1))), r1 = as.integer(r),
        h0v = as.double(tau1), hg1 = as.double(tau0),
        maxiter1 = as.integer(maxiter), stop1 = as.double(10^-3),
        Posbeta1 = as.integer(0), LogLik1 = as.double(0)
    )
    ## Posbeta=1 if the regression coefficients are positifs
    ## Posbeta=0 if the regression coefficients are positifs or negatifs
    ## K=number of covariates (prognostic factors)
    ## k=nbr total of selected covariates
    ## ng=nbr of selected genes
    ## Ps matrix of significant covariates
    ## K1 is the number of features that we are not selecting
    ## ng is the number of feature we want to select
    ## k=1+K1+ng
    ## alpha and psi are hyperparameters for sigma2
    ## y is the dependent var
    ## Ps1 is the matrix of selected covariates including the constant
    ## betaMax is the estimated beta
    ## r1 is the r paramater in the pMOM
    ## hv, h1v ann hg1 should be equal, hov should be large (hov), it's the tau for the constant as the constant should be noninformative
    ## stop is the error term to stop with the algorithm and can be set as stop1=10^-3
    ## maxiter1 is the number of iterations to maximize beta
    ##
    # LogLikNonLocal(int *k1, int *K1, int *ng1,int* N1,double *alpha1,double * psi1,double * y, double* Ps1, double * betaMax,int* r1,double* hv,double* h1v,double* h0v,double* hg1,int* maxiter1, double* stop1,_Bool *Posbeta1,double* LogLik1)
    list(LogLik = result$LogLik1, beta = result$betaMax)
}