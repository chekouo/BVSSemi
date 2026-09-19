#' @importFrom stats rnorm rbinom rgamma dgamma runif pnorm cor sd quantile
#' @importFrom truncnorm rtruncnorm
#' @importFrom gear solve_chol
NULL

YLatent2 <- function(Yobs, X, Xcov, betaR) {
  n <- nrow(X)
  ll <- rep(0, n)
  ll[Yobs != 0] <- -Inf
  uu <- rep(0, n)
  uu[Yobs == 0] <- Inf
  Xb <- cbind(rep(1, n), Xcov, X) %*% betaR
  Ys <- rtruncnorm(n, a = ll, b = uu, mean = as.vector(Xb), sd = 1)
  return(Ys)
}

SampleGammaCombProb <- function(N2, Gamma, U, y2, X, Xcov, tau2, Bigtau2, nu, sigma2) {
  if (is.null(Xcov)) {
    pc <- 0
  } else {
    pc <- ncol(Xcov)
  }
  logEX <- nu
  GammaNew <- proposalGam(Gamma)
  sumX <- sum(log(1 + exp(logEX)))
  logprior_old <- sum(Gamma * logEX) - sumX
  logprior_new <- sum(GammaNew * logEX) - sumX
  loglikOldF <- loglik(U, X, Xcov, Gamma, 1, tau2, Bigtau2)
  loglikOld <- loglikOldF$logl
  cholMat <- loglikOldF$cholMat

  betaMeanBin <- loglikOldF$betaMean
  loglikOldFR <- loglik(y2, X[N2, , drop = FALSE], Xcov[N2, , drop = FALSE], Gamma, sigma2, tau2, Bigtau2)
  cholMatCont <- loglikOldFR$cholMat
  betaMeanCont <- loglikOldFR$betaMean
  loglikOld <- loglikOld + loglikOldFR$logl
  uSu <- loglikOldFR$uSu
  loglikNewF <- loglik(U, X, Xcov, GammaNew, 1, tau2, Bigtau2)
  loglikNew <- loglikNewF$logl
  loglikNewFR <- loglik(y2, X[N2, , drop = FALSE], Xcov[N2, , drop = FALSE], GammaNew, sigma2, tau2, Bigtau2)
  loglikNew <- loglikNew + loglikNewFR$logl
  logratio <- loglikNew + logprior_new - (loglikOld + logprior_old)
  u1 <- runif(1, 0, 1)
  if (log(u1) < logratio) {
    Gamma <- GammaNew
    uSu <- loglikNewFR$uSu
    betaMeanCont <- loglikNewFR$betaMean
    betaMeanBin <- loglikNewF$betaMean
    cholMat <- loglikNewF$cholMat
    cholMatCont <- loglikNewFR$cholMat
  }
  ## betaBin is unaffected by sigma2 staleness (the binary/probit part is
  ## always evaluated at fixed sigma2 = 1), so it is safe to draw here.
  ## betaCont must NOT be drawn from the sigma2 passed in above (the value
  ## from the previous MCMC iteration): the caller resamples sigma2 from
  ## its fresh full conditional using uSu below, then draws betaCont from
  ## betaMeanCont/cholMatCont using that updated sigma2.
  betaBin <- DrawBeta(betaMeanBin, cholMat, 1, Gamma, pc)
  return(list(Gamma = Gamma, uSu = uSu, betaBin = betaBin,
              betaMeanCont = betaMeanCont, cholMatCont = cholMatCont))
}


proposalGam <- function(gamma) {
  prop <- gamma
  p <- length(gamma)
  u <- runif(1, 0, 1)
  id1 <- which(gamma == 1)
  L <- length(id1)
  if (u < 0.5) { ## Add/Delete
    l <- sample.int(p, 1)
    prop[l] <- 1 - gamma[l]
  } else if (L > 0 && L < p) { ## Swap
    id2 <- setdiff(1:p, id1)
    l1 <- sample(id1, 1)
    l2 <- sample(id2, 1)
    prop[l1] <- 0
    prop[l2] <- 1
  } ## else: Swap infeasible at the boundary; propose no change
  return(prop)
}

### Draw regression coefficients beta | Gamma, sigma2 from the Gaussian
### posterior N(betaMean, sigma2 * Mat^-1), given Mat's Cholesky factor
DrawBeta <- function(betaMean, cholMat, sigma2, Gamma, pc) {
  p <- length(Gamma)
  beta <- rep(0, p + pc + 1)
  wh <- which(Gamma == 1)
  pp <- sum(Gamma == 1)
  UU <- rnorm(pp + pc + 1)
  Bet <- betaMean + sqrt(sigma2) * backsolve(cholMat, UU)
  beta[1:(1 + pc)] <- Bet[1:(1 + pc)]
  if (pp >= 1) {
    beta[wh + 1 + pc] <- Bet[(2 + pc):(pp + 1 + pc)]
  }
  beta
}


loglikLinear <- function(y = y, X = X, Xcov = Xcov, gamma = gamma, asigma, bsigma, 
                tau2 = tau2, Bigtau2 = Bigtau2) {
  n <- length(y)
  Xcov1 <- Xcov
  if (is.null(Xcov)) {
    pc <- 0
  } else if (is.vector(Xcov)) {
    pc <- 1
    nn <- length(Xcov)
    Xcov1 <- matrix(Xcov, nn, 1)
  } else {
    pc <- ncol(Xcov)
  }
  XX <- cbind(rep(1, n), Xcov1, X[, gamma == 1])
  p <- sum(gamma == 1) + pc + 1
  u <- t(y) %*% XX
  if (p == 1) {
    XX <- matrix(XX, n, 1)
    Mat <- 1 / Bigtau2 + t(XX) %*% XX
  } else if (p == pc + 1) {
    Mat <- diag(c(rep(1 / Bigtau2, 1 + pc))) + t(XX) %*% XX
  } else {
    Mat <- diag(c(rep(1 / Bigtau2, 1 + pc), rep(1 / tau2, p - pc - 1))) + t(XX) %*% XX
  }


  #######
  CholMat <- chol(Mat)
  logdet <- sum(log(diag(CholMat)^2))
  betaMean <- solve_chol(CholMat, t(u))
  uSu <- sum(y^2) - u %*% betaMean


  #######
  #   loglik1 <- -(0.5 / (sigma2)) * (uSu) - 0.5 * logdet - 0.5 * n * log(sigma2) - 0.5 * (pc + 1) * log(Bigtau2) - 0.5 * (p - pc - 1) * log(tau2) - 0.5 * n * log(2 * pi)
  loglik1 <- -(n/2+asigma)*log(.5* (uSu)+bsigma) - 0.5 * logdet  - 0.5 * (pc + 1) * log(Bigtau2) - 0.5 * (p - pc - 1) * log(tau2) - 0.5 * n * log(2 * pi)
  return(list(logl = loglik1, uSu = uSu, betaMean = betaMean, cholMat = CholMat))
}


loglik <- function(y = y, X = X, Xcov = Xcov, gamma = gamma, sigma2 = sigma2, tau2 = tau2, Bigtau2 = Bigtau2) {
  n <- length(y)
  Xcov1 <- Xcov
  if (is.null(Xcov)) {
    pc <- 0
  } else if (is.vector(Xcov)) {
    pc <- 1
    nn <- length(Xcov)
    Xcov1 <- matrix(Xcov, nn, 1)
  } else {
    pc <- ncol(Xcov)
  }
  XX <- cbind(rep(1, n), Xcov1, X[, gamma == 1])
  p <- sum(gamma == 1) + pc + 1
  u <- t(y) %*% XX
  if (p == 1) {
    XX <- matrix(XX, n, 1)
    Mat <- 1 / Bigtau2 + t(XX) %*% XX
  } else if (p == pc + 1) {
    Mat <- diag(c(rep(1 / Bigtau2, 1 + pc))) + t(XX) %*% XX
  } else {
    Mat <- diag(c(rep(1 / Bigtau2, 1 + pc), rep(1 / tau2, p - pc - 1))) + t(XX) %*% XX
  }


  #######
  CholMat <- chol(Mat)
  logdet <- sum(log(diag(CholMat)^2))
  betaMean <- solve_chol(CholMat, t(u))
  uSu <- sum(y^2) - u %*% betaMean


  #######
  loglik1 <- -(0.5 / (sigma2)) * (uSu) - 0.5 * logdet - 0.5 * n * log(sigma2) - 0.5 * (pc + 1) * log(Bigtau2) - 0.5 * (p - pc - 1) * log(tau2) - 0.5 * n * log(2 * pi)
  return(list(logl = loglik1, uSu = uSu, betaMean = betaMean, cholMat = CholMat))
}

### sample from sigma2

Sigma2 <- function(n, uSu, aa, ba) {
  return(1 / rgamma(1, shape = .5 * n + aa, rate = 0.5 * uSu + ba))
  # return (1/rltrgamma(n=1, shape=.5*n+aa, rate=0.5*uSu+ba,trunc=2/3))
  # return (1)
}

### Sample Gamma
## prior proba
SampleGamma <- function(GammaM1 = GammaM1, y = y, X = X, Xcov = Xcov, pc = pc, sigma2 = sigma2, tau2 = tau2,
                        Bigtau2 = Bigtau2, theta = theta, GammaM2 = GammaM2, nu = nu) {
  GammaOutput <- GammaM1
  Xcov1 <- Xcov
  if (!is.null(Xcov)) {
    Xcov1 <- as.matrix(Xcov, nrow = length(y), ncol = pc)
  }
  p <- length(GammaM1)
  logEX <- nu + theta * GammaM2
  GammaNew <- proposalGam(GammaM1)
  logprior_old <- sum(GammaM1 * logEX) #-sum(log(1+exp(logEX)))
  logprior_new <- sum(GammaNew * logEX) #-sum(log(1+exp(logEX)))
    loglikOldF <- loglik(y = y, X = X, Xcov = Xcov1, gamma = GammaM1, sigma2 = sigma2, tau2 = tau2, Bigtau2 = Bigtau2)
    loglikOld <- loglikOldF$logl
    uSu <- loglikOldF$uSu
    betaMean <- loglikOldF$betaMean
    cholMat <- loglikOldF$cholMat
    loglikNewF <- loglik(y = y, X = X, Xcov = Xcov1, gamma = GammaNew, sigma2 = sigma2, tau2 = tau2, Bigtau2 = Bigtau2)
    loglikNew <- loglikNewF$logl

  logratio <- loglikNew + logprior_new - (loglikOld + logprior_old)
  u2 <- runif(1, 0, 1)
  if (log(u2) < logratio) {
    GammaOutput <- GammaNew
   
      uSu <- loglikNewF$uSu
      betaMean <- loglikNewF$betaMean
      cholMat <- loglikNewF$cholMat
    
  }
    beta <- rep(0, p + pc + 1)
    wh <- which(GammaOutput == 1)
    pp <- sum(GammaOutput == 1)
    UU <- rnorm(pp + pc + 1)
    Bet <- betaMean + sqrt(sigma2) * backsolve(cholMat, UU)
    beta[1:(1 + pc)] <- Bet[1:(1 + pc)]
    if (pp >= 1) {
      beta[wh + 1 + pc] <- Bet[(2 + pc):(pp + 1 + pc)]
    }
    return(list(GammaM1 = GammaOutput, uSu = uSu, beta = beta))
 
}

## sample theta

SampleTheta <- function(theta = theta, nu1 = nu1, nu2 = nu2, Gamma1 = Gamma1, Gamma2 = Gamma2, alpha1 = alpha1, beta1 = beta1, varpropo = 1) {
  p <- length(Gamma1)
  thetaOutput <- theta
  NormaCost <- 1 + exp(nu1 + nu2 + theta) + exp(nu1) + exp(nu2)
  logEX <- theta * sum(Gamma1 * Gamma2) - p * log(NormaCost)
  alphanew <- theta^2 / varpropo
  betanew <- theta / varpropo
  thetaProp <- rgamma(1, alphanew, rate = betanew)
  NormaCostnew <- 1 + exp(nu1 + nu2 + thetaProp) + exp(nu1) + exp(nu2)
  logEXnew <- thetaProp * sum(Gamma1 * Gamma2) - p * log(NormaCostnew)
  logratio <- logEXnew + dgamma(thetaProp, shape = alpha1, rate = beta1, log = T) + dgamma(theta, shape = thetaProp^2 / varpropo, rate = thetaProp / varpropo, log = T) - logEX - dgamma(theta, shape = alpha1, rate = beta1, log = T) - dgamma(thetaProp, alphanew, rate = betanew, log = T)
  u3 <- runif(1, 0, 1)
  acceptTheta <- 0
  if (log(u3) < logratio) {
    thetaOutput <- thetaProp
    acceptTheta <- 1
  }
  return(list(theta = thetaOutput, accept = acceptTheta))
}


SampleGammaLinear <- function(GammaM1 = GammaM1, y = y, X = X, Xcov = Xcov, pc = pc, asigma, bsigma, tau2 = tau2,
                        Bigtau2 = Bigtau2, theta = theta, GammaM2 = GammaM2, nu = nu) {
  GammaOutput <- GammaM1
  Xcov1 <- Xcov
  if (!is.null(Xcov)) {
    Xcov1 <- as.matrix(Xcov, nrow = length(y), ncol = pc)
  }
  logEX <- nu + theta * GammaM2
  GammaNew <- proposalGam(GammaM1)
  logprior_old <- sum(GammaM1 * logEX) #-sum(log(1+exp(logEX)))
  logprior_new <- sum(GammaNew * logEX) #-sum(log(1+exp(logEX)))
    loglikOldF <- loglikLinear(y = y, X = X, Xcov = Xcov1, gamma = GammaM1, asigma, bsigma, tau2 = tau2, Bigtau2 = Bigtau2)
    loglikOld <- loglikOldF$logl
    uSu <- loglikOldF$uSu
    betaMean <- loglikOldF$betaMean
    cholMat <- loglikOldF$cholMat
    loglikNewF <- loglikLinear(y = y, X = X, Xcov = Xcov1, gamma = GammaNew, asigma, bsigma, tau2 = tau2, Bigtau2 = Bigtau2)
    loglikNew <- loglikNewF$logl

  logratio <- loglikNew + logprior_new - (loglikOld + logprior_old)

  u2 <- runif(1, 0, 1)
  if (log(u2) < logratio) {
    GammaOutput <- GammaNew
      uSu <- loglikNewF$uSu
      betaMean <- loglikNewF$betaMean
      cholMat <- loglikNewF$cholMat
    
  }
  
    return(list(GammaM1 = GammaOutput, uSu = uSu, betaMean = betaMean, cholMat = cholMat))

}
