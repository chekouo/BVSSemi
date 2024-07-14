library(truncnorm)

library(gear)

YLatent2 <- function(Yobs, X, Xcov, betaR) {
  n <- nrow(X)
  ll <- rep(0, n)
  ll[Yobs != 0] <- -Inf
  uu <- rep(0, n)
  uu[Yobs == 0] <- Inf
  Xb <- cbind(rep(1, n), Xcov, X) %*% betaR
  Ys <- sapply(1:n, function(x) rtruncnorm(1, a = ll[x], b = uu[x], mean = Xb[x], sd = 1))
  return(Ys)
}

SampleGammaCombProb <- function(N2, Gamma, U, Y, X, Xcov, tau2, Bigtau2, nu, sigma2) {
  p <- length(Gamma)
  if (is.null(Xcov)) {
    pc <- 0
  } else {
    pc <- ncol(Xcov)
  }
  beta <- rep(0, p + pc + 1)
  logEX <- nu
  GammaNew <- proposalGam(Gamma)
  sumX <- sum(log(1 + exp(logEX)))
  logprior_old <- sum(Gamma * logEX) - sumX
  logprior_new <- sum(GammaNew * logEX) - sumX
  loglikOldF <- loglik(U, X, Xcov, Gamma, 1, tau2, Bigtau2)
  loglikOld <- loglikOldF$logl
  cholMat <- loglikOldF$cholMat
  betaMean <- loglikOldF$betaMean
  loglikOldFR <- loglik(Y[N2], X[N2, ], Xcov[N2, ], Gamma, sigma2, tau2, Bigtau2)
  loglikOld <- loglikOld + loglikOldFR$logl
  uSu <- loglikOldFR$uSu
  loglikNewF <- loglik(U, X, Xcov, GammaNew, 1, tau2, Bigtau2)
  # loglikNewF=loglikLogit(Y1, X, Xcov,GammaNew, sigma2r,tau2,Bigtau2,R)
  loglikNew <- loglikNewF$logl
  loglikNewFR <- loglik(Y[N2], X[N2, ], Xcov[N2, ], GammaNew, sigma2, tau2, Bigtau2)
  loglikNew <- loglikNew + loglikNewFR$logl
  logratio <- loglikNew + logprior_new - (loglikOld + logprior_old)
  u1 <- runif(1, 0, 1)
  if (log(u1) < logratio) {
    Gamma <- GammaNew
    uSu <- loglikNewFR$uSu
    betaMean <- loglikNewF$betaMean
    cholMat <- loglikNewF$cholMat
  }
  beta <- rep(0, p + pc + 1)
  wh <- which(Gamma == 1)
  pp <- sum(Gamma == 1)
  UU <- rnorm(pp + pc + 1)
  Bet <- betaMean + backsolve(cholMat, UU)
  beta[1:(1 + pc)] <- Bet[1:(1 + pc)]
  if (pp >= 1) {
    beta[wh + 1 + pc] <- Bet[(2 + pc):(pp + 1 + pc)]
  }
  return(list(Gamma = Gamma, uSu = uSu, beta = beta))
}

# library(gear)
proposalGam <- function(gamma) {
  prop <- gamma
  p <- length(gamma)
  u <- runif(1, 0, 1)
  id1 <- which(gamma == 1)
  # print(id1)
  L <- length(id1)
  if ((u < 0.5) || ((L == 0) || (L == p))) { ## Add/Delete
    l <- sample.int(p, 1)
    prop[l] <- 1 - gamma[l]
  } else { ## Swap
    id2 <- setdiff(1:p, id1)
    l1 <- sample.int(length(id1), 1)
    l2 <- sample.int(length(id2), 1)
    prop[l1] <- 0
    prop[l2] <- 1
  }
  return(prop)
}


loglikLinear <- function(y = y, X = X, Xcov = Xcov, gamma = gamma, asigma, bsigma, tau2 = tau2, Bigtau2 = Bigtau2) {
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
                        Bigtau2 = Bigtau2, theta = theta, GammaM2 = GammaM2, nu = nu, Local = TRUE) {
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
  if (Local == TRUE) {
    loglikOldF <- loglik(y = y, X = X, Xcov = Xcov1, gamma = GammaM1, sigma2 = sigma2, tau2 = tau2, Bigtau2 = Bigtau2)
    loglikOld <- loglikOldF$logl
    uSu <- loglikOldF$uSu
    betaMean <- loglikOldF$betaMean
    cholMat <- loglikOldF$cholMat
    loglikNewF <- loglik(y = y, X = X, Xcov = Xcov1, gamma = GammaNew, sigma2 = sigma2, tau2 = tau2, Bigtau2 = Bigtau2)
    loglikNew <- loglikNewF$logl
  } else {
    loglikOld <- LogLikNonLocal(sel = which(GammaM1 == 1), X = X, y, alpha = 0.1, lambda = 0.1, tau = 1, tau0 = 10^3, r = 1, maxiter = 40)$LogLik
    loglikNew <- LogLikNonLocal(sel = which(GammaNew == 1), X = X, y, alpha = 0.1, lambda = 0.1, tau = 1, tau0 = 10^3, r = 1, maxiter = 40)$LogLik
    #print(loglikOld)
  }
  logratio <- loglikNew + logprior_new - (loglikOld + logprior_old)
  u2 <- runif(1, 0, 1)
  if (log(u2) < logratio) {
    GammaOutput <- GammaNew
    if (Local == TRUE) {
      uSu <- loglikNewF$uSu
      betaMean <- loglikNewF$betaMean
      cholMat <- loglikNewF$cholMat
    }
  }
  if (Local == TRUE) {
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
  } else {
    return(list(GammaM1 = GammaOutput))
  }
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


SampleGammaLinear <- function(GammaM1 = GammaM1, y = y, X = X, Xcov = Xcov, pc = pc, asigma, bsigma, sigma2=sigma2, tau2 = tau2,
                        Bigtau2 = Bigtau2, theta = theta, GammaM2 = GammaM2, nu = nu, Local = TRUE) {
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
  if (Local == TRUE) {
    loglikOldF <- loglikLinear(y = y, X = X, Xcov = Xcov1, gamma = GammaM1, asigma, bsigma, tau2 = tau2, Bigtau2 = Bigtau2)
    loglikOld <- loglikOldF$logl
    uSu <- loglikOldF$uSu
    betaMean <- loglikOldF$betaMean
    cholMat <- loglikOldF$cholMat
    loglikNewF <- loglikLinear(y = y, X = X, Xcov = Xcov1, gamma = GammaNew, asigma, bsigma, tau2 = tau2, Bigtau2 = Bigtau2)
    loglikNew <- loglikNewF$logl
  } else {
    loglikOld <- LogLikNonLocal(sel = which(GammaM1 == 1), X = X, Xcov=Xcov1, y=y, tau0 = tau2, tau1=Bigtau2)$LogLik
    loglikNew <- LogLikNonLocal(sel = which(GammaNew == 1), X = X, Xcov=Xcov1, y=y, tau0 = tau2, tau1=Bigtau2)$LogLik
    #print(loglikOld)
   # print(loglikNew)
  }
  logratio <- loglikNew + logprior_new - (loglikOld + logprior_old)
  u2 <- runif(1, 0, 1)
  if (log(u2) < logratio) {
    GammaOutput <- GammaNew
    if (Local == TRUE) {
      uSu <- loglikNewF$uSu
      betaMean <- loglikNewF$betaMean
      cholMat <- loglikNewF$cholMat
    }
  }
  if (Local == TRUE) {
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
  } else {
    return(list(GammaM1 = GammaOutput))
  }
}
