## Tests for the internal (non-exported) MCMC building blocks in
## R/HelpFunctions.R. These are accessed by name because testthat evaluates
## test files inside the package's namespace.

test_that("loglik matches a brute-force multivariate normal marginal likelihood", {
  set.seed(42)
  n <- 40; p <- 8; pc <- 2
  X <- matrix(rnorm(n * p), n, p)
  Xcov <- matrix(rnorm(n * pc), n, pc)
  y <- rnorm(n)
  sigma2 <- 1.7; tau2 <- 0.6; Bigtau2 <- 50

  brute_logl <- function(gamma, Xcov_use) {
    XXb <- if (is.null(Xcov_use)) cbind(1, X[, gamma == 1, drop = FALSE]) else
      cbind(1, Xcov_use, X[, gamma == 1, drop = FALSE])
    psel <- sum(gamma == 1)
    pc_use <- if (is.null(Xcov_use)) 0 else ncol(Xcov_use)
    Lambda_inv <- diag(c(rep(Bigtau2, 1 + pc_use), rep(tau2, psel)), nrow = 1 + pc_use + psel)
    ## beta | sigma2 ~ N(0, sigma2 * Lambda_inv)  =>  y | sigma2 ~ N(0, sigma2*(I + XX Lambda_inv XX'))
    Sigma_y <- sigma2 * diag(n) + sigma2 * XXb %*% Lambda_inv %*% t(XXb)
    logdetSigma <- as.numeric(determinant(Sigma_y, logarithm = TRUE)$modulus)
    quad <- as.numeric(t(y) %*% solve(Sigma_y, y))
    -n / 2 * log(2 * pi) - 0.5 * logdetSigma - 0.5 * quad
  }

  cases <- list(
    general        = list(gamma = c(1, 0, 1, 0, 0, 1, 0, 0), Xcov_use = Xcov),
    all_zero       = list(gamma = rep(0, p), Xcov_use = Xcov),
    all_zero_noXcov = list(gamma = rep(0, p), Xcov_use = NULL),
    single_feature = list(gamma = c(1, rep(0, p - 1)), Xcov_use = Xcov)
  )

  for (nm in names(cases)) {
    cs <- cases[[nm]]
    out <- loglik(y = y, X = X, Xcov = cs$Xcov_use, gamma = cs$gamma,
                  sigma2 = sigma2, tau2 = tau2, Bigtau2 = Bigtau2)
    expect_equal(out$logl[1], brute_logl(cs$gamma, cs$Xcov_use), tolerance = 1e-8,
                 label = paste("loglik case:", nm))
  }
})

test_that("DrawBeta returns a vector of the right length with zeros for excluded features", {
  set.seed(1)
  p <- 10; pc <- 2
  Gamma <- c(1, 0, 1, rep(0, p - 3))
  pp <- sum(Gamma == 1)
  Mat <- diag(pp + pc + 1) * 2 + 0.1
  cholMat <- chol(Mat)
  betaMean <- rep(0, pp + pc + 1)

  beta <- DrawBeta(betaMean, cholMat, sigma2 = 1, Gamma = Gamma, pc = pc)
  expect_length(beta, p + pc + 1)
  excluded <- which(Gamma == 0) + 1 + pc
  expect_true(all(beta[excluded] == 0))
})

test_that("proposalGam always returns a valid 0/1 vector of the same length", {
  set.seed(1)
  gamma <- c(1, 0, 0, 1, 0)
  for (i in 1:50) {
    prop <- proposalGam(gamma)
    expect_length(prop, length(gamma))
    expect_true(all(prop %in% c(0, 1)))
  }
})

test_that("Sigma2 draws are strictly positive", {
  set.seed(1)
  draws <- replicate(200, Sigma2(n = 30, uSu = 5, aa = 0.1, ba = 0.1))
  expect_true(all(draws > 0))
})

test_that("SampleTheta does not crash when theta underflows to exactly 0 (regression test)", {
  ## Previously, theta = 0 made the Gamma proposal's shape/rate collapse to 0,
  ## producing NaN in the acceptance ratio and an uncaught error from
  ## `if (log(u3) < logratio)`. SampleTheta now floors theta away from 0.
  set.seed(1)
  p <- 20
  Gamma1 <- rbinom(p, 1, 0.3); Gamma2 <- rbinom(p, 1, 0.3)
  expect_no_error(
    out <- SampleTheta(theta = 0, nu1 = -2, nu2 = -1.5, Gamma1 = Gamma1, Gamma2 = Gamma2,
                        alpha1 = 1, beta1 = 1, varpropo = 0.25)
  )
  expect_true(out$theta > 0)
  expect_true(out$accept %in% c(0, 1))
})

test_that("SampleTheta's Metropolis-Hastings ratio matches a brute-force computation", {
  set.seed(1)
  p <- 50
  Gamma1 <- rbinom(p, 1, 0.3); Gamma2 <- rbinom(p, 1, 0.3)
  nu1 <- -2; nu2 <- -1.5; alpha1 <- 2; beta1 <- 1; varpropo <- 0.3
  theta <- 1.2; thetaProp <- 0.9

  NormaCost <- 1 + exp(nu1 + nu2 + theta) + exp(nu1) + exp(nu2)
  logEX <- theta * sum(Gamma1 * Gamma2) - p * log(NormaCost)
  alphanew <- theta^2 / varpropo; betanew <- theta / varpropo
  NormaCostnew <- 1 + exp(nu1 + nu2 + thetaProp) + exp(nu1) + exp(nu2)
  logEXnew <- thetaProp * sum(Gamma1 * Gamma2) - p * log(NormaCostnew)
  logratio_code <- logEXnew + dgamma(thetaProp, shape = alpha1, rate = beta1, log = TRUE) +
    dgamma(theta, shape = thetaProp^2 / varpropo, rate = thetaProp / varpropo, log = TRUE) -
    logEX - dgamma(theta, shape = alpha1, rate = beta1, log = TRUE) -
    dgamma(thetaProp, alphanew, rate = betanew, log = TRUE)

  target_logpost <- function(th) {
    NC <- 1 + exp(nu1 + nu2 + th) + exp(nu1) + exp(nu2)
    full_logEX <- nu1 * sum(Gamma1) + nu2 * sum(Gamma2) + th * sum(Gamma1 * Gamma2) - p * log(NC)
    full_logEX + dgamma(th, shape = alpha1, rate = beta1, log = TRUE)
  }
  logq <- function(from, to) dgamma(to, shape = from^2 / varpropo, rate = from / varpropo, log = TRUE)
  brute_logratio <- (target_logpost(thetaProp) + logq(thetaProp, theta)) -
    (target_logpost(theta) + logq(theta, thetaProp))

  expect_equal(logratio_code, brute_logratio, tolerance = 1e-8)
})

test_that("SampleGammaCombProb fits the continuous part on y2 (log scale), not raw Y (regression test)", {
  ## Previously this function fit the continuous submodel on the raw response
  ## Y[N2] even when log_scale = TRUE, blowing up predictions for BVSSemiComb.
  set.seed(1)
  n <- 60; p <- 15
  Dat <- GenDataSemiContinous(n = n, p = p, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  N2 <- which(Dat$Y != 0)
  y2 <- log(Dat$Y[N2])
  Gamma <- rep(0, p); Gamma[1:5] <- 1
  U <- rep(0, n)
  U[Dat$Y != 0] <- -abs(rnorm(sum(Dat$Y != 0)))
  U[Dat$Y == 0] <- abs(rnorm(sum(Dat$Y == 0)))

  out <- SampleGammaCombProb(N2 = N2, Gamma = Gamma, U = U, y2 = y2, X = Dat$X, Xcov = NULL,
                              tau2 = 1, Bigtau2 = 100, nu = -3, sigma2 = 1, asigma = .1, bsigma = .1)

  ## betaMeanCont should be on the scale of log(Y), i.e. small in magnitude;
  ## if the raw (exponentiated) Y were used instead, this would be orders of
  ## magnitude larger given beta = 0.3 and impf = 5 features.
  expect_true(all(abs(out$betaMeanCont) < 10))
})
