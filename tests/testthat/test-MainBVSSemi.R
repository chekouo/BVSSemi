## Small, fast MCMC runs used purely to smoke-test MainBVSSemi's plumbing
## (output shapes, error handling), not to check statistical convergence.

make_small_data <- function(seed = 1, n = 80, p = 20, impf = 5, log_scale = TRUE) {
  GenDataSemiContinous(n = n, p = p, sd = 1, impf = impf, beta = 0.3,
                        percentOverlap = "Full", seed = seed, log_scale = log_scale)
}

test_that("MainBVSSemi runs for all three methods and returns the documented fields", {
  Dat <- make_small_data()
  expected_fields <- c(
    "prob.Z.Cont", "prob.Z.Bin",
    "Gamma.cont.draws", "Gamma.bin.draws", "logPost.draws",
    "Gamma.cont.models", "Gamma.bin.models", "logPost.models",
    "beta.Cont.mode", "beta.Bin.mode", "sigma2.mean", "U.mean",
    "sigma2.draws", "theta.draws", "AcceptanceRateTheta", "pc", "p",
    "log_scale", "X.center", "X.scale", "Xcov.center", "Xcov.scale"
  )
  ## Gam1Mean/Gam2Mean are accumulated by summing Gamma / n_keep over many
  ## iterations, which can land a hair outside [0, 1] due to floating-point
  ## summation error (e.g. 1.0000000000000007); tolerate that, not a defect.
  eps <- 1e-8

  for (meth in c("BVSSemiMRF", "BVSSemiIndep", "BVSSemiComb")) {
    fit <- quietly(MainBVSSemi(Method = meth, Y = Dat$Y, X = Dat$X, seed = 1,
                                mcmcsample = 200, burnin = 100, thin = 5))
    expect_named(fit, expected_fields, ignore.order = TRUE, label = meth)
    expect_length(fit$prob.Z.Cont, ncol(Dat$X))
    expect_length(fit$prob.Z.Bin, ncol(Dat$X))
    expect_true(all(fit$prob.Z.Cont >= -eps & fit$prob.Z.Cont <= 1 + eps))
    expect_true(all(fit$prob.Z.Bin >= -eps & fit$prob.Z.Bin <= 1 + eps))
    expect_true(all(fit$sigma2.draws > 0))
    expect_true(fit$sigma2.mean > 0)
    expect_length(fit$sigma2.mean, 1)

    nModels <- nrow(fit$Gamma.cont.models)
    expect_equal(nrow(fit$Gamma.bin.models), nModels)
    expect_equal(length(fit$logPost.models), nModels)
    expect_equal(ncol(fit$beta.Cont.mode), ncol(Dat$X) + 1)
    expect_equal(ncol(fit$beta.Bin.mode), ncol(Dat$X) + 1)
    expect_equal(nrow(fit$beta.Cont.mode), nModels)
    expect_equal(nrow(fit$beta.Bin.mode), nModels)
    expect_length(fit$U.mean, nrow(Dat$X))
    expect_true(all(is.finite(fit$logPost.models)))
  }
})

test_that("BVSSemiComb shares the same Gamma between its continuous and binary models", {
  Dat <- make_small_data()
  fit <- quietly(MainBVSSemi(Method = "BVSSemiComb", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 200, burnin = 100))
  expect_equal(fit$Gamma.cont.draws, fit$Gamma.bin.draws)
  expect_equal(fit$Gamma.cont.models, fit$Gamma.bin.models)
})

test_that("only BVSSemiMRF estimates theta and an acceptance rate", {
  Dat <- make_small_data()
  fitMRF <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                                 mcmcsample = 200, burnin = 100))
  expect_true(is.numeric(fitMRF$AcceptanceRateTheta))
  expect_true(fitMRF$AcceptanceRateTheta >= 0 && fitMRF$AcceptanceRateTheta <= 1)

  fitIndep <- quietly(MainBVSSemi(Method = "BVSSemiIndep", Y = Dat$Y, X = Dat$X, seed = 1,
                                   mcmcsample = 200, burnin = 100))
  expect_null(fitIndep$AcceptanceRateTheta)
  expect_true(all(fitIndep$theta.draws == 0))
})

test_that("log_scale = TRUE rejects nonpositive nonzero Y values", {
  Dat <- make_small_data(log_scale = FALSE)
  Dat$Y[Dat$Y != 0][1] <- -1  ## force an invalid nonzero value
  expect_error(
    quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, mcmcsample = 50, burnin = 20)),
    "log_scale = TRUE requires"
  )
})

test_that("nu2bin only warns for BVSSemiComb when explicitly supplied (regression test)", {
  Dat <- make_small_data()
  expect_warning(
    quietly(MainBVSSemi(Method = "BVSSemiComb", Y = Dat$Y, X = Dat$X, nu2bin = -2,
                         mcmcsample = 50, burnin = 20)),
    "nu2bin is ignored"
  )
  expect_no_warning(
    quietly(MainBVSSemi(Method = "BVSSemiComb", Y = Dat$Y, X = Dat$X,
                         mcmcsample = 50, burnin = 20))
  )
  expect_no_warning(
    quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, nu2bin = -2,
                         mcmcsample = 50, burnin = 20))
  )
})

test_that("MainBVSSemi accepts forced-in covariates (Xcov)", {
  Dat <- make_small_data()
  Xcov <- matrix(rnorm(nrow(Dat$X) * 2), nrow(Dat$X), 2)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiIndep", Y = Dat$Y, X = Dat$X, Xcov = Xcov,
                              seed = 1, mcmcsample = 100, burnin = 50))
  expect_equal(fit$pc, 2)
  expect_length(fit$Xcov.center, 2)
  expect_length(fit$Xcov.scale, 2)
  expect_equal(ncol(fit$beta.Cont.mode), ncol(Dat$X) + 2 + 1)
  expect_equal(ncol(fit$beta.Bin.mode), ncol(Dat$X) + 2 + 1)
})
