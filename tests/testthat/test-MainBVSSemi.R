## Small, fast MCMC runs used purely to smoke-test MainBVSSemi's plumbing
## (output shapes, error handling), not to check statistical convergence.

make_small_data <- function(seed = 1, n = 80, p = 20, impf = 5, log_scale = TRUE) {
  GenDataSemiContinous(n = n, p = p, sd = 1, impf = impf, beta = 0.3,
                        percentOverlap = "Full", seed = seed, log_scale = log_scale)
}

test_that("MainBVSSemi runs for all three methods and returns the documented fields", {
  Dat <- make_small_data()
  expected_fields <- c(
    "prob.Z.Cont", "prob.Z.Bin", "beta.Cont.draws", "beta.Bin.draws",
    "sigma2.draws", "theta.draws", "AcceptanceRateTheta", "pc", "p",
    "log_scale", "X.center", "X.scale", "Xcov.center", "Xcov.scale"
  )

  for (meth in c("BVSSemiMRF", "BVSSemiIndep", "BVSSemiComb")) {
    fit <- quietly(MainBVSSemi(Method = meth, Y = Dat$Y, X = Dat$X, seed = 1,
                                mcmcsample = 200, burnin = 100, thin = 5))
    expect_named(fit, expected_fields, ignore.order = TRUE, label = meth)
    expect_length(fit$prob.Z.Cont, ncol(Dat$X))
    expect_length(fit$prob.Z.Bin, ncol(Dat$X))
    expect_true(all(fit$prob.Z.Cont >= 0 & fit$prob.Z.Cont <= 1))
    expect_true(all(fit$prob.Z.Bin >= 0 & fit$prob.Z.Bin <= 1))
    expect_true(all(fit$sigma2.draws > 0))
    expect_equal(ncol(fit$beta.Cont.draws), ncol(Dat$X) + 1)
    expect_equal(ncol(fit$beta.Bin.draws), ncol(Dat$X) + 1)
    expect_equal(nrow(fit$beta.Cont.draws), nrow(fit$beta.Bin.draws))
  }
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
  expect_equal(ncol(fit$beta.Cont.draws), ncol(Dat$X) + 2 + 1)
})
