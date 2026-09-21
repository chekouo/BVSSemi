test_that("CVPredictBVSSemi requires nmodels", {
  Dat <- GenDataSemiContinous(n = 60, p = 15, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  expect_error(
    CVPredictBVSSemi(Y = Dat$Y, X = Dat$X, K = 3, Method = "BVSSemiMRF",
                      mcmcsample = 100, burnin = 50),
    "nmodels is required"
  )
})

test_that("CVPredictBVSSemi returns one row of metrics per fold", {
  Dat <- GenDataSemiContinous(n = 100, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  K <- 3
  perf <- quietly(CVPredictBVSSemi(Y = Dat$Y, X = Dat$X, K = K, seed = 1, nmodels = 5,
                                    Method = "BVSSemiMRF", mcmcsample = 150, burnin = 75))

  expect_s3_class(perf, "data.frame")
  expect_equal(nrow(perf), K)
  expect_true(all(c("RMSE_combined", "MAE_combined", "AUC_binary") %in% names(perf)))
  expect_true(all(perf$RMSE_combined >= 0, na.rm = TRUE))
})

test_that("ncores > 1 (PSOCK cluster) gives identical results to sequential (regression test)", {
  ## Regression test for a fix where runFold()'s closure carried a stale
  ## reference to the caller's environment (e.g. Dat$Y) through lazy
  ## argument promises, making PSOCK workers fail with
  ## "object 'Dat' not found". Also guards against the parallel path
  ## silently diverging from the sequential one.
  skip_on_cran() ## spawns real worker processes; keep this off CRAN's test machines
  Dat <- GenDataSemiContinous(n = 100, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  K <- 4
  perf_seq <- quietly(CVPredictBVSSemi(Y = Dat$Y, X = Dat$X, K = K, seed = 1, ncores = 1, nmodels = 5,
                                        Method = "BVSSemiMRF", mcmcsample = 150, burnin = 75))
  perf_par <- quietly(CVPredictBVSSemi(Y = Dat$Y, X = Dat$X, K = K, seed = 1, ncores = 2, nmodels = 5,
                                        Method = "BVSSemiMRF", mcmcsample = 150, burnin = 75))

  expect_s3_class(perf_par, "data.frame")
  expect_equal(nrow(perf_par), K)
  expect_equal(perf_seq, perf_par)
})

test_that("ncores > 1 also works with Xcov and extra MainBVSSemi hyperparameters", {
  skip_on_cran()
  n <- 80; p <- 15
  Dat <- GenDataSemiContinous(n = n, p = p, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  Xcov <- matrix(rnorm(n * 2), n, 2)
  perf <- quietly(CVPredictBVSSemi(Y = Dat$Y, X = Dat$X, Xcov = Xcov, K = 3, seed = 5, ncores = 2,
                                    nmodels = 5, Method = "BVSSemiIndep", nu1cont = -4,
                                    mcmcsample = 100, burnin = 50))

  expect_s3_class(perf, "data.frame")
  expect_equal(nrow(perf), 3)
})
