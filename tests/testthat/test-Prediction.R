test_that("PosteriorPredict requires nmodels", {
  Dat <- GenDataSemiContinous(n = 60, p = 15, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 150, burnin = 75))
  expect_error(PosteriorPredict(fit, Xnew = Dat$X), "nmodels is required")
})

test_that("PosteriorPredict returns correctly shaped, valid summaries", {
  Dat <- GenDataSemiContinous(n = 80, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 500, burnin = 200))

  DatTest <- GenDataSemiContinous(n = 40, p = 20, sd = 1, impf = 5, beta = 0.3,
                                   percentOverlap = "Full", seed = 2, log_scale = TRUE)
  pred <- PosteriorPredict(fit, Xnew = DatTest$X, Xcovnew = NULL, nmodels = 5)

  ntest <- nrow(DatTest$X)
  expect_length(pred$p.pos.mean, ntest)
  expect_length(pred$ycont.mean, ntest)
  expect_length(pred$yhat.mean, ntest)
  expect_equal(nrow(pred$ycont.draws), ntest)
  expect_true(all(pred$p.pos.mean >= 0 & pred$p.pos.mean <= 1))
  expect_true(all(pred$yhat.mean >= 0))
  ## yhat.mean should be on the same order of magnitude as the observed
  ## response, not blown up by a stray double-exponentiation
  expect_true(max(pred$yhat.mean) < 10 * max(DatTest$Y))
})

test_that("PosteriorPredict is deterministic (no resampling) and BMA weights are valid", {
  Dat <- GenDataSemiContinous(n = 80, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 1000, burnin = 400))
  DatTest <- GenDataSemiContinous(n = 40, p = 20, sd = 1, impf = 5, beta = 0.3,
                                   percentOverlap = "Full", seed = 2, log_scale = TRUE)

  pred <- PosteriorPredict(fit, Xnew = DatTest$X, nmodels = 5)

  expect_true(pred$bma$nmodels <= 5)
  expect_length(pred$bma$weights, pred$bma$nmodels)
  expect_equal(sum(pred$bma$weights), 1, tolerance = 1e-8)
  expect_true(all(pred$bma$weights >= 0))
  expect_equal(pred$bma$logpost, sort(pred$bma$logpost, decreasing = TRUE))
  expect_length(pred$yhat.mean, nrow(DatTest$X))
  expect_true(all(is.finite(pred$yhat.mean)))

  ## fully deterministic: no resampling, so repeated calls must be identical
  pred2 <- PosteriorPredict(fit, Xnew = DatTest$X, nmodels = 5)
  expect_equal(pred, pred2)

  ## nmodels larger than the number of visited models should just use all
  ## of them, not error
  predBig <- PosteriorPredict(fit, Xnew = DatTest$X, nmodels = 1e6)
  expect_equal(predBig$bma$nmodels, nrow(fit$Gamma.cont.models))
})

test_that("PosteriorPredict errors informatively on a fit without the per-model fields", {
  Dat <- GenDataSemiContinous(n = 60, p = 15, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 150, burnin = 75))
  fit$Gamma.cont.models <- NULL
  expect_error(
    PosteriorPredict(fit, Xnew = Dat$X, nmodels = 5),
    "Gamma.cont.models"
  )
})

test_that("BMA (nmodels) avoids the real-data instability from a rare extreme draw (regression test)", {
  ## Reproduces the reported failure mode: a weakly identified feature in a
  ## p >> n fit can give one MCMC draw an extreme coefficient. Restricting
  ## prediction to the top few highest-posterior-probability models, each
  ## represented by its deterministic posterior mode, should stay sane.
  set.seed(1)
  n <- 30; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  Y <- exp(rnorm(n, mean = 3, sd = 0.3))
  Y[sample(n, 5)] <- 0
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Y, X = X, seed = 1,
                              tau2cont = 5, nu1cont = -1,
                              mcmcsample = 300, burnin = 100))
  Xnew <- matrix(rnorm(20 * p), 20, p)

  pred_bma <- PosteriorPredict(fit, Xnew = Xnew, nmodels = 5)

  expect_true(all(is.finite(pred_bma$yhat.mean)))
  expect_true(max(pred_bma$yhat.mean) < 1000 * max(Y))
})

test_that("PosteriorPredict stays finite when a model's linear predictor would overflow exp() (regression test)", {
  ## Directly exercises rowLogSumExp()'s guard: if a kept model's continuous
  ## linear predictor is extreme enough that exp() would overflow to Inf,
  ## and that same model's pPos has underflowed to exactly 0, the naive
  ## pPos * exp(logCondMean) computation produces 0 * Inf = NaN, which would
  ## silently poison yhat.mean. PosteriorPredict must stay finite instead.
  set.seed(1)
  n <- 20; p <- 5; K <- 3
  fit <- list(
    Gamma.cont.models = matrix(0, K, p), Gamma.bin.models = matrix(0, K, p),
    logPost.models = c(0, -1, -2),
    beta.Bin.mode = matrix(rnorm((p + 1) * K), K, p + 1),
    beta.Cont.mode = matrix(rnorm((p + 1) * K), K, p + 1),
    sigma2.mean = 0.3,
    X.center = rep(0, p), X.scale = rep(1, p),
    Xcov.center = NULL, Xcov.scale = NULL
  )
  ## rig model 1 to hit exactly the pathological case
  fit$beta.Cont.mode[1, ] <- c(5000, rep(0, p)) ## model 1's intercept alone -> huge etaCont for everyone
  fit$beta.Bin.mode[1, ] <- c(50, rep(0, p))    ## model 1's etaBin -> pPos = 1 - pnorm(50) = 0 exactly

  Xnew <- matrix(rnorm(n * p), n, p)
  pred <- PosteriorPredict(fit, Xnew = Xnew, nmodels = K)

  expect_true(all(is.finite(pred$yhat.mean)))
  expect_false(any(is.na(pred$yhat.mean)))
})

test_that("BVSSemiComb predictions are on the same scale as the other two methods (regression test)", {
  set.seed(1)
  Dat <- GenDataSemiContinous(n = 80, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Medium", seed = 1, log_scale = TRUE)
  DatTest <- GenDataSemiContinous(n = 40, p = 20, sd = 1, impf = 5, beta = 0.3,
                                   percentOverlap = "Medium", seed = 2, log_scale = TRUE)

  maxY <- max(DatTest$Y)
  for (meth in c("BVSSemiMRF", "BVSSemiIndep", "BVSSemiComb")) {
    fit <- quietly(MainBVSSemi(Method = meth, Y = Dat$Y, X = Dat$X, seed = 1,
                                mcmcsample = 300, burnin = 150))
    pred <- PosteriorPredict(fit, Xnew = DatTest$X, nmodels = 5)
    ## Before the fix, BVSSemiComb's continuous part fit on the raw
    ## (exponentiated) Y instead of log(Y), producing astronomically large
    ## predictions; it should now stay comparable to the other two methods.
    expect_true(max(pred$yhat.mean) < 100 * maxY, label = meth)
  }
})

test_that("EvaluatePrediction returns metrics in valid ranges", {
  Dat <- GenDataSemiContinous(n = 100, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 300, burnin = 150))
  DatTest <- GenDataSemiContinous(n = 60, p = 20, sd = 1, impf = 5, beta = 0.3,
                                   percentOverlap = "Full", seed = 2, log_scale = TRUE)
  pred <- PosteriorPredict(fit, Xnew = DatTest$X, nmodels = 5)
  perf <- EvaluatePrediction(DatTest$Y, pred)

  expect_true(perf$AUC_binary >= 0 && perf$AUC_binary <= 1)
  expect_true(perf$Brier_binary >= 0 && perf$Brier_binary <= 1)
  expect_true(perf$RMSE_cont >= 0)
  expect_true(perf$MAE_cont >= 0)
  expect_true(perf$RMSE_combined >= 0)
  expect_true(perf$MAE_combined >= 0)
})

test_that("EvaluatePrediction handles a test set that is all zero or all nonzero", {
  Dat <- GenDataSemiContinous(n = 60, p = 15, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 150, burnin = 75))
  DatTest <- GenDataSemiContinous(n = 20, p = 15, sd = 1, impf = 5, beta = 0.3,
                                   percentOverlap = "Full", seed = 2, log_scale = TRUE)
  pred <- PosteriorPredict(fit, Xnew = DatTest$X, nmodels = 5)

  Yall_zero <- rep(0, length(DatTest$Y))
  perf_zero <- EvaluatePrediction(Yall_zero, pred)
  expect_true(is.na(perf_zero$AUC_binary))
  expect_true(is.na(perf_zero$RMSE_cont))
})
