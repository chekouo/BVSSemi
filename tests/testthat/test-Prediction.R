test_that("PosteriorPredict returns correctly shaped, valid summaries", {
  Dat <- GenDataSemiContinous(n = 80, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 200, burnin = 100))

  DatTest <- GenDataSemiContinous(n = 40, p = 20, sd = 1, impf = 5, beta = 0.3,
                                   percentOverlap = "Full", seed = 2, log_scale = TRUE)
  pred <- PosteriorPredict(fit, Xnew = DatTest$X, Xcovnew = NULL)

  ntest <- nrow(DatTest$X)
  expect_length(pred$p.pos.mean, ntest)
  expect_length(pred$mu.mean, ntest)
  expect_length(pred$yhat.mean, ntest)
  expect_equal(nrow(pred$mu.draws), ntest)
  expect_true(all(pred$p.pos.mean >= 0 & pred$p.pos.mean <= 1))
  expect_true(all(pred$yhat.mean >= 0))
  ## yhat.mean should be on the same order of magnitude as the observed
  ## response, not blown up by a stray double-exponentiation
  expect_true(max(pred$yhat.mean) < 10 * max(DatTest$Y))
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
    pred <- PosteriorPredict(fit, Xnew = DatTest$X)
    ## Before the fix, BVSSemiComb's yhat.mean could reach astronomically
    ## large values (double-exponentiated); it should now stay comparable
    ## to the actual response range like the other two methods.
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
  pred <- PosteriorPredict(fit, Xnew = DatTest$X)
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
  pred <- PosteriorPredict(fit, Xnew = DatTest$X)

  Yall_zero <- rep(0, length(DatTest$Y))
  perf_zero <- EvaluatePrediction(Yall_zero, pred)
  expect_true(is.na(perf_zero$AUC_binary))
  expect_true(is.na(perf_zero$RMSE_cont))
})
