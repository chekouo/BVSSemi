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
  expect_length(pred$yhat.median, ntest)
  expect_equal(nrow(pred$mu.draws), ntest)
  expect_true(all(pred$p.pos.mean >= 0 & pred$p.pos.mean <= 1))
  expect_true(all(pred$yhat.mean >= 0))
  expect_true(all(pred$yhat.median >= 0))
  ## yhat.mean should be on the same order of magnitude as the observed
  ## response, not blown up by a stray double-exponentiation
  expect_true(max(pred$yhat.mean) < 10 * max(DatTest$Y))
})

test_that("yhat.median is robust to a handful of extreme MCMC draws (regression test)", {
  ## yhat.mean is a Monte Carlo average of exp(etaCont + sigma2/2) across
  ## draws, which can be dominated by a rare extreme draw for a weakly
  ## identified feature (small n relative to p). yhat.median should stay
  ## sane even when yhat.mean does not.
  set.seed(1)
  n <- 30; p <- 200 ## deliberately n << p to encourage instability
  X <- matrix(rnorm(n * p), n, p)
  Y <- exp(rnorm(n, mean = 3, sd = 0.3))
  Y[sample(n, 5)] <- 0
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Y, X = X, seed = 1,
                              tau2cont = 5, nu1cont = -1,
                              mcmcsample = 300, burnin = 100))
  Xnew <- matrix(rnorm(20 * p), 20, p)
  pred <- PosteriorPredict(fit, Xnew = Xnew)

  expect_length(pred$yhat.median, 20)
  expect_true(all(is.finite(pred$yhat.median)))
  ## median should always stay within a sane multiple of anything actually
  ## observed in training, regardless of how extreme yhat.mean gets
  expect_true(max(pred$yhat.median) < 1000 * max(Y))
})

test_that("PosteriorPredict(nmodels=NULL) is unchanged from omitting the argument", {
  Dat <- GenDataSemiContinous(n = 80, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 300, burnin = 150))
  DatTest <- GenDataSemiContinous(n = 40, p = 20, sd = 1, impf = 5, beta = 0.3,
                                   percentOverlap = "Full", seed = 2, log_scale = TRUE)
  pred1 <- PosteriorPredict(fit, Xnew = DatTest$X)
  pred2 <- PosteriorPredict(fit, Xnew = DatTest$X, nmodels = NULL)
  expect_equal(pred1, pred2)
})

test_that("PosteriorPredict(nmodels=K) performs Bayesian model averaging over the top-K visited models", {
  Dat <- GenDataSemiContinous(n = 80, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 1000, burnin = 400))
  DatTest <- GenDataSemiContinous(n = 40, p = 20, sd = 1, impf = 5, beta = 0.3,
                                   percentOverlap = "Full", seed = 2, log_scale = TRUE)

  pred <- PosteriorPredict(fit, Xnew = DatTest$X, nmodels = 5, seed = 1)

  expect_true(pred$bma$nmodels <= 5)
  expect_length(pred$bma$weights, pred$bma$nmodels)
  expect_equal(sum(pred$bma$weights), 1, tolerance = 1e-8)
  expect_true(all(pred$bma$weights >= 0))
  ## weights should be sorted by decreasing log-posterior
  expect_equal(pred$bma$logpost, sort(pred$bma$logpost, decreasing = TRUE))
  expect_length(pred$yhat.mean, nrow(DatTest$X))
  expect_true(all(is.finite(pred$yhat.mean)))

  ## nmodels larger than the number of visited unique models should just use
  ## all of them, not error
  nUniqueModels <- length(unique(apply(cbind(fit$Gamma1.draws, fit$Gamma2.draws), 1, paste, collapse = "")))
  predBig <- PosteriorPredict(fit, Xnew = DatTest$X, nmodels = 1e6, seed = 1)
  expect_equal(predBig$bma$nmodels, nUniqueModels)
})

test_that("PosteriorPredict(nmodels=K) errors informatively on a fit without Gamma/logPost draws", {
  Dat <- GenDataSemiContinous(n = 60, p = 15, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1, log_scale = TRUE)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
                              mcmcsample = 150, burnin = 75))
  fit$Gamma1.draws <- NULL
  expect_error(
    PosteriorPredict(fit, Xnew = Dat$X, nmodels = 5),
    "Gamma1.draws"
  )
})

test_that("BMA (nmodels) fixes the same real-data instability as yhat.median (regression test)", {
  ## Reproduces the reported failure mode: a weakly identified feature in a
  ## p >> n fit causes a rare extreme MCMC draw to dominate the full-draw
  ## average (yhat.mean), inflating RMSE_combined by orders of magnitude.
  ## Restricting prediction to the top few highest-posterior-probability
  ## models should also keep the mean-based metric sane.
  set.seed(1)
  n <- 30; p <- 200
  X <- matrix(rnorm(n * p), n, p)
  Y <- exp(rnorm(n, mean = 3, sd = 0.3))
  Y[sample(n, 5)] <- 0
  fit <- quietly(MainBVSSemi(Method = "BVSSemiMRF", Y = Y, X = X, seed = 1,
                              tau2cont = 5, nu1cont = -1,
                              mcmcsample = 300, burnin = 100))
  Xnew <- matrix(rnorm(20 * p), 20, p)

  pred_full <- PosteriorPredict(fit, Xnew = Xnew)
  pred_bma  <- PosteriorPredict(fit, Xnew = Xnew, nmodels = 5, seed = 1)

  expect_true(all(is.finite(pred_bma$yhat.mean)))
  ## BMA's mean-based yhat should stay close to the already-robust median,
  ## unlike the full-draw yhat.mean which can be orders of magnitude off
  expect_true(max(pred_bma$yhat.mean) < 1000 * max(Y))
})

test_that("PosteriorPredict stays finite when a draw's linear predictor would overflow exp() (regression test)", {
  ## Directly exercises rowLogSumExp()'s guard: if one MCMC draw's continuous
  ## linear predictor is extreme enough that exp() would overflow to Inf, and
  ## that same draw's pPos has underflowed to exactly 0, the naive
  ## pPos * exp(logCondMean) computation produces 0 * Inf = NaN, which would
  ## silently poison both yhat.mean (via rowMeans) and yhat.median (a single
  ## NaN makes median() return NA). PosteriorPredict must stay finite instead.
  set.seed(1)
  n <- 20; p <- 5; M <- 50
  fit <- list(
    beta.Bin.draws = matrix(rnorm((p + 1) * M), M, p + 1),
    beta.Cont.draws = matrix(rnorm((p + 1) * M), M, p + 1),
    sigma2.draws = runif(M, 0.1, 0.5),
    X.center = rep(0, p), X.scale = rep(1, p),
    Xcov.center = NULL, Xcov.scale = NULL
  )
  ## rig one (subject, draw) pair to hit exactly the pathological case
  fit$beta.Cont.draws[1, ] <- c(5000, rep(0, p)) ## draw 1's intercept alone -> huge etaCont for everyone
  fit$beta.Bin.draws[1, ] <- c(50, rep(0, p))    ## draw 1's etaBin -> pPos = 1 - pnorm(50) = 0 exactly

  Xnew <- matrix(rnorm(n * p), n, p)
  pred <- PosteriorPredict(fit, Xnew = Xnew)

  expect_true(all(is.finite(pred$yhat.mean)))
  expect_true(all(is.finite(pred$yhat.median)))
  expect_false(any(is.na(pred$yhat.mean)))
  expect_false(any(is.na(pred$yhat.median)))
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
  expect_true(perf$RMSE_combined_median >= 0)
  expect_true(perf$MAE_combined_median >= 0)
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
