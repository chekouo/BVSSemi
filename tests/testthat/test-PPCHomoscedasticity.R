## Semicontinuous data whose nonzero log-values have a given residual-SD
## profile: sd_fun(eta) maps the true linear predictor to the residual SD.
make_hetero_data <- function(sd_fun, n = 300, p = 10, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  eta <- 1.5 * X[, 1] + 1.5 * X[, 2]
  logY <- eta + rnorm(n, sd = sd_fun(eta))
  nonzero <- rbinom(n, 1, 0.6) == 1
  list(Y = ifelse(nonzero, exp(logY), 0), X = X)
}

test_that("PPCHomoscedasticity returns well-formed output", {
  Dat <- make_hetero_data(function(eta) rep(1, length(eta)), n = 120)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiIndep", Y = Dat$Y, X = Dat$X,
                              seed = 1, mcmcsample = 400, burnin = 200, thin = 2))
  seedBefore <- .Random.seed
  hc <- PPCHomoscedasticity(fit, plot = FALSE)
  expect_identical(.Random.seed, seedBefore) ## caller's RNG state untouched
  n2 <- sum(Dat$Y != 0)
  expect_equal(length(hc$T.obs) + hc$n.draws.dropped, nrow(fit$fitted.cont.draws))
  expect_length(hc$T.rep, length(hc$T.obs))
  expect_true(all(abs(c(hc$T.obs, hc$T.rep)) <= 1))
  expect_equal(hc$p.value, mean(hc$T.rep >= hc$T.obs))
  for (v in c("fitted.mean", "resid.mean", "resid.lower", "resid.upper")) {
    expect_length(hc[[v]], n2)
  }
  expect_true(all(hc$resid.lower <= hc$resid.mean & hc$resid.mean <= hc$resid.upper))
  expect_equal(hc$resid.mean, fit$y.cont - hc$fitted.mean)
  ## same seed, same result; plotting does not change the numbers
  pdf(NULL); on.exit(dev.off())
  expect_equal(PPCHomoscedasticity(fit, plot = TRUE), hc)
})

test_that("PPCHomoscedasticity flags a funnel but not constant variance", {
  homo <- make_hetero_data(function(eta) rep(1, length(eta)))
  hetero <- make_hetero_data(function(eta) 0.2 * exp(0.6 * eta))
  run <- function(D) {
    fit <- quietly(MainBVSSemi(Method = "BVSSemiIndep", Y = D$Y, X = D$X,
                                seed = 1, mcmcsample = 600, burnin = 300, thin = 2))
    PPCHomoscedasticity(fit, plot = FALSE)$p.value
  }
  expect_gt(run(homo), 0.05)
  expect_lt(run(hetero), 0.05) ## |residual| grows with fitted -> PPP near 0
})

test_that("PPCHomoscedasticity rejects fits without the stored draws", {
  expect_error(PPCHomoscedasticity(list(sigma2.draws = 1)), "fitted.cont.draws")
})
