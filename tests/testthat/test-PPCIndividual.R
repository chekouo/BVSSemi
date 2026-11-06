make_ind_data <- function(n = 150, p = 10, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  logY <- 1.5 * X[, 1] + 1.5 * X[, 2] + rnorm(n)
  nonzero <- rbinom(n, 1, 0.6) == 1
  list(Y = ifelse(nonzero, exp(logY), 0), X = X)
}

test_that("PPCIndividual returns well-formed p-values matching their exact value", {
  Dat <- make_ind_data()
  fit <- quietly(MainBVSSemi(Method = "BVSSemiIndep", Y = Dat$Y, X = Dat$X,
                              seed = 1, mcmcsample = 2200, burnin = 200, thin = 1))
  expect_equal(fit$index.cont, which(Dat$Y != 0))
  seedBefore <- .Random.seed
  ind <- PPCIndividual(fit, plot = FALSE)
  expect_identical(.Random.seed, seedBefore) ## caller's RNG state untouched
  n2 <- sum(Dat$Y != 0)
  expect_length(ind$p.value, n2)
  expect_equal(names(ind$p.value), as.character(which(Dat$Y != 0)))
  expect_true(all(ind$p.value >= 0 & ind$p.value <= 1))
  ## Monte Carlo p_i vs exact mean_s P(y_rep >= y_i | mu^(s), sigma^(s))
  z <- sweep(-fit$fitted.cont.draws, 2, fit$y.cont, "+") / sqrt(fit$sigma2.draws)
  exact <- colMeans(pnorm(z, lower.tail = FALSE))
  expect_lt(max(abs(ind$p.value - exact)), 0.05)
  expect_equal(ind$underpredicted, ind$index[ind$p.value < 0.05])
  expect_equal(ind$overpredicted, ind$index[ind$p.value > 0.95])
  expect_equal(ind$prop.extreme, mean(ind$p.value < 0.05 | ind$p.value > 0.95))
  pdf(NULL); on.exit(dev.off())
  expect_equal(PPCIndividual(fit, plot = TRUE), ind) ## plotting changes nothing
})

test_that("PPCIndividual flags a planted outlier as underpredicted", {
  Dat <- make_ind_data()
  i <- which(Dat$Y != 0)[5]
  Dat$Y[i] <- Dat$Y[i] * exp(6) ## +6 on the log scale, ~6 residual SDs
  fit <- quietly(MainBVSSemi(Method = "BVSSemiIndep", Y = Dat$Y, X = Dat$X,
                              seed = 1, mcmcsample = 600, burnin = 300, thin = 2))
  ind <- PPCIndividual(fit, plot = FALSE)
  expect_true(i %in% ind$underpredicted)
  expect_lt(ind$p.value[[as.character(i)]], 0.05)
})

test_that("PPCIndividual validates its inputs", {
  expect_error(PPCIndividual(list(sigma2.draws = 1)), "index.cont")
  Dat <- make_ind_data(n = 60)
  fit <- quietly(MainBVSSemi(Method = "BVSSemiIndep", Y = Dat$Y, X = Dat$X,
                              seed = 1, mcmcsample = 100, burnin = 50, thin = 5))
  expect_error(PPCIndividual(fit, lower = 0.9, upper = 0.1, plot = FALSE), "lower")
})
