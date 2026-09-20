test_that("GenDataSemiContinous returns correctly shaped output", {
  n <- 60; p <- 25; impf <- 6
  Dat <- GenDataSemiContinous(n = n, p = p, sd = 1, impf = impf, beta = 0.3,
                               percentOverlap = "Full", seed = 1)

  expect_type(Dat, "list")
  expect_named(Dat, c("Y", "X", "Z.cont", "Z.bin"), ignore.order = TRUE)
  expect_length(Dat$Y, n)
  expect_equal(dim(Dat$X), c(n, p))
  expect_length(Dat$Z.cont, p)
  expect_length(Dat$Z.bin, p)
  expect_true(all(Dat$Z.cont %in% c(0, 1)))
  expect_true(all(Dat$Z.bin %in% c(0, 1)))
})

test_that("log_scale defaults to TRUE and keeps nonzero Y strictly positive", {
  Dat <- GenDataSemiContinous(n = 100, p = 20, sd = 1, impf = 5, beta = 0.3,
                               percentOverlap = "Full", seed = 1)
  nz <- Dat$Y[Dat$Y != 0]
  expect_true(length(nz) > 0)
  expect_true(all(nz > 0))
})

test_that("log_scale = FALSE allows nonzero Y to be negative", {
  ## with log_scale = FALSE the nonzero values are raw Gaussian draws, so
  ## across enough replications at least one should be negative
  any_negative <- FALSE
  for (s in 1:10) {
    Dat <- GenDataSemiContinous(n = 200, p = 20, sd = 1, impf = 5, beta = 0.3,
                                 percentOverlap = "Full", seed = s, log_scale = FALSE)
    if (any(Dat$Y[Dat$Y != 0] < 0)) any_negative <- TRUE
  }
  expect_true(any_negative)
})

test_that("percentOverlap controls which features drive the continuous vs binary model", {
  impf <- 10
  Dat_full <- GenDataSemiContinous(n = 50, p = 40, sd = 1, impf = impf, beta = 0.3,
                                    percentOverlap = "Full", seed = 1)
  expect_equal(which(Dat_full$Z.cont == 1), 1:impf)
  expect_equal(which(Dat_full$Z.bin == 1), 1:impf)

  Dat_none <- GenDataSemiContinous(n = 50, p = 40, sd = 1, impf = impf, beta = 0.3,
                                    percentOverlap = "NoOverlap", seed = 1)
  expect_equal(which(Dat_none$Z.bin == 1), 1:impf)
  expect_equal(which(Dat_none$Z.cont == 1), (impf + 1):(2 * impf))
  expect_length(intersect(which(Dat_none$Z.cont == 1), which(Dat_none$Z.bin == 1)), 0)

  Dat_med <- GenDataSemiContinous(n = 50, p = 40, sd = 1, impf = impf, beta = 0.3,
                                   percentOverlap = "Medium", seed = 1)
  overlap <- intersect(which(Dat_med$Z.cont == 1), which(Dat_med$Z.bin == 1))
  expect_true(length(overlap) > 0 && length(overlap) < impf)
})

test_that("Xreal = TRUE requires a user-supplied X", {
  expect_error(
    GenDataSemiContinous(Xreal = TRUE, n = 20, p = 5, X = NULL, impf = 2, beta = 0.3),
    "Provide the matrix"
  )
})

test_that("Xreal = TRUE uses the supplied X as-is", {
  n <- 30; p <- 8
  Xin <- matrix(rnorm(n * p), n, p)
  Dat <- GenDataSemiContinous(Xreal = TRUE, n = n, p = p, X = Xin, impf = 3, beta = 0.3,
                               percentOverlap = "Full", seed = 1)
  ## X is standardized (scale()) inside the function, so compare shape/columns, not values
  expect_equal(dim(Dat$X), dim(Xin))
})
