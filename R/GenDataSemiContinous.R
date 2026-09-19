#' Generate Simulated Semicontinuous Data
#'
#' Generates simulated data as described in the reference manuscript: a
#' semicontinuous response \code{Y} (a mix of exact zeros and continuous
#' positive values) driven by a matrix of features \code{X}, where the set
#' of features important to the binary (zero/nonzero) part of the response
#' and the set important to the continuous (magnitude) part can be made to
#' fully overlap, partially overlap, or not overlap at all.
#'
#' @param Xreal Logical. If \code{TRUE}, the matrix of features \code{X}
#'   must be supplied by the user and is used as-is to simulate the
#'   semicontinuous response. If \code{FALSE} (the default), \code{X} is
#'   simulated internally as an \eqn{n \times p} matrix of independent
#'   standard normal entries.
#' @param n Number of individuals (observations) to simulate.
#' @param p Number of features in \code{X}.
#' @param X An optional numeric matrix of features of dimension
#'   \eqn{n \times p}, required when \code{Xreal = TRUE}. Ignored (and
#'   simulated instead) when \code{Xreal = FALSE}.
#' @param sd Standard deviation of the error term added to the continuous
#'   (nonzero) part of the response.
#' @param impf Number of features that are truly important to the binary
#'   part of the response. When \code{percentOverlap = "Full"}, these are
#'   also the features that are important to the continuous part.
#' @param beta The common regression effect size assigned to each important
#'   feature in the data-generating model.
#' @param percentOverlap Character string controlling the amount of overlap
#'   between the sets of features important to the continuous and binary
#'   models: \code{"Full"} (the same \code{impf} features drive both
#'   models), \code{"Medium"} (half of the important features are shared
#'   between the two models), or \code{"NoOverlap"} (the two models have no
#'   important features in common).
#' @param seed Integer seed used to initialize the random number generator.
#' @param log_scale Logical. If \code{TRUE}, the nonzero values of \code{Y}
#'   are exponentiated after being generated, so they are guaranteed to be
#'   strictly positive. Use this to simulate data suitable for fitting with
#'   \code{\link{MainBVSSemi}}'s \code{log_scale = TRUE} option (its
#'   default), which requires all nonzero response values to be positive.
#'   Defaults to \code{FALSE}.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{Y}}{A semicontinuous response vector of length \eqn{n}.}
#'     \item{\code{X}}{The (scaled) matrix of features of dimension
#'       \eqn{n \times p} used to generate \code{Y}.}
#'     \item{\code{Z.cont}}{A binary vector of length \eqn{p} indicating
#'       which features are truly important in the continuous model.}
#'     \item{\code{Z.bin}}{A binary vector of length \eqn{p} indicating
#'       which features are truly important in the binary model.}
#'   }
#'
#' @references
#' Samuel Babatunde, Tolulope Sajobi and Thierry Chekouo (2026),
#' \emph{A Bayesian Variable Selection for Semicontinuous Response data:
#' Application to cardiovascular disease}, submitted.
#'
#' @seealso \code{\link{MainBVSSemi}}
#'
#' @export
#'
#' @examples
#' library(BVSSemi)
#' Dat <- GenDataSemiContinous(n = 500, p = 500, sd = 1, impf = 20,
#'                              beta = 0.3, percentOverlap = "Full", seed = 1,
#'                              log_scale = TRUE)
#' str(Dat)
GenDataSemiContinous <- function(Xreal = FALSE, n = n, p = p, X = NULL, sd = 1, impf = 20, beta = 1,
                                 percentOverlap = "Full", seed = 1, log_scale = TRUE) {
  set.seed(seed)
  # impf is the number of important covariates. if impf=10 then the 10 first features are important
  # percentOverlap is the percentage of overlap important features between model 1 and model 2: "Full", "Medium" and "NoOverlap"
  Gamma1 <- rep(0, p) ### continuous model
  Gamma2 <- rep(0, p)
  Gamma2[1:impf] <- 1 ## binary model
  if (Xreal == T) {
    if (is.null(X)) {
      stop("Provide the matrix of the set of features")
    }
  } else {
    X <- matrix(rnorm(n * p), n, p)
  }
  X <- scale(X)
  Ystar <- X[, 1:impf] %*% c(rep(beta, impf)) + rnorm(n)
  Y <- 1 - (Ystar > 0)
  N1 <- which(Y == 0)
  N2 <- setdiff(1:n, N1)
  error <- sd * rnorm(length(N2))
  if (percentOverlap == "Full") {
    Y[N2] <- X[N2, 1:impf] %*% c(rep(beta, impf)) + error
    Gamma1[1:impf] <- 1
  }
  halffeat <- floor(impf / 2)
  if (percentOverlap == "Medium") {
    Y[N2] <- X[N2, 1:(halffeat)] %*% c(rep(beta, halffeat)) + X[N2, (1 + impf):(impf + impf - halffeat)] %*% c(rep(beta, impf - halffeat)) + error
    Gamma1[1:halffeat] <- 1
    Gamma1[(1 + impf):(impf + impf - halffeat)] <- 1
  }
  if (percentOverlap == "NoOverlap") {
    Y[N2] <- X[N2, (1 + impf):(2 * impf)] %*% c(rep(beta, impf)) + error
    Gamma1[(1 + impf):(2 * impf)] <- 1
  }
  if (log_scale) {
    Y[N2] <- exp(Y[N2])
  }
  return(list(Y = Y, X = X, Z.cont = Gamma1, Z.bin = Gamma2))
}
