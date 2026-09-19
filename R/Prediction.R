#' Posterior Predictive Distribution for New Subjects
#'
#' Computes posterior predictive summaries for new subjects from a fitted
#' BVSSemi semicontinuous regression model (see \code{\link{MainBVSSemi}}),
#' which jointly models the probability that the outcome is nonzero (a
#' probit occurrence model) and its magnitude given it is nonzero (a linear
#' model). Predictions are obtained by averaging over the retained
#' post-burn-in MCMC draws of the regression coefficients and residual
#' variance. \code{Xnew} and \code{Xcovnew} are standardized using the
#' column means and standard deviations stored in \code{fit} (as computed
#' by \code{\link{MainBVSSemi}} on the training data), so they should be
#' passed on their original (unstandardized) scale.
#'
#' @param fit A list returned by \code{\link{MainBVSSemi}} (or an equivalent
#'   BVSSemi fitting function) containing at least the following components:
#'   \describe{
#'     \item{\code{beta.Bin.draws}}{An M x (1 + pc + p) matrix of posterior
#'       draws of the occurrence-model coefficients (intercept, forced-in
#'       covariates: pc, then predictors:p), one row per retained MCMC
#'       iteration.}
#'     \item{\code{beta.Cont.draws}}{An M x (1 + pc + p) matrix of posterior
#'       draws of the continuous (magnitude) model coefficients, with the
#'       same column layout as \code{beta.Bin.draws}.}
#'     \item{\code{sigma2.draws}}{A numeric vector of length M with the
#'       posterior draws of the residual variance of the continuous model.}
#'     \item{\code{X.center}, \code{X.scale}}{Length-p vectors of the column
#'       means and standard deviations used to standardize \code{X} when
#'       fitting the model, used to standardize \code{Xnew} the same way.}
#'     \item{\code{Xcov.center}, \code{Xcov.scale}}{Length-pc vectors of the
#'       column means and standard deviations used to standardize
#'       \code{Xcov} when fitting the model (\code{NULL} when \code{Xcov}
#'       was not used), used to standardize \code{Xcovnew} the same way.}
#'   }
#' @param Xnew A numeric matrix (or an object coercible via
#'   \code{as.matrix}) of dimension n x p containing the covariates, on
#'   their original (unstandardized) scale, for the n new subjects, with
#'   columns in the same order used to fit the model.
#' @param Xcovnew An optional numeric matrix or vector of dimension
#'   n x pc containing forced-in covariates, on their original
#'   (unstandardized) scale, for the new subjects, matching the
#'   \code{Xcov} used when fitting the model. Defaults to \code{NULL}
#'   when no such covariates were used.
#' @param log_scale Logical. If \code{TRUE} (the default, matching the
#'   default of \code{\link{MainBVSSemi}}), the continuous part of the
#'   model is assumed to have been fit on the log scale, and the
#'   conditional mean on the original scale is obtained via a log-normal
#'   bias correction, \code{exp(mu + sigma2 / 2)}. If \code{FALSE}, the
#'   continuous linear predictor is used directly as the conditional mean.
#'   Must match the \code{log_scale} value used when fitting \code{fit}
#'   with \code{\link{MainBVSSemi}}.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{p.pos.mean}}{Length-n vector of posterior mean
#'       probabilities that the outcome is nonzero for each new subject.}
#'     \item{\code{mu.mean}}{Length-n vector of posterior mean linear
#'       predictors from the continuous (magnitude) model.}
#'     \item{\code{mu.draws}}{An n x M matrix of the continuous linear
#'       predictor for every posterior draw and new subject.}
#'     \item{\code{sigma2.draws}}{The length-M vector of residual variance
#'       draws, passed through from \code{fit}.}
#'     \item{\code{yhat.mean}}{Length-n vector of posterior mean predicted
#'       outcomes, E[Y | x] = P(Y > 0 | x) * E[Y | Y > 0, x], averaged over
#'       draws.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Simulate a training set with GenDataSemiContinous(): a semicontinuous
#' ## response Y (n = 500) driven by p = 500 features, of which 20 are truly
#' ## important to both the continuous and binary parts of the model
#' ## (percentOverlap = "Full"), with no forced-in covariates (Xcov = NULL).
#' ## log_scale = TRUE exponentiates the nonzero values of Y so they are
#' ## strictly positive, matching MainBVSSemi's default log_scale = TRUE
#' Dat <- GenDataSemiContinous(n = 500, p = 500, sd = 1, impf = 20,
#'                              beta = 0.3, percentOverlap = "Full", seed = 1,
#'                              log_scale = TRUE)
#'
#' ## Fit a BVSSemi model on the training data
#' fit <- MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X,
#'                     mcmcsample = 10000, burnin = 5000)
#'
#' ## Simulate an independent test set of new subjects from the same
#' ## data-generating model (same true important features, new draws of X)
#' DatTest <- GenDataSemiContinous(n = 100, p = 500, sd = 1, impf = 20,
#'                                  beta = 0.3, percentOverlap = "Full", seed = 2,
#'                                  log_scale = TRUE)
#'
#' ## Predict for the new subjects
#' pred <- PosteriorPredict(fit, Xnew = DatTest$X, Xcovnew = NULL)
#'
#' ## Posterior mean predicted outcome for each new subject
#' str(pred$yhat.mean)
#' }
PosteriorPredict <- function(fit, Xnew, Xcovnew = NULL, log_scale = TRUE) {
  Xnew <- as.matrix(Xnew)
  n <- nrow(Xnew)
  Xnew <- sweep(sweep(Xnew, 2, fit$X.center, "-"), 2, fit$X.scale, "/")
  if (!is.null(Xcovnew)) {
    if (is.vector(Xcovnew)) {
      Xcovnew <- matrix(Xcovnew, n, 1)
    } else {
      Xcovnew <- as.matrix(Xcovnew)
    }
    Xcovnew <- sweep(sweep(Xcovnew, 2, fit$Xcov.center, "-"), 2,
                     fit$Xcov.scale, "/")
  }
  XX <- cbind(rep(1, n), Xcovnew, Xnew)
  etaBin  <- XX %*% t(fit$beta.Bin.draws)     # n x M
  etaCont <- XX %*% t(fit$beta.Cont.draws)    # n x M
  pPos <- 1 - pnorm(etaBin)                    # P(y>0 | x, draw m)
  sigma2 <- fit$sigma2.draws                    # length M
 
  if (log_scale) {
    condMean <- exp(sweep(etaCont, 2, sigma2 / 2, "+"))
  } else {
    condMean <- etaCont
  }
  yhatDraws <- pPos * condMean
 
  list(p.pos.mean = rowMeans(pPos), mu.mean = rowMeans(etaCont), mu.draws = etaCont,
       sigma2.draws = sigma2, yhat.mean = rowMeans(yhatDraws))
}


 
