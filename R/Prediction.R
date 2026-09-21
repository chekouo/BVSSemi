#' Posterior Predictive Distribution for New Subjects via Bayesian Model Averaging
#'
#' Computes posterior predictive summaries for new subjects from a fitted
#' BVSSemi semicontinuous regression model (see \code{\link{MainBVSSemi}}),
#' which jointly models the probability that the outcome is nonzero (a
#' probit occurrence model) and its magnitude given it is nonzero (a linear
#' model). Predictions are formed by Bayesian model averaging (BMA) over the
#' \code{nmodels} distinct \code{(Gamma.cont, Gamma.bin)} feature-selection
#' models visited during MCMC with the highest estimated log-posterior
#' (\code{fit$logPost.models}), each represented by its deterministic
#' posterior-mode coefficients (\code{fit$beta.Cont.mode},
#' \code{fit$beta.Bin.mode}, computed once per unique model as
#' post-processing in \code{\link{MainBVSSemi}}, not from any single noisy
#' MCMC draw) and a single shared residual variance (\code{fit$sigma2.mean}),
#' combined by a weighted average using each kept model's normalized
#' posterior probability. \code{Xnew} and \code{Xcovnew} are standardized
#' using the column means and standard deviations stored in \code{fit} (as
#' computed by \code{\link{MainBVSSemi}} on the training data), so they
#' should be passed on their original (unstandardized) scale.
#'
#' @param fit A list returned by \code{\link{MainBVSSemi}} containing at
#'   least the following components:
#'   \describe{
#'     \item{\code{Gamma.cont.models}, \code{Gamma.bin.models}}{Matrices (one
#'       row per unique visited model) of the continuous- and binary-model
#'       inclusion indicators for that model.}
#'     \item{\code{logPost.models}}{Length-\eqn{K} vector (\eqn{K} = number
#'       of unique visited models) of each model's estimated log-posterior.}
#'     \item{\code{beta.Cont.mode}, \code{beta.Bin.mode}}{\eqn{K} x
#'       (1 + pc + p) matrices of each model's posterior-mode coefficients.}
#'     \item{\code{sigma2.mean}}{A single number: the posterior mean residual
#'       variance, shared across every model.}
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
#' @param nmodels Integer, required. Bayesian model averaging is performed
#'   over the top \code{nmodels} visited models by estimated log-posterior
#'   (\code{fit$logPost.models}); if fewer than \code{nmodels} models were
#'   visited, all of them are used.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{p.pos.mean}}{Length-n vector of the BMA-weighted
#'       probability that the outcome is nonzero for each new subject.}
#'     \item{\code{ycont.mean}}{Length-n vector of the BMA-weighted
#'       continuous (magnitude) model linear predictor.}
#'     \item{\code{ycont.draws}}{An n x \code{nmodels} matrix of the
#'       continuous (magnitude) model linear predictor for every kept model
#'       (unweighted; used e.g. by \code{\link{EvaluatePrediction}}'s
#'       posterior-interval coverage check).}
#'     \item{\code{sigma2.draws}}{The single shared \code{fit$sigma2.mean},
#'       repeated once per kept model (length \code{nmodels}) for shape
#'       compatibility with \code{ycont.draws}.}
#'     \item{\code{yhat.mean}}{Length-n vector of the BMA-weighted predicted
#'       outcome, \eqn{E[Y | x] = \sum_k w_k P(Y > 0 | x, model_k) E[Y | Y >
#'       0, x, model_k]}, where \eqn{w_k} are the normalized BMA weights.}
#'     \item{\code{bma}}{A list with \code{nmodels} (models actually kept,
#'       \code{<= nmodels} if fewer were visited), \code{weights} (their
#'       normalized BMA weights, summing to 1) and \code{logpost} (their
#'       estimated log-posterior values), named by a string key identifying
#'       each kept \code{(Gamma.cont, Gamma.bin)} model.}
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
#' ## Predict for the new subjects via BMA over the top 10 visited models
#' pred <- PosteriorPredict(fit, Xnew = DatTest$X, Xcovnew = NULL, nmodels = 10)
#'
#' ## Posterior mean predicted outcome for each new subject
#' str(pred$yhat.mean)
#' }
PosteriorPredict <- function(fit, Xnew, Xcovnew = NULL, log_scale = TRUE, nmodels) {
  if (missing(nmodels)) {
    stop("nmodels is required: PosteriorPredict performs Bayesian model averaging over ",
         "the top nmodels visited models (see fit$logPost.models); there is no ",
         "full-MCMC-draw averaging fallback.")
  }
  if (is.null(fit$Gamma.cont.models) || is.null(fit$Gamma.bin.models) || is.null(fit$logPost.models) ||
      is.null(fit$beta.Cont.mode) || is.null(fit$beta.Bin.mode) || is.null(fit$sigma2.mean)) {
    stop("nmodels requires Gamma.cont.models/Gamma.bin.models/logPost.models/beta.Cont.mode/",
         "beta.Bin.mode/sigma2.mean in fit; refit with the current MainBVSSemi().")
  }

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

  ## Keep the top nmodels visited models by estimated log-posterior, and
  ## turn their log-posteriors into normalized Bayesian model averaging
  ## weights (a softmax, numerically stabilized by subtracting the max).
  K <- min(nmodels, length(fit$logPost.models))
  keep <- order(fit$logPost.models, decreasing = TRUE)[seq_len(K)]
  lp <- fit$logPost.models[keep]
  w <- exp(lp - max(lp))
  w <- w / sum(w)
  keyOf <- function(i) paste(c(fit$Gamma.cont.models[i, ], fit$Gamma.bin.models[i, ]), collapse = "")
  bma <- list(nmodels = K, weights = stats::setNames(w, vapply(keep, keyOf, character(1))),
              logpost = stats::setNames(lp, vapply(keep, keyOf, character(1))))

  XX <- cbind(rep(1, n), Xcovnew, Xnew)
  etaBin  <- XX %*% t(fit$beta.Bin.mode[keep, , drop = FALSE])   # n x K
  etaCont <- XX %*% t(fit$beta.Cont.mode[keep, , drop = FALSE])  # n x K
  pPos <- 1 - pnorm(etaBin)                                       # P(y>0 | x, model k)
  ## A single shared residual variance across all K models (not
  ## model-specific), repeated for shape compatibility with etaCont/pPos.
  sigma2 <- rep(fit$sigma2.mean, K)

  if (log_scale) {
    ## Work in log-space (log(weight) + log(pPos) + logCondMean) for as long
    ## as possible, rather than forming weight * pPos * exp(logCondMean)
    ## directly: if a model's linear predictor is extreme enough that exp()
    ## would overflow to Inf, multiplying it by a same-model pPos that has
    ## underflowed to exactly 0 produces 0 * Inf = NaN, which would silently
    ## poison yhat.mean. Staying in log-space until the final exponentiation
    ## makes that case degenerate to 0 (a near-zero occurrence probability
    ## wins over an unbounded magnitude) instead of NaN.
    logCondMean <- sweep(etaCont, 2, sigma2 / 2, "+")
    logPW <- log(pPos) + logCondMean
    logPW[is.nan(logPW)] <- -Inf ## 0 * Inf coincidence: treat as no contribution
    logW <- sweep(logPW, 2, log(w), "+")
    yhat.mean <- exp(rowLogSumExp(logW)) ## weights already sum to 1
  } else {
    yhat.mean <- as.vector((pPos * etaCont) %*% w)
  }

  list(p.pos.mean = as.vector(pPos %*% w), ycont.mean = as.vector(etaCont %*% w),
       ycont.draws = etaCont, sigma2.draws = sigma2, yhat.mean = yhat.mean, bma = bma)
}

### Numerically stable row-wise log(sum(exp(.))), robust to rows that are
### entirely -Inf (which would otherwise produce NaN from -Inf - (-Inf)).
rowLogSumExp <- function(M) {
  rmax <- apply(M, 1, max)
  out <- rep(-Inf, nrow(M))
  ok <- is.finite(rmax)
  if (any(ok)) {
    out[ok] <- rmax[ok] + log(rowSums(exp(M[ok, , drop = FALSE] - rmax[ok])))
  }
  out
}
