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
#' @param nmodels Optional integer. If supplied, switches from averaging
#'   over every retained MCMC draw to Bayesian model averaging (BMA) over
#'   only the \code{nmodels} distinct \code{(Gamma1, Gamma2)} models visited
#'   during MCMC with the highest estimated log-posterior (via
#'   \code{fit$logPost.draws}, requires a \code{fit} from the current
#'   \code{\link{MainBVSSemi}}). This avoids the plain full-draw average
#'   being dominated by a rare, poorly identified draw (see
#'   \code{yhat.mean}), by excluding low-posterior-probability models
#'   entirely rather than just being robust to their influence like
#'   \code{yhat.median} is. Defaults to \code{NULL} (use every retained
#'   draw, the original behavior).
#' @param nresample Only used when \code{nmodels} is supplied: the total
#'   number of pseudo-draws to build from the \code{nmodels} kept models,
#'   allocated across them proportional to their normalized BMA weight and
#'   resampled (with replacement) from each model's own visited MCMC draws
#'   (preserving within-model parameter uncertainty, not just between-model
#'   uncertainty). Defaults to \code{2000}.
#' @param seed Only used when \code{nmodels} is supplied: optional integer
#'   seed for the resampling above, for reproducibility.
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
#'       draws. When \code{log_scale = TRUE}, this is a Monte Carlo average
#'       of \code{exp(etaCont + sigma2 / 2)} across draws, which can be
#'       dominated by rare, extreme draws for weakly identified features
#'       (small training \code{n} relative to \code{p}); see
#'       \code{yhat.median}.}
#'     \item{\code{yhat.median}}{Length-n vector of posterior median
#'       predicted outcomes (median across draws of
#'       \code{P(Y > 0 | x, draw) * E[Y | Y > 0, x, draw]}), a more robust
#'       point estimate than \code{yhat.mean} when a handful of draws are
#'       extreme.}
#'     \item{\code{bma}}{Only present when \code{nmodels} was supplied: a
#'       list with \code{nmodels} (models actually kept, \code{<= nmodels}
#'       if fewer were visited), \code{weights} (their normalized BMA
#'       weights) and \code{logpost} (their estimated log-posterior values),
#'       named by a string key identifying each kept
#'       \code{(Gamma1, Gamma2)} model.}
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
PosteriorPredict <- function(fit, Xnew, Xcovnew = NULL, log_scale = TRUE,
                              nmodels = NULL, nresample = 2000, seed = NULL) {
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

  bma <- NULL
  if (!is.null(nmodels)) {
    if (is.null(fit$Gamma1.draws) || is.null(fit$Gamma2.draws) || is.null(fit$logPost.draws)) {
      stop("nmodels requires Gamma1.draws/Gamma2.draws/logPost.draws in fit; refit with the current MainBVSSemi().")
    }
    if (!is.null(seed)) set.seed(seed)

    ## Group retained MCMC draws by their (Gamma1, Gamma2) model, and
    ## estimate each visited model's log-posterior as the average of its
    ## per-visit logPost.draws (Rao-Blackwellizing away the extra noise from
    ## nuisance variables -- e.g. the binary submodel's latent U, resampled
    ## every iteration -- that the per-visit value also reflects).
    keys <- apply(cbind(fit$Gamma1.draws, fit$Gamma2.draws), 1, paste, collapse = "")
    modelIdx <- split(seq_along(keys), keys)
    logPostModel <- vapply(modelIdx, function(idx) mean(fit$logPost.draws[idx]), numeric(1))

    K <- min(nmodels, length(modelIdx))
    keepNames <- names(sort(logPostModel, decreasing = TRUE))[seq_len(K)]
    lp <- logPostModel[keepNames]
    w <- exp(lp - max(lp))
    w <- w / sum(w) ## normalized Bayesian model averaging weights (softmax of log-posterior)

    ## Build nresample pseudo-draws: allocate them across the K kept models
    ## proportional to w (Multinomial), then within each chosen model
    ## resample (with replacement) from the actual MCMC draws that visited
    ## it. This preserves both between-model (BMA) and within-model
    ## (parameter) uncertainty, and lets the unchanged code below treat these
    ## pseudo-draws exactly like ordinary MCMC draws.
    modelCounts <- as.vector(rmultinom(1, nresample, w))
    resampleIdx <- function(idx, cnt) {
      if (cnt == 0) return(integer(0))
      if (length(idx) == 1) return(rep(idx, cnt)) ## avoid sample()'s length-1 "sample from 1:idx" gotcha
      sample(idx, cnt, replace = TRUE)
    }
    sel <- unlist(Map(resampleIdx, modelIdx[keepNames], modelCounts), use.names = FALSE)

    fit <- list(beta.Bin.draws = fit$beta.Bin.draws[sel, , drop = FALSE],
                beta.Cont.draws = fit$beta.Cont.draws[sel, , drop = FALSE],
                sigma2.draws = fit$sigma2.draws[sel])
    bma <- list(nmodels = K, weights = w, logpost = lp)
  }

  XX <- cbind(rep(1, n), Xcovnew, Xnew)
  etaBin  <- XX %*% t(fit$beta.Bin.draws)     # n x M
  etaCont <- XX %*% t(fit$beta.Cont.draws)    # n x M
  pPos <- 1 - pnorm(etaBin)                    # P(y>0 | x, draw m)
  sigma2 <- fit$sigma2.draws                    # length M
 
  if (log_scale) {
    ## Work in log-space (log(pPos) + logCondMean) for as long as possible,
    ## rather than forming pPos * exp(logCondMean) directly: if a rare draw's
    ## linear predictor is extreme enough that exp() would overflow to Inf,
    ## multiplying it by a same-draw pPos that has underflowed to exactly 0
    ## produces 0 * Inf = NaN, which would silently poison both yhat.mean and
    ## yhat.median (a single NaN makes median() return NA). Stay in log-space
    ## until the final exponentiation so that case degenerates to 0 (a
    ## near-zero occurrence probability wins over an unbounded magnitude)
    ## instead of NaN. This does not change results in the ordinary regime
    ## where exp() doesn't overflow (verified to match the naive computation
    ## exactly there) -- it only guards against literal overflow.
    logCondMean <- sweep(etaCont, 2, sigma2 / 2, "+")
    logW <- log(pPos) + logCondMean
    logW[is.nan(logW)] <- -Inf ## 0 * Inf coincidence: treat as no contribution
    yhatDraws <- exp(logW)
    yhat.mean <- exp(rowLogSumExp(logW) - log(ncol(logW)))
  } else {
    condMean <- etaCont
    yhatDraws <- pPos * condMean
    yhat.mean <- rowMeans(yhatDraws)
  }

  ## yhat.mean is a Monte Carlo average of exp(etaCont + sigma2/2) across
  ## draws; when a candidate feature is only weakly identified (small n
  ## relative to p), a handful of draws can have an extreme etaCont for a
  ## given subject and dominate that average (a log-normal-mean/Jensen's-gap
  ## effect), producing wildly inflated point predictions even though most
  ## draws are reasonable. yhat.median is far more robust to this and is
  ## usually the better point estimate in p >> n settings.
  result <- list(p.pos.mean = rowMeans(pPos), mu.mean = rowMeans(etaCont), mu.draws = etaCont,
                 sigma2.draws = sigma2, yhat.mean = yhat.mean,
                 yhat.median = apply(yhatDraws, 1, median))
  if (!is.null(bma)) result$bma <- bma
  result
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


 
