#' Evaluate Posterior Predictive Performance for a Semicontinuous Response
#'
#' Computes performance metrics comparing a held-out semicontinuous response
#' \code{Ytest} to the posterior predictive summaries returned by
#' \code{\link{PosteriorPredict}}: discrimination and calibration metrics
#' for the binary (zero/nonzero) part, error and coverage metrics for the
#' continuous (magnitude) part conditional on the response being nonzero,
#' and overall error metrics on the combined (occurrence x magnitude)
#' prediction.
#'
#' @param Ytest A numeric vector of held-out semicontinuous response values
#'   (a mix of exact zeros and continuous positive values), on their
#'   original (unstandardized, un-logged) scale.
#' @param pred A list returned by \code{\link{PosteriorPredict}}, containing
#'   at least \code{p.pos.mean}, \code{ycont.mean}, \code{ycont.draws},
#'   \code{sigma2.draws} and \code{yhat.mean} for the same subjects as
#'   \code{Ytest}.
#' @param coverage_level Nominal coverage level of the posterior predictive
#'   interval used to assess calibration of the continuous part (e.g.
#'   \code{0.95} for a 95 percent interval). Defaults to \code{0.95}.
#' @param log_scale Logical. Must match the \code{log_scale} value used to
#'   fit the model with \code{\link{MainBVSSemi}} and to generate
#'   \code{pred} with \code{\link{PosteriorPredict}}. If \code{TRUE} (the
#'   default), \code{pred$ycont.mean} and \code{pred$ycont.draws} are on the
#'   \code{log(Y)} scale, so the nonzero values of \code{Ytest} are
#'   log-transformed before being compared to them (and must therefore be
#'   strictly positive); if \code{FALSE}, \code{Ytest} is compared to them
#'   on its original scale. \code{RMSE_combined}/\code{MAE_combined} are
#'   unaffected, since \code{pred$yhat.mean} is always on the original
#'   scale of \code{Ytest}.
#'
#' @return A named list of performance metrics:
#'   \describe{
#'     \item{\code{AUC_binary}}{Area under the ROC curve for predicting
#'       whether the response is nonzero, using \code{pred$p.pos.mean} as
#'       the score. \code{NA} when \code{Ytest} is all zero or all
#'       nonzero.}
#'     \item{\code{Brier_binary}}{Brier score (mean squared error) of
#'       \code{pred$p.pos.mean} against the nonzero indicator.}
#'     \item{\code{RMSE_cont}, \code{MAE_cont}}{Root-mean-squared and mean
#'       absolute error of \code{pred$ycont.mean} against the (nonzero, and
#'       optionally log-transformed) observations, restricted to subjects
#'       with \code{Ytest != 0}. \code{NA} when fewer than 2 such subjects
#'       are present.}
#'     \item{\code{Corr_cont}}{Pearson correlation between
#'       \code{pred$ycont.mean} and the observations among subjects with
#'       \code{Ytest != 0}.}
#'     \item{\code{Coverage<XX>_cont}}{Empirical coverage of the
#'       \code{coverage_level} posterior predictive interval for the
#'       continuous part, among subjects with \code{Ytest != 0}, where
#'       \code{<XX>} is \code{coverage_level} expressed as a percentage
#'       (e.g. \code{Coverage95_cont}).}
#'     \item{\code{RMSE_combined}, \code{MAE_combined}}{Root-mean-squared
#'       and mean absolute error of \code{pred$yhat.mean} against
#'       \code{Ytest}, over all subjects, on the original response scale.}
#'   }
#'
#' @seealso \code{\link{MainBVSSemi}}, \code{\link{PosteriorPredict}},
#'   \code{\link{CVPredictBVSSemi}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Dat <- GenDataSemiContinous(n = 500, p = 500, sd = 1, impf = 20,
#'                              beta = 0.3, percentOverlap = "Full", seed = 1,
#'                              log_scale = TRUE)
#' fit <- MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X,
#'                     mcmcsample = 10000, burnin = 5000)
#'
#' DatTest <- GenDataSemiContinous(n = 100, p = 500, sd = 1, impf = 20,
#'                                  beta = 0.3, percentOverlap = "Full", seed = 2,
#'                                  log_scale = TRUE)
#' pred <- PosteriorPredict(fit, Xnew = DatTest$X, Xcovnew = NULL)
#'
#' perf <- EvaluatePrediction(DatTest$Y, pred)
#' str(perf)
#' }
## ---- prediction performance metrics ----
EvaluatePrediction <- function(Ytest, pred, coverage_level = 0.95, log_scale = TRUE) {
  pos <- Ytest != 0
  ind <- as.numeric(pos)
  out <- list()
  if (length(unique(ind)) > 1) {
    ## simple trapezoidal AUC (avoids an extra package dependency)
    ord <- order(pred$p.pos.mean)
    ind_o <- ind[ord]; scr_o <- pred$p.pos.mean[ord]
    n1 <- sum(ind == 1); n0 <- sum(ind == 0)
    ranks <- rank(pred$p.pos.mean)
    out$AUC_binary <- (sum(ranks[ind == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  } else out$AUC_binary <- NA
  out$Brier_binary <- mean((pred$p.pos.mean - ind)^2)

  if (sum(pos) >= 2) {
    ## pred$ycont.mean/ycont.draws are on the scale the continuous model was
    ## fit on (log(Y) when log_scale = TRUE), so obs must be put on that
    ## same scale before comparing, or RMSE_cont/MAE_cont mix scales
    if (log_scale) {
      if (any(Ytest[pos] <= 0)) {
        stop("log_scale = TRUE requires all nonzero values of Ytest to be positive")
      }
      obs <- log(Ytest[pos])
    } else {
      obs <- Ytest[pos]
    }
    muhat <- pred$ycont.mean[pos]
    resid <- obs - muhat
    out$RMSE_cont <- sqrt(mean(resid^2))
    out$MAE_cont  <- mean(abs(resid))
    out$Corr_cont <- if (sd(muhat) > 0) cor(obs, muhat) else NA
    M <- ncol(pred$ycont.draws)
    eps <- matrix(rnorm(sum(pos) * M), sum(pos), M)
    predDraws <- pred$ycont.draws[pos, , drop = FALSE] + sweep(eps, 2, sqrt(pred$sigma2.draws), "*")
    lo <- apply(predDraws, 1, quantile, probs = (1 - coverage_level) / 2)
    hi <- apply(predDraws, 1, quantile, probs = 1 - (1 - coverage_level) / 2)
    out[[paste0("Coverage", round(coverage_level * 100), "_cont")]] <- mean(obs >= lo & obs <= hi)
  } else {
    out$RMSE_cont <- out$MAE_cont <- out$Corr_cont <- NA
  }
  residAll <- Ytest - pred$yhat.mean
  out$RMSE_combined <- sqrt(mean(residAll^2))
  out$MAE_combined  <- mean(abs(residAll))
  out
}