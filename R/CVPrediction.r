#' K-Fold Cross-Validated Prediction Performance for BVSSemi
#'
#' Runs K-fold cross-validation of the BVSSemi two-part semicontinuous
#' model: for each fold, fits \code{\link{MainBVSSemi}} on the training
#' subjects, predicts the held-out subjects with
#' \code{\link{PosteriorPredict}}, and scores those predictions with
#' \code{\link{EvaluatePrediction}}. This provides an out-of-sample estimate
#' of predictive performance for a given choice of \code{MainBVSSemi}
#' hyperparameters, passed through \code{...}.
#'
#' @param Y A numeric vector giving the semicontinuous response (a mix of
#'   exact zeros and continuous positive values) for all \eqn{n} subjects.
#' @param X A numeric matrix of candidate features of dimension
#'   \eqn{n \times p}, subject to variable selection.
#' @param Xcov An optional numeric matrix or vector of covariates that are
#'   not subject to variable selection, always included in both models.
#'   Defaults to \code{NULL}.
#' @param K Number of cross-validation folds. Defaults to \code{5}.
#' @param log_scale Logical. Passed through unchanged to
#'   \code{\link{MainBVSSemi}}, \code{\link{PosteriorPredict}} and
#'   \code{\link{EvaluatePrediction}} in every fold, so the continuous
#'   model is fit, predicted, and scored consistently on the same
#'   (\code{log(Y)} when \code{TRUE}) scale. Defaults to \code{TRUE},
#'   matching the default of \code{\link{MainBVSSemi}}; requires all
#'   nonzero values of \code{Y} to be strictly positive.
#' @param seed Integer base seed. Fold \code{k} is fit with seed
#'   \code{seed + k}; the fold assignment itself is also determined by
#'   \code{seed} via \code{set.seed(seed)}.
#' @param ... Additional arguments passed on to \code{\link{MainBVSSemi}}
#'   for every fold (e.g. \code{Method}, \code{mcmcsample}, \code{burnin},
#'   \code{nu1cont}, \code{nu2bin}, etc.).
#'
#' @return A data frame with \code{K} rows (one per fold) and one column
#'   per metric returned by \code{\link{EvaluatePrediction}} (e.g.
#'   \code{AUC_binary}, \code{Brier_binary}, \code{RMSE_cont},
#'   \code{MAE_cont}, \code{Corr_cont}, a \code{Coverage<XX>_cont} column,
#'   \code{RMSE_combined} and \code{MAE_combined}).
#'
#' @seealso \code{\link{MainBVSSemi}}, \code{\link{PosteriorPredict}},
#'   \code{\link{EvaluatePrediction}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Dat <- GenDataSemiContinous(n = 300, p = 200, sd = 1, impf = 20,
#'                              beta = 0.3, percentOverlap = "Full", seed = 1,
#'                              log_scale = TRUE)
#' perf <- CVPredictBVSSemi(Y = Dat$Y, X = Dat$X, K = 5,
#'                           Method = "BVSSemiMRF",
#'                           mcmcsample = 5000, burnin = 1000)
#' colMeans(perf)
#' }
## ---- K-fold cross-validated prediction driver ----
CVPredictBVSSemi <- function(Y, X, Xcov = NULL, K = 5, log_scale = TRUE, seed = 1, ...) {
  set.seed(seed)
  n <- length(Y)
  folds <- sample(rep(1:K, length.out = n))
  results <- vector("list", K)
  for (k in 1:K) {
    te <- which(folds == k); tr <- which(folds != k)
    Xcov_tr <- if (is.null(Xcov)) NULL else Xcov[tr, , drop = FALSE]
    Xcov_te <- if (is.null(Xcov)) NULL else Xcov[te, , drop = FALSE]
    fit <- MainBVSSemi(Y = Y[tr], X = X[tr, , drop = FALSE], Xcov = Xcov_tr, seed = seed + k,
                        log_scale = log_scale, ...)
    pred <- PosteriorPredict(fit, Xnew = X[te, , drop = FALSE], Xcovnew = Xcov_te, log_scale = log_scale)
    results[[k]] <- EvaluatePrediction(Y[te], pred, log_scale = log_scale)
  }
  perf <- do.call(rbind, lapply(results, as.data.frame))
  perf
}