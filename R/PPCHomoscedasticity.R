#' Posterior Predictive Check of Homoscedasticity for the Continuous Model
#'
#' Checks the constant-variance assumption of the continuous (magnitude)
#' part of a model fitted by \code{\link{MainBVSSemi}}, using Bayesian
#' residuals and a posterior predictive p-value (PPP-value) based on a
#' discrepancy statistic.
#'
#' For every retained MCMC draw \eqn{s = 1, \dots, S}, with fitted values
#' \eqn{\mu^{(s)} = X\beta^{(s)}} for the nonzero subjects (stored in
#' \code{fit$fitted.cont.draws}) and residual variance \eqn{\sigma^{2(s)}}
#' (\code{fit$sigma2.draws}):
#' \enumerate{
#'   \item the residuals are \eqn{\epsilon^{(s)} = y - \mu^{(s)}}, where
#'     \eqn{y} is the nonzero response on the scale of the fit
#'     (\code{fit$y.cont}, i.e. \code{log(Y[Y != 0])} when the model was fit
#'     with \code{log_scale = TRUE});
#'   \item the discrepancy is the rank correlation between the absolute
#'     residuals and the fitted values,
#'     \eqn{T(y, \beta^{(s)}) = \mathrm{corr}(|\epsilon^{(s)}|, \mu^{(s)})};
#'   \item a replicated dataset \eqn{y^{rep,(s)} \sim N(\mu^{(s)},
#'     \sigma^{2(s)} I)} is drawn, and the same statistic computed,
#'     \eqn{T(y^{rep,(s)}, \beta^{(s)}) = \mathrm{corr}(|y^{rep,(s)} -
#'     \mu^{(s)}|, \mu^{(s)})};
#'   \item the PPP-value is
#'     \eqn{\frac{1}{S}\sum_s 1\{T(y^{rep,(s)}, \beta^{(s)}) \ge
#'     T(y, \beta^{(s)})\}}.
#' }
#' A PPP-value near 0.5 indicates that the observed residual spread is
#' consistent with constant variance. Values near 0 indicate that the
#' absolute residuals grow with the fitted values more than the model
#' allows (a "funnel" opening to the right); values near 1 indicate that
#' they shrink with the fitted values.
#'
#' Draws whose fitted values are all equal (e.g. a draw with no selected
#' features and no forced-in covariates, so that only the intercept is
#' fitted) give an undefined correlation and are dropped; their number is
#' returned in \code{n.draws.dropped}.
#'
#' @param fit A list returned by \code{\link{MainBVSSemi}}, containing
#'   \code{y.cont}, \code{fitted.cont.draws} and \code{sigma2.draws}.
#' @param method Correlation used for the discrepancy statistic, passed to
#'   \code{\link[stats]{cor}}: \code{"spearman"} (the default, robust to
#'   outliers and nonlinearity) or \code{"pearson"}.
#' @param level Probability level of the pointwise credible intervals of
#'   the residuals shown in the residual plot. Defaults to \code{0.95}.
#' @param plot Logical. If \code{TRUE} (the default), draws two panels: the
#'   posterior mean residuals \eqn{\bar\epsilon} against the posterior mean
#'   fitted values \eqn{\bar\mu}, with each point's \code{level} credible
#'   interval of its residual as a vertical bar; and the replicated
#'   discrepancies \eqn{T(y^{rep,(s)}, \beta^{(s)})} against the observed
#'   ones \eqn{T(y, \beta^{(s)})}, with the line \eqn{y = x} (the PPP-value
#'   is the share of points on or above it).
#' @param seed Integer seed for drawing the replicated datasets. The
#'   caller's random-number state is restored on exit.
#'
#' @return Invisibly, a list with:
#'   \describe{
#'     \item{\code{T.obs}, \code{T.rep}}{Vectors (one entry per retained
#'       draw used) of the observed and replicated discrepancies.}
#'     \item{\code{p.value}}{The PPP-value, \code{mean(T.rep >= T.obs)}.}
#'     \item{\code{p.value.two.sided}}{\eqn{2\min\{P(T^{rep} \ge T),
#'       P(T^{rep} \le T)\}}, capped at 1.}
#'     \item{\code{fitted.mean}, \code{resid.mean}}{Length-\eqn{n_2}
#'       vectors of the posterior mean fitted values \eqn{\bar\mu} and
#'       residuals \eqn{\bar\epsilon} of the nonzero subjects.}
#'     \item{\code{resid.lower}, \code{resid.upper}}{Length-\eqn{n_2}
#'       vectors of the bounds of each residual's \code{level} credible
#'       interval.}
#'     \item{\code{method}, \code{level}}{The arguments used.}
#'     \item{\code{n.draws.dropped}}{Number of draws dropped because their
#'       correlation was undefined (see Details).}
#'   }
#'
#' @seealso \code{\link{MainBVSSemi}} (whose \code{PPC.bin} and
#'   \code{PPC.cont} components give posterior predictive checks of the
#'   proportion of nonzero values and of the mean and median of the
#'   nonzero values).
#'
#' @importFrom graphics par plot segments points abline mtext
#' @importFrom stats median
#' @export
#'
#' @examples
#' \dontrun{
#' Dat <- GenDataSemiContinous(n = 500, p = 500, sd = 1, impf = 20,
#'                              beta = 0.3, percentOverlap = "Full", seed = 1,
#'                              log_scale = TRUE)
#' fit <- MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X,
#'                     mcmcsample = 10000, burnin = 5000)
#' hc <- PPCHomoscedasticity(fit)
#' hc$p.value
#' }
PPCHomoscedasticity <- function(fit, method = c("spearman", "pearson"),
                                level = 0.95, plot = TRUE, seed = 1) {
  method <- match.arg(method)
  if (is.null(fit$y.cont) || is.null(fit$fitted.cont.draws) || is.null(fit$sigma2.draws)) {
    stop("fit must contain y.cont, fitted.cont.draws and sigma2.draws; ",
         "refit with the current MainBVSSemi().")
  }
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number strictly between 0 and 1")
  }
  y <- fit$y.cont
  mu <- fit$fitted.cont.draws          # S x n2
  sigma2 <- fit$sigma2.draws           # length S
  S <- nrow(mu)
  n2 <- length(y)

  ## Draw the replicates without disturbing the caller's RNG state
  if (exists(".Random.seed", envir = globalenv())) {
    oldSeed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", oldSeed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)

  resid <- sweep(-mu, 2, y, "+")       # S x n2: y - mu^(s), row-wise
  Tobs <- Trep <- rep(NA_real_, S)
  for (s in seq_len(S)) {
    m <- mu[s, ]
    yrep <- m + rnorm(n2, sd = sqrt(sigma2[s]))
    if (stats::sd(m) > 0) {
      Tobs[s] <- stats::cor(abs(resid[s, ]), m, method = method)
      Trep[s] <- stats::cor(abs(yrep - m), m, method = method)
    }
  }
  ok <- !is.na(Tobs) & !is.na(Trep)
  if (!any(ok)) {
    stop("the discrepancy is undefined for every draw: the fitted values ",
         "are constant across subjects in all retained draws (no selected ",
         "features and no forced-in covariates)")
  }
  Tobs <- Tobs[ok]
  Trep <- Trep[ok]
  pval <- mean(Trep >= Tobs)

  alpha <- (1 - level) / 2
  out <- list(
    T.obs = Tobs, T.rep = Trep, p.value = pval,
    p.value.two.sided = min(1, 2 * min(pval, mean(Trep <= Tobs))),
    fitted.mean = colMeans(mu), resid.mean = colMeans(resid),
    resid.lower = apply(resid, 2, stats::quantile, probs = alpha, names = FALSE),
    resid.upper = apply(resid, 2, stats::quantile, probs = 1 - alpha, names = FALSE),
    method = method, level = level, n.draws.dropped = sum(!ok)
  )

  if (plot) {
    oldPar <- par(mfrow = c(1, 2))
    on.exit(par(oldPar), add = TRUE)
    scaleLab <- if (isTRUE(fit$log_scale)) " (log scale)" else ""
    plot(out$fitted.mean, out$resid.mean, type = "n",
         ylim = range(out$resid.lower, out$resid.upper),
         xlab = paste0("Posterior mean fitted value", scaleLab),
         ylab = "Posterior mean residual",
         main = "Residuals vs fitted")
    segments(out$fitted.mean, out$resid.lower, out$fitted.mean, out$resid.upper,
             col = "grey70")
    points(out$fitted.mean, out$resid.mean, pch = 19, cex = 0.6)
    abline(h = 0, lty = 2, col = "red")
    mtext(sprintf("bars: %g%% credible intervals", 100 * level), side = 3,
          line = 0.2, cex = 0.8)

    lim <- range(Tobs, Trep)
    plot(Tobs, Trep, xlim = lim, ylim = lim, pch = 19, cex = 0.5,
         col = "grey40",
         xlab = expression(T(y, beta^(s))),
         ylab = expression(T(y^rep, beta^(s))),
         main = "Discrepancy: corr(|residual|, fitted)")
    abline(0, 1, lty = 2, col = "red")
    mtext(sprintf("PPP-value = %.3f", pval), side = 3, line = 0.2, cex = 0.8)
  }
  invisible(out)
}
