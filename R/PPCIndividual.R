#' Individual Posterior Predictive P-Values for the Continuous Model
#'
#' Computes, for each subject with a nonzero response, a posterior
#' predictive p-value checking whether its observed value lies within the
#' range predicted by the continuous (magnitude) part of a model fitted by
#' \code{\link{MainBVSSemi}}, and plots them.
#'
#' For every retained MCMC draw \eqn{s = 1, \dots, S}, with fitted values
#' \eqn{\mu_i^{(s)} = x_i^\top\beta^{(s)}} (stored in
#' \code{fit$fitted.cont.draws}) and residual variance
#' \eqn{\sigma^{2(s)}} (\code{fit$sigma2.draws}), a replicate
#' \eqn{y_i^{rep,(s)} \sim N(\mu_i^{(s)}, \sigma^{2(s)})} is drawn for each
#' nonzero subject \eqn{i}. The individual p-value is
#' \deqn{p_i = \frac{1}{S}\sum_{s=1}^S 1\{y_i^{rep,(s)} \ge y_i\},}
#' where \eqn{y_i} is the observed nonzero response on the scale of the fit
#' (\code{fit$y.cont}, i.e. \code{log(Y[Y != 0])} when the model was fit
#' with \code{log_scale = TRUE}).
#'
#' If the model fits well, the \eqn{p_i} should be spread roughly uniformly
#' over \eqn{(0, 1)} (in practice somewhat concentrated towards the middle,
#' since each \eqn{y_i} also informed the posterior). Values close to 0
#' flag subjects whose observed value is larger than the model predicts
#' (underpredicted); values close to 1 flag subjects whose observed value
#' is smaller than predicted (overpredicted).
#'
#' @param fit A list returned by \code{\link{MainBVSSemi}}, containing
#'   \code{y.cont}, \code{index.cont}, \code{fitted.cont.draws} and
#'   \code{sigma2.draws}.
#' @param lower,upper Thresholds below/above which a p-value is flagged as
#'   extreme, drawn as horizontal reference lines in the plot. Default to
#'   \code{0.05} and \code{0.95}.
#' @param plot Logical. If \code{TRUE} (the default), plots \eqn{p_i}
#'   against each subject's position in \code{Y}, with horizontal lines at
#'   \code{lower} and \code{upper}; flagged subjects are drawn in red and
#'   labelled with their position in \code{Y}.
#' @param seed Integer seed for drawing the replicates. The caller's
#'   random-number state is restored on exit.
#'
#' @return Invisibly, a list with:
#'   \describe{
#'     \item{\code{p.value}}{Length-\eqn{n_2} vector of the individual
#'       p-values \eqn{p_i}, named by each subject's position in \code{Y}.}
#'     \item{\code{index}}{Positions in \code{Y} of the nonzero subjects
#'       (\code{fit$index.cont}), in the order of \code{p.value}.}
#'     \item{\code{underpredicted}, \code{overpredicted}}{Positions in
#'       \code{Y} of the subjects with \eqn{p_i < } \code{lower} and
#'       \eqn{p_i > } \code{upper}, respectively.}
#'     \item{\code{prop.extreme}}{Proportion of nonzero subjects flagged
#'       (outside \code{[lower, upper]}).}
#'     \item{\code{lower}, \code{upper}}{The thresholds used.}
#'   }
#'
#' @seealso \code{\link{PPCHomoscedasticity}}, and the \code{PPC.bin} and
#'   \code{PPC.cont} components of \code{\link{MainBVSSemi}}'s output.
#'
#' @importFrom graphics plot points abline text
#' @export
#'
#' @examples
#' \dontrun{
#' Dat <- GenDataSemiContinous(n = 500, p = 500, sd = 1, impf = 20,
#'                              beta = 0.3, percentOverlap = "Full", seed = 1,
#'                              log_scale = TRUE)
#' fit <- MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X,
#'                     mcmcsample = 10000, burnin = 5000)
#' ind <- PPCIndividual(fit)
#' ind$prop.extreme
#' ind$underpredicted
#' }
PPCIndividual <- function(fit, lower = 0.05, upper = 0.95, plot = TRUE, seed = 1) {
  if (is.null(fit$y.cont) || is.null(fit$index.cont) ||
      is.null(fit$fitted.cont.draws) || is.null(fit$sigma2.draws)) {
    stop("fit must contain y.cont, index.cont, fitted.cont.draws and sigma2.draws; ",
         "refit with the current MainBVSSemi().")
  }
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != 1 ||
      length(upper) != 1 || !(0 <= lower && lower < upper && upper <= 1)) {
    stop("lower and upper must be single numbers with 0 <= lower < upper <= 1")
  }
  y <- fit$y.cont
  mu <- fit$fitted.cont.draws          # S x n2
  sigma <- sqrt(fit$sigma2.draws)      # length S
  S <- nrow(mu)
  n2 <- length(y)
  idx <- fit$index.cont

  ## Draw the replicates without disturbing the caller's RNG state
  if (exists(".Random.seed", envir = globalenv())) {
    oldSeed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", oldSeed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)

  ## Count, draw by draw, how often y_i^rep,(s) >= y_i (no S x n2 replicate
  ## matrix is kept)
  count <- rep(0, n2)
  for (s in seq_len(S)) {
    yrep <- mu[s, ] + rnorm(n2, sd = sigma[s])
    count <- count + (yrep >= y)
  }
  pval <- stats::setNames(count / S, idx)

  out <- list(p.value = pval, index = idx,
              underpredicted = idx[pval < lower],
              overpredicted = idx[pval > upper],
              prop.extreme = mean(pval < lower | pval > upper),
              lower = lower, upper = upper)

  if (plot) {
    extreme <- pval < lower | pval > upper
    plot(idx, pval, ylim = c(0, 1), pch = 19, cex = 0.6,
         col = ifelse(extreme, "red", "grey40"),
         xlab = "Subject (position in Y)",
         ylab = expression(p[i] == P(y[i]^rep >= y[i] ~ "|" ~ y)),
         main = sprintf("Individual posterior predictive p-values (%.1f%% outside [%g, %g])",
                        100 * out$prop.extreme, lower, upper),
         cex.main = 0.95)
    abline(h = c(lower, upper), lty = 2, col = "red")
    if (any(extreme)) {
      text(idx[extreme], pval[extreme], labels = idx[extreme], pos = 4,
           cex = 0.6, col = "red", xpd = TRUE)
    }
  }
  invisible(out)
}
