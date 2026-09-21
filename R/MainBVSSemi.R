#' Bayesian Variable Selection for a Semicontinuous Response
#'
#' An MCMC algorithm to perform Bayesian variable selection for a
#' semicontinuous response via a two-part model: a continuous model (a
#' linear model for the nonzero response values) and a binary model (a
#' probit model for whether the response is zero or nonzero). Three
#' selection strategies are available: \code{"BVSSemiMRF"}, which encourages
#' common selection of important features between the two models via a
#' Markov random field prior linking the two models' inclusion indicators;
#' \code{"BVSSemiComb"}, which forces the same set of selected features to
#' be used in both the continuous and binary models; and
#' \code{"BVSSemiIndep"}, which assumes the sets of selected features for the
#' two models need not coincide and fits them independently. The algorithm
#' returns, among other things, the marginal posterior probability of
#' inclusion of each feature in each model. Both \code{X} and \code{Xcov}
#' are standardized column-wise (centered and scaled to unit standard
#' deviation) before fitting; the per-column means and standard deviations
#' used are returned so that \code{\link{PosteriorPredict}} can apply the
#' same standardization to new data.
#'
#' @param Method One of \code{"BVSSemiMRF"}, \code{"BVSSemiComb"} or
#'   \code{"BVSSemiIndep"}. Defaults to \code{"BVSSemiMRF"}.
#' @param Y A numeric vector giving the semicontinuous response (a mix of
#'   exact zeros and continuous positive values).
#' @param X A numeric matrix of candidate features of dimension
#'   \eqn{n \times p}, subject to variable selection. Standardized
#'   internally before fitting (see Description).
#' @param Xcov An optional numeric matrix or vector of covariates that are
#'   not subject to variable selection (e.g. clinical/demographic variables
#'   such as sex, age, race, etc.), always included in both models. Defaults
#'   to \code{NULL}. Standardized internally before fitting (see
#'   Description).
#' @param seed Integer seed used to initialize the random number generator
#'   for the MCMC algorithm.
#' @param atheta Shape (hyper)parameter of the gamma prior distribution of
#'   \code{theta}, used only by the \code{"BVSSemiMRF"} method. The
#'   parameter \code{theta} measures the strength of borrowing between the
#'   two models and encourages common feature selection between them.
#' @param btheta Rate (hyper)parameter of the gamma prior distribution of
#'   \code{theta}, used only by the \code{"BVSSemiMRF"} method.
#' @param tau2cont Variance (hyper)parameter of the normal prior
#'   distribution of the regression effects in the continuous model.
#' @param tau2bin Variance (hyper)parameter of the normal prior distribution
#'   of the regression effects in the binary model.
#' @param nu1cont Log-odds of the prior probability of feature inclusion in
#'   the continuous model. For \code{Method = "BVSSemiComb"}, the continuous
#'   and binary models share a single inclusion indicator, and
#'   \code{nu1cont} alone controls its prior log-odds (\code{nu2bin} is
#'   ignored in that case).
#' @param nu2bin Log-odds of the prior probability of feature inclusion in
#'   the binary model. Ignored when \code{Method = "BVSSemiComb"} (see
#'   \code{nu1cont}).
#' @param varpropTheta Variance of the Metropolis-Hastings proposal
#'   distribution for \code{theta}. It should be tuned to give an
#'   acceptance rate of roughly 20-60 percent. Defaults to \code{.25}.
#' @param Bigtau2 Variance (hyper)parameter of the normal prior distribution
#'   of the regression effects of features that are not subject to variable
#'   selection (e.g. the intercept and covariates in \code{Xcov}).
#' @param asigma Shape parameter of the inverse-gamma prior distribution of
#'   the continuous-model residual variance, \code{sigma2}.
#' @param bsigma Scale parameter of the inverse-gamma prior distribution of
#'   \code{sigma2}.
#' @param mcmcsample Total number of MCMC iterations to run. Must be larger
#'   than \code{burnin}.
#' @param burnin Number of initial MCMC iterations to discard as burn-in.
#' @param thin Thinning interval applied to the post-burn-in draws: every
#'   \code{thin}-th retained iteration is stored in the returned coefficient
#'   and variance draws. Defaults to \code{5}.
#' @param log_scale Logical. If \code{TRUE} (the default), the continuous
#'   model is fit to \code{log(Y)} restricted to the nonzero (positive)
#'   values of \code{Y}, instead of to \code{Y} itself. Requires every
#'   nonzero value of \code{Y} to be strictly positive; an error is raised
#'   otherwise. Pass the same value to \code{log_scale} in
#'   \code{\link{PosteriorPredict}} when predicting from the fitted model.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{prob.Z.Cont}}{Length-\eqn{p} vector of marginal posterior
#'       inclusion probabilities of each feature in the continuous model.}
#'     \item{\code{prob.Z.Bin}}{Length-\eqn{p} vector of marginal posterior
#'       inclusion probabilities of each feature in the binary model.}
#'     \item{\code{Gamma.cont.draws}, \code{Gamma.bin.draws}}{Matrices (one
#'       row per retained draw, one column per feature) of the continuous-
#'       and binary-model inclusion indicators for that draw. Identical to
#'       each other for \code{"BVSSemiComb"} (a single shared indicator).}
#'     \item{\code{logPost.draws}}{Length-\eqn{M} vector (\eqn{M} = number of
#'       retained draws) of each draw's total log-posterior of its discrete
#'       model \code{(Gamma.cont, Gamma.bin)} (and \code{theta} for
#'       \code{"BVSSemiMRF"}), up to a model-independent constant.}
#'     \item{\code{Gamma.cont.models}, \code{Gamma.bin.models}}{Matrices (one
#'       row per \emph{unique} visited model, one column per feature) of the
#'       distinct \code{(Gamma.cont, Gamma.bin)} configurations visited
#'       during MCMC -- i.e. \code{Gamma.cont.draws}/\code{Gamma.bin.draws}
#'       deduplicated.}
#'     \item{\code{logPost.models}}{Length-\eqn{K} vector (\eqn{K} = number
#'       of unique visited models) of each model's log-posterior, recomputed
#'       once as post-processing after the MCMC loop from the fixed training
#'       response \code{y2} (continuous part) and \code{U.mean} (binary
#'       part, below) rather than averaged from any single noisy iteration's
#'       value, plus the properly normalized Gamma-inclusion log-prior
#'       evaluated at the posterior mean \code{theta} (\code{"BVSSemiMRF"}
#'       only). Used by \code{\link{PosteriorPredict}}'s \code{nmodels}
#'       argument (Bayesian model averaging over the
#'       highest-posterior-probability visited models).}
#'     \item{\code{beta.Cont.mode}, \code{beta.Bin.mode}}{\eqn{K} x
#'       (1 + pc + p) matrices (one row per unique visited model) of the
#'       continuous- and binary-model posterior-mode coefficients for that
#'       model, computed once as post-processing after the MCMC loop (not
#'       retained per draw): the continuous part from the fixed training
#'       response, the binary part using \code{U.mean} (below) in place of
#'       any single iteration's latent draw. Used by
#'       \code{\link{PosteriorPredict}}'s \code{nmodels} argument.}
#'     \item{\code{sigma2.mean}}{A single number: the posterior mean of
#'       \code{sigma2.draws}, used as a common residual variance for every
#'       model by \code{\link{PosteriorPredict}}'s \code{nmodels} argument
#'       (unlike \code{beta.Cont.mode}/\code{beta.Bin.mode}, this does not
#'       vary by model).}
#'     \item{\code{U.mean}}{Length-\eqn{n} vector: the posterior mean of the
#'       binary submodel's latent augmentation variable \code{U} across all
#'       retained draws, used (in place of any single draw's noisy \code{U})
#'       to compute \code{beta.Bin.mode} and \code{logPost.models} for every
#'       visited model.}
#'     \item{\code{sigma2.draws}}{Thinned post-burn-in posterior draws of the
#'       continuous-model residual variance, \code{sigma2}.}
#'     \item{\code{theta.draws}}{Thinned post-burn-in posterior draws of
#'       \code{theta} from the \code{"BVSSemiMRF"} method (otherwise a
#'       vector of zeros, since \code{theta} is fixed at 0).}
#'     \item{\code{AcceptanceRateTheta}}{MCMC acceptance rate of \code{theta}
#'       for the \code{"BVSSemiMRF"} method, or \code{NULL} for the other
#'       two methods (which do not sample \code{theta}).}
#'     \item{\code{pc}}{Number of forced-in covariates in \code{Xcov}.}
#'     \item{\code{p}}{Number of candidate features in \code{X}.}
#'     \item{\code{log_scale}}{The \code{log_scale} argument used to fit the
#'       model, passed through for use by \code{\link{PosteriorPredict}}.}
#'     \item{\code{X.center}}{Length-\eqn{p} vector of the column means of
#'       \code{X} used to standardize it before fitting.}
#'     \item{\code{X.scale}}{Length-\eqn{p} vector of the column standard
#'       deviations of \code{X} used to standardize it before fitting.}
#'     \item{\code{Xcov.center}}{Length-\eqn{pc} vector of the column means
#'       of \code{Xcov} used to standardize it before fitting, or
#'       \code{NULL} when \code{Xcov} was not supplied.}
#'     \item{\code{Xcov.scale}}{Length-\eqn{pc} vector of the column
#'       standard deviations of \code{Xcov} used to standardize it before
#'       fitting, or \code{NULL} when \code{Xcov} was not supplied.}
#'   }
#'
#' @references
#' Samuel Babatunde, Tolulope Sajobi and Thierry Chekouo (2026),
#' \emph{A Bayesian Variable Selection for Semicontinuous Response data:
#' Application to cardiovascular disease}, submitted.
#'
#' @seealso \code{\link{GenDataSemiContinous}}, \code{\link{PosteriorPredict}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(BVSSemi)
#' ## log_scale = TRUE exponentiates the nonzero values of Y so they are
#' ## strictly positive, matching MainBVSSemi's default log_scale = TRUE
#' Dat <- GenDataSemiContinous(n = 500, p = 500, sd = 1, impf = 20,
#'                              beta = .3, percentOverlap = "Full", seed = 1,
#'                              log_scale = TRUE)
#' str(Dat)
#'
#' result <- MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X, seed = 1,
#'                        atheta = 1, btheta = 1, tau2cont = 0.05, tau2bin = 1,
#'                        nu1cont = -4, nu2bin = -4, varpropTheta = .25,
#'                        Bigtau2 = 100, mcmcsample = 40000, burnin = 1000,
#'                        thin = 5)
#' str(result)
#'
#' library(pROC)
#' AUC1 <- as.numeric(auc(roc(as.factor(Dat$Z.cont), result$prob.Z.Cont)))
#' AUC2 <- as.numeric(auc(roc(as.factor(Dat$Z.bin), result$prob.Z.Bin)))
#' print(AUC1); print(AUC2)
#' }
MainBVSSemi <- function(
    Method = "BVSSemiMRF", Y, X, Xcov = NULL, seed = 1,
    atheta = 1, btheta = 1, tau2cont = 1, tau2bin = 0.5, nu1cont = -3, nu2bin = -3, varpropTheta = .25,
    Bigtau2 = 100, asigma = .1, bsigma = .1, mcmcsample = 10000, burnin = 5000, thin=5,
    log_scale = TRUE) {
  if (Method == "BVSSemiComb" && !missing(nu2bin)) {
    warning("nu2bin is ignored when Method = \"BVSSemiComb\": the continuous ",
            "and binary models share a single inclusion indicator, whose prior ",
            "log-odds is controlled by nu1cont alone.")
  }
  tau21 <- tau2cont; tau22 <- tau2bin
  nu1 <- nu1cont; nu2 <- nu2bin
  set.seed(seed)
  n <- length(Y)
  N2 <- which(Y != 0)
  if (log_scale) {
    if (any(Y[N2] <= 0)) {
      stop("log_scale = TRUE requires all nonzero values of Y to be positive")
    }
    y2 <- log(Y[N2])
  } else {
    y2 <- Y[N2]
  }
  if (is.null(Xcov)) {
    pc <- 0
    Xcov.center <- NULL
    Xcov.scale <- NULL
  } else {
    if (is.vector(Xcov)) {
      Xcov <- matrix(Xcov, n, 1)
    } else {
      Xcov <- as.matrix(Xcov)
    }
    pc <- ncol(Xcov)
    Xcov.center <- colMeans(Xcov)
    Xcov.scale <- apply(Xcov, 2, sd)
    Xcov.scale[Xcov.scale == 0] <- 1
    Xcov <- sweep(sweep(Xcov, 2, Xcov.center, "-"), 2, Xcov.scale, "/")
  }
  X <- as.matrix(X)
  p <- ncol(X)
  X.center <- colMeans(X)
  X.scale <- apply(X, 2, sd)
  X.scale[X.scale == 0] <- 1
  X <- sweep(sweep(X, 2, X.center, "-"), 2, X.scale, "/")

  aa <- asigma; ba <- bsigma
  # aa=ba=0.001

  #####
  alpha1 <- atheta;  beta1 <- btheta

  ### Initialization
  Gamma1 <- rbinom(p, size = 1, prob = 1 / (1 + exp(-nu1))) ### Regression model
  Gamma2 <- rbinom(p, size = 1, prob = 1 / (1 + exp(-nu2))) # Binary model
  sigma2 <- .1
  if (Method == "BVSSemiMRF") {
    theta <- .5
  } else {
    theta <- 0
  }
  U <- rep(0, n)
  U[Y != 0] <- -abs(rnorm(sum(Y != 0)))
  U[Y == 0] <- abs(rnorm(sum(Y == 0)))

  #### Output values
  NN <- mcmcsample
  acceptTheta <- 0
  Gam1Mean <- Gam2Mean <- rep(0, p)
  Ubar <- rep(0, n) ## posterior mean of the latent U, used in place of any
                    ## single iteration's noisy U to compute betaBin's
                    ## posterior mode per model after the loop (see below)

   n_keep <- NN - burnin
  n_thin <- floor(n_keep / thin)
  sigma2Draws   <- rep(0, n_thin)
  thetaDraws    <- rep(0, n_thin)
  GammaContDraws <- matrix(0, n_thin, p)
  GammaBinDraws  <- matrix(0, n_thin, p)
  logPostDraws  <- rep(0, n_thin)
  keep_i <- 0
  NR <- length(N2)

  for (s in 1:NN) {
    # set.seed(s)
    if ((Method == "BVSSemiIndep") || (Method == "BVSSemiMRF")) {
      ### Sample Gamma1
       Xcov2 <- if (is.null(Xcov)) NULL else Xcov[N2, , drop = FALSE]
      Gamma1F <- SampleGammaLinear(
        GammaM1 = Gamma1, y = y2, X = X[N2, ,drop = FALSE], Xcov = Xcov2, pc = pc, asigma = asigma,
        bsigma = bsigma, tau2 = tau21,Bigtau2 = Bigtau2,
        theta = theta, GammaM2 = Gamma2, nu = nu1
      )
      Gamma1 <- Gamma1F$GammaM1
      sigma2 <- Sigma2(NR, Gamma1F$uSu, aa, ba)   ## sigma2 | Gamma1 (conjugate IG)
      ### Sample Gamma2 binary model
      Gamma2F <- SampleGamma(
        GammaM1 = Gamma2, y = U, X = X, Xcov = Xcov, pc = pc, sigma2 = 1, tau2 = tau22,
        Bigtau2 = Bigtau2, theta = theta, GammaM2 = Gamma1, nu = nu2
      )
      Gamma2 <- Gamma2F$GammaM1
      betaBin <- Gamma2F$beta ## still needed to update U below; not stored anywhere
      contLogl <- Gamma1F$logl
      binLogl <- Gamma2F$logl
    } else if (Method == "BVSSemiComb") {
      GammaF <- SampleGammaCombProb(
        N2 = N2, Gamma = Gamma1, U = U, y2 = y2, X = X, Xcov = Xcov, tau2 = tau21,
        Bigtau2 = Bigtau2, nu = nu1, sigma2 = sigma2, asigma = asigma, bsigma = bsigma
      )
      betaBin <- GammaF$betaBin ## still needed to update U below; not stored anywhere
      Gamma1 <- GammaF$Gamma
      sigma2 <- Sigma2(NR, GammaF$uSu, aa, ba)
      Gamma2 <- Gamma1
      contLogl <- GammaF$loglCont
      binLogl <- GammaF$loglBin
    }

     if (Method == "BVSSemiMRF") {
      ### Sample theta
      thetaF <- SampleTheta(theta, nu1, nu2, Gamma1, Gamma2, alpha1, beta1, varpropo = varpropTheta)
      theta <- thetaF$theta
      # thetaMean=thetaMean+theta/NN
      acceptTheta <- acceptTheta + thetaF$accept
    } else {
      theta <- 0
    }
    ## Per-iteration total log-posterior of the discrete model (Gamma1,
    ## Gamma2 [, theta]), up to a Gamma/theta-independent constant -- used
    ## only for cross-iteration model-posterior bookkeeping (e.g. Bayesian
    ## model averaging in PosteriorPredict), never for the MCMC transitions
    ## themselves (those are handled by the Sample* functions above).
    logPriorVal <- if (Method == "BVSSemiComb") logPriorShared(Gamma1, nu1) else logPriorMRF(Gamma1, Gamma2, theta, nu1, nu2)
    logPostTotal <- contLogl + binLogl + logPriorVal
    if (s > burnin) {
      Gam1Mean <- Gam1Mean + Gamma1 / n_keep
      Gam2Mean <- Gam2Mean + Gamma2 / n_keep
      Ubar <- Ubar + U / n_keep
      if ((s - burnin) %% thin == 0 && keep_i < n_thin) {
        keep_i <- keep_i + 1
        sigma2Draws[keep_i] <- sigma2
        thetaDraws[keep_i] <- theta
        GammaContDraws[keep_i, ] <- Gamma1
        GammaBinDraws[keep_i, ] <- Gamma2
        logPostDraws[keep_i] <- logPostTotal
      }
    }

    U <- YLatent2(Y, X, Xcov, betaBin)

    if (s %% (NN / 5) == 1) {
      print(paste("Number of mcmc samples is =", s))
      #  print(beta[1:30])
    }
   
  }
  
  GammaContDraws <- GammaContDraws[seq_len(keep_i), , drop = FALSE]
  GammaBinDraws <- GammaBinDraws[seq_len(keep_i), , drop = FALSE]
  logPostDraws <- logPostDraws[seq_len(keep_i)]
  thetaBar <- mean(thetaDraws[seq_len(keep_i)]) ## representative theta for the
                                                 ## log-prior of each model below
  sigma2Bar <- mean(sigma2Draws[seq_len(keep_i)]) ## single shared posterior
                                                   ## mean sigma2 across all models

  ## ---- Post-processing (after the MCMC loop): deterministic posterior
  ## mode beta, and log-posterior, for each UNIQUE (Gamma.cont, Gamma.bin)
  ## model visited, computed once per unique model rather than once per draw
  ## (far fewer computations, since a model is often revisited many times).
  ## The continuous part uses the fixed training response y2 (as during
  ## sampling); the binary part uses Ubar -- the posterior mean of the
  ## latent U across the whole chain -- in place of any single iteration's
  ## noisy U, since U is resampled every iteration even for a fixed
  ## Gamma.bin. This is the sole source of coefficients for prediction (see
  ## PosteriorPredict's nmodels), by design: no per-draw beta is retained.
  keys <- apply(cbind(GammaContDraws, GammaBinDraws), 1, paste, collapse = "")
  modelIdx <- split(seq_along(keys), keys)
  nModels <- length(modelIdx)
  firstOfEach <- vapply(modelIdx, `[`, integer(1), 1)
  GammaContModels <- GammaContDraws[firstOfEach, , drop = FALSE]
  GammaBinModels <- GammaBinDraws[firstOfEach, , drop = FALSE]

  Xcov2 <- if (is.null(Xcov)) NULL else Xcov[N2, , drop = FALSE]
  tau2Bin <- if (Method == "BVSSemiComb") tau21 else tau22
  betaContModel <- matrix(0, nModels, p + pc + 1)
  betaBinModel  <- matrix(0, nModels, p + pc + 1)
  logPostModels <- rep(0, nModels)
  for (m in seq_len(nModels)) {
    g1 <- GammaContModels[m, ]
    g2 <- GammaBinModels[m, ]
    contF <- loglikLinear(y2, X[N2, , drop = FALSE], Xcov2, g1, asigma, bsigma, tau21, Bigtau2)
    betaContModel[m, ] <- ExpandBeta(contF$betaMean, g1, pc)
    binF <- loglik(Ubar, X, Xcov, g2, 1, tau2Bin, Bigtau2)
    betaBinModel[m, ] <- ExpandBeta(binF$betaMean, g2, pc)
    logPriorVal <- if (Method == "BVSSemiComb") logPriorShared(g1, nu1) else logPriorMRF(g1, g2, thetaBar, nu1, nu2)
    logPostModels[m] <- contF$logl + binF$logl + logPriorVal
  }

  AcceptanceRateTheta <- if (Method == "BVSSemiMRF") acceptTheta / NN else NULL
  return(list(prob.Z.Cont = Gam1Mean, prob.Z.Bin = Gam2Mean,
  Gamma.cont.draws = GammaContDraws,
  Gamma.bin.draws = GammaBinDraws,
  logPost.draws = logPostDraws,
  Gamma.cont.models = GammaContModels,
  Gamma.bin.models = GammaBinModels,
  logPost.models = logPostModels,
  beta.Cont.mode = betaContModel,
  beta.Bin.mode = betaBinModel,
  sigma2.mean = sigma2Bar,
  U.mean = Ubar,
  sigma2.draws = sigma2Draws[seq_len(keep_i)],
  theta.draws = thetaDraws[seq_len(keep_i)],
  AcceptanceRateTheta = AcceptanceRateTheta, pc = pc, p = p,
  log_scale = log_scale, X.center = X.center, X.scale = X.scale,
  Xcov.center = Xcov.center, Xcov.scale = Xcov.scale))
}
