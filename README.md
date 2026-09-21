# BVSSemi

`BVSSemi` is an R package implementing an MCMC algorithm for **Bayesian
variable selection for a semicontinuous response** — a response variable
that is a mixture of exact zeros and continuous positive values (e.g.
healthcare costs, biomarker concentrations, or other outcomes with an excess
of zeros).

The model has two linked parts:

- a **binary (occurrence) part** — a probit model for whether the response
  is zero or nonzero;
- a **continuous (magnitude) part** — a linear model for the value of the
  response given that it is nonzero (optionally fit on the log scale).

Each part has its own set of candidate predictors and its own variable
selection indicator, and the package provides three strategies for linking
the two parts together:

| Method           | Description |
|------------------|-------------|
| `"BVSSemiMRF"`   | (default) Links the two parts through a Markov random field prior with a learned interaction parameter `theta`: selecting a feature in one part makes it more likely to be selected in the other, without forcing them to match. |
| `"BVSSemiComb"`  | Forces the two parts to share exactly the same selected features. |
| `"BVSSemiIndep"` | Fits the two parts completely independently. |

The main output of the algorithm is the **marginal posterior probability of
inclusion** of each candidate feature in each part of the model, along with
posterior draws of the regression coefficients and residual variance that
can be used for out-of-sample prediction.

License: GPL (>= 2.0)

For more information please contact tchekouo@umn.edu

## Installation

Using the `devtools` package:

```r
install.packages("devtools")
library(devtools)
install_github("chekouo/BVSSemi")
```

## Quick start

```r
library(BVSSemi)

## Simulate a semicontinuous response: n = 500 subjects, p = 500 candidate
## features, of which 20 are truly important to both the continuous and
## binary parts of the model (percentOverlap = "Full"). log_scale = TRUE
## exponentiates the nonzero response values so they are strictly positive,
## as required to fit with MainBVSSemi's (default) log_scale = TRUE.
Dat <- GenDataSemiContinous(n = 500, p = 500, sd = 1, impf = 20, beta = 0.3,
                             percentOverlap = "Full", seed = 1,
                             log_scale = TRUE)

## Fit the MRF-linked two-part model
fit <- MainBVSSemi(Method = "BVSSemiMRF", Y = Dat$Y, X = Dat$X,
                    mcmcsample = 10000, burnin = 5000)

## Marginal posterior inclusion probability of each feature
head(fit$prob.Z.Cont)  # continuous (magnitude) part
head(fit$prob.Z.Bin)   # binary (occurrence) part

## Predict for new subjects (Bayesian model averaging over the top 10
## visited feature-selection models) and evaluate performance
pred <- PosteriorPredict(fit, Xnew = Dat$X, nmodels = 10)
EvaluatePrediction(Dat$Y, pred)

## K-fold cross-validated predictive performance
CVPredictBVSSemi(Y = Dat$Y, X = Dat$X, K = 5, nmodels = 10, Method = "BVSSemiMRF",
                  mcmcsample = 10000, burnin = 5000)
```

## Main functions

| Function | Purpose |
|----------|---------|
| `MainBVSSemi()` | Runs the MCMC algorithm and returns marginal posterior inclusion probabilities, posterior draws of the visited feature-selection models and (for `"BVSSemiMRF"`) the MRF interaction parameter `theta`, plus each visited model's posterior-mode coefficients and residual variance (for `PosteriorPredict()`'s Bayesian model averaging). |
| `GenDataSemiContinous()` | Simulates semicontinuous data under the model described in the reference below, with configurable overlap between the important features of the two parts. |
| `PosteriorPredict()` | Computes posterior predictive summaries (probability of a nonzero response, predicted magnitude, predicted mean) for new subjects from a fitted model. |
| `EvaluatePrediction()` | Scores predictions against a held-out response: AUC and Brier score for the binary part; RMSE, MAE, correlation, and interval coverage for the continuous part; combined RMSE/MAE overall. |
| `CVPredictBVSSemi()` | Runs K-fold cross-validation of the full fit-predict-evaluate pipeline and returns one row of performance metrics per fold. |

See the PDF manual (`BVSSemi-manual.pdf`) or each function's help page
(e.g. `?MainBVSSemi`) for full argument and return-value documentation.

## Reference

Samuel Babatunde, Tolulope Sajobi and Thierry Chekouo (2026), *A Bayesian
Variable Selection for Semicontinuous Response Data: Application to
Cardiovascular Disease*, submitted.
