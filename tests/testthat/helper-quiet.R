## MainBVSSemi() prints MCMC progress via print(); silence it in tests so
## test output only shows testthat's own reporting. Relies on R's lazy
## evaluation of `expr` (a promise), so it isn't evaluated until after sink()
## has redirected stdout.
quietly <- function(expr) {
  sink(tempfile())
  on.exit(sink())
  force(expr)
}
