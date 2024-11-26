
library(DeclareDesign)
library(ordbetareg)
library(fabricatr)
library(future)
library(marginaleffects)
library(broom)
library(glmmTMB)

set.seed(20241125)

options(parallelly.fork.enable = TRUE) # required for use in RStudio

# parallelism
plan(multicore,workers=4)

# helper function for glmmTMB ordbetareg fit

tidy_avg_slopes <- function(x) {
  tidy(avg_slopes(x))
}

# create two designs with identical treatment effects/expected values
# first uses rordbeta to generate proportion in [0,1]
# second uses rordbeta to generate proportion in [0,1], then
# dichotomizes to 0 or 1 by rounding at 0.5
# compare to see which has greater power given same number of obsevations
# & check for bias

# first, a simulated proportion using rordbeta (ordered beta distribution)
# see https://osf.io/preprints/socarxiv/2sx6y
# cutpoints = number of observations at the bound (i.e. 0 or 1)
# phi = dispersion, a value of 1 means relatively "flat"

design_rordbeta <-
  declare_model(
    N = 100,
    potential_outcomes(Y ~ rordbeta(N, mu = .5 + .05*Z,phi = 1,
                                    cutpoints=c(-3,3)
    ))
  ) +
  declare_inquiry(ATE = 0.05) +
  declare_assignment(Z = complete_ra(N, m = 50)) +
  declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
  declare_estimator(Y ~ Z, .method = glmmTMB::glmmTMB, inquiry = "ATE",
                    family=glmmTMB::ordbeta,
                    .summary= tidy_avg_slopes,
                    term="Z")

# equivalent to dichotomizing (if a proportion)
design_dichot <-
  declare_model(
    N = 100,
    potential_outcomes(Y ~ as.numeric(rordbeta(N, mu = .5 + .05*Z,phi = 1,
                                               cutpoints=c(-3,3))>0.5))
  ) +
  declare_inquiry(ATE = 0.05) +
  declare_assignment(Z = complete_ra(N, m = 50)) +
  declare_measurement(Y = reveal_outcomes(Y ~ Z)) +
  declare_estimator(Y ~ Z, .method = lm_robust, inquiry = "ATE")

# NB: DeclareDesign is using lm_robust as it's default estimation
# However, it should be unbiased for the mean/expected value for the
# binomial/dichotomized model

diagnosands <-
  declare_diagnosands(bias = mean(estimate - estimand),
                      power = mean(p.value <= 0.05))

# compare in terms of bias on the ATE & Power

diagnose_design(design_rordbeta, diagnosands = diagnosands)
diagnose_design(design_dichot, diagnosands = diagnosands)

# rordbeta has greater power than dichotomization (about 50% more power)
# dichotomized response also has bias - bias is positive, possibly because ATE is positive and
# greater than 0.5
