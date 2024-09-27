# ordbetareg v0.7.3.1

- removed `transformr` from imports because of many dependencies
- fixed bug if the intercept was specified and centered due to `brms` code change
- fixed bug that prevented the `make_stancode` option from working. changed that option 
to `return_stancode`

# ordbetareg v0.7.2

- `glmmTMB` no longer has standard errors from the `marginaleffects` package. Updated vignette to reflect this.

# ordbetareg v0.7.1

- Removed dependency on `faux` as `faux` was taken off CRAN. This also means the 
`sim_ordbeta` function cannot generate correlated covariates.
- Fixed a bug with the `transform` option to the `slopes` function in `marginaleffects`
in the vignette.
- `marginaleffects` now supports standard errors for ordered beta regression models
 in glmmTMB, so vignette was updated to reflect this.

# ordbetareg v0.7.0

- Allow intercepts to receive separate priors from main effects (permitting zero-ing out intercepts).
- Changed log-likelihood calculation to match `brms` documentation and permit point-wise 
`loo` calculation.
- Added animated density plots for continuous responses.
- Improved plot formatting and added theme and label options.
- Updated vignette to include information about `glmmTMB` as an alternative for estimation
and new plot functions.
- Added distribution `rordbeta` and `dordbeta` functions.
- Removed legacy `0 + Intercept` function code. Now accepts any kind of formula.

# ordbetareg v0.5.0

- Fixed bug disabling extra priors for phi regression/modeling.
- Added `pp_check_ordbetareg` function for accurate posterior predictive plots for discrete and continuous outcomes.
- Updated vignette and documentation.

# ordbetareg v0.3.0

- Fixed bugs relating to processing of bounded outcomes. 
- Allow for user-specific bounds to be set for normalization of outcome to 0/1 (issue #2).
- Enable multiple imputation and multivariate responses with `brms`.
- Fix documentation issues (issues #4 and #5).

# ordbetareg v0.2.1

- Added RNG code for the `induced_dirichlet` prior on cutpoints to allow for 
  the `sample_prior` option to work with `brm`, allowing for Bayes factors
  and other post-estimation tools with the `bridgesampling` package (issue #1).
- Fixed error in terms of correctly accounting for reverse transform of the
  ordered cutpoint vector.
- Made the `ldpf` function more stable by using logs instead of direct
  multiplication.
- Fixed typos in the vignette.

# ordbetareg v0.2

- Initial release.
