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
