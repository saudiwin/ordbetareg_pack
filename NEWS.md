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
