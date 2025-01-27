# Ordered Beta Regression


This website contains information about the R package `ordbetareg` and
general information about other implementations of the ordered beta
regression model in other statistical packages. If you just want to get
started with the R package `ordbetareg`, see the
[vignette](vignettes/package_introduction.html).

## About the Model

Proportion data–percentiles, slider scales, visual-analog scales–have
long puzzled statisticians because they combine two things: a continuous
measure (the proportion) and a discrete measure (the top value of 1 or
100% and the bottom value of 0 or 0%). Conventional statistical models
like the oft-used OLS regression implicitly assume these boundaries do
not exist; this means an OLS regression can predict to absurd values
like 115% of patients being sick or -5% of legislators being elected.

To address this issue, I developed the ordered beta regression model
that combines two things: beta regression, which is defined for any
bounded continuous scale, and ordered logit, which works for discrete
categories. By doing this, you can fit an ordered beta regression for
any percentile/proportion for *both* the middle continuous part and the
bounds of the scale. This model can also predict within this scale as
well.

I describe the model and how to do this type of [modeling in this blog
post](https://www.robertkubinec.com/post/limited_dvs/). For more detail
on the inner workings of the model, I refer you to the paper:

Kubinec, Robert. 2023. “Ordered Beta Regression: A Parsimonious,
Well-Fitting Model for Continuous Data with Lower and Upper Bounds.”
Political Analysis 31(4): 519–36. doi:
<http://doi.org/10.1017/pan.2022.20>.

For an ungated version of the article, see the pre-publication draft on
OSF:

<https://osf.io/preprints/socarxiv/2sx6y>

## Support ordbetareg

If you’ve enjoyed using the package, consider buying some `ordbetareg`
swag! Click the image below to go to the online store (I receive a few
dollars from each order).

## How to Use the Model

At present, ordered beta regression is available in Stan
(<a href="#sec-stan" class="quarto-xref">Section 3.5</a>), R
(<a href="#sec-r" class="quarto-xref">Section 3.1</a>), Python
(<a href="#sec-python" class="quarto-xref">Section 3.2</a>), Stata
(<a href="#sec-stata" class="quarto-xref">Section 3.3</a>), and Julia
(<a href="#sec-julia" class="quarto-xref">Section 3.4</a>). R has the
strongest support for the model with multiple implementations, but it is
perfectly usable in other statistical frameworks. Below I briefly
describe the available packages.

### R

Ordered beta regression is currently available in R in two different
packages:
[`ordbetareg`](https://cran.r-project.org/web/packages/ordbetareg/index.html)
and
[`glmmTMB`](https://cran.r-project.org/web/packages/glmmTMB/index.html).
The primary difference between these two packages is estimation method
and the number of ordered beta-specific functions. `ordbetareg` is a
package that estimates a Bayesian implementation as described in the
paper above by using `brms`, a powerful R package for Bayesian
regression. `ordbetareg` is maintained by me and you can see the full
source code here on
[Github](https://github.com/saudiwin/ordbetareg_pack) (and report [an
issue with the
package](https://github.com/saudiwin/ordbetareg_pack/issues) if you have
one). The package includes auxiliary functions like power analysis and
plots that are specific to proportion responses—check out [the package
vignette](https://cran.r-project.org/web/packages/ordbetareg/index.html)
for more details.

`glmmTMB`, by contrast, is a general purpose regression package that
specializes in mixed multilevel models. It implements a broad array of
regression models using maximum likelihood. This means that a model
estimated in `glmmTMB` will almost certainly be faster than `ordbetareg`
(although if you set up `ordbetareg` to use
[`cmdstanr`](https://mc-stan.org/cmdstanr/) it can get pretty fast). On
the other hand, maximum likelihood estimation has some limitations.
Arguably the most important one is that it can be tricky to fit a model
that doesn’t have observations at the bounds, i.e. either 0 (0%) or 1
(100%). The Bayesian implementation in `ordbetareg` has no problem doing
this.

That being said, I think `glmmTMB` is a fine package and believe it
works fine for many scholars’ problems. `brms`, which `ordbetareg` is
based on, offers a lot more features, including native support for
multiple imputation, time series modeling and even latent
variable/factor analysis modeling, but it does require more setup and
the models are generally slower.

For either package, I highly recommend using
[`marginaleffects`](https://marginaleffects.com/), and in particular the
`avg_slopes` function, to convert the coefficient estimates in the
package to the bounded scale (i.e. between 0 and 1). This function
allows you to then interpret your model estimates as the effect of a
covariate in terms of percentage change/proportion change in the bounded
response/outcome. `marginaleffects` is available [on
CRAN](https://cran.r-project.org/web/packages/marginaleffects/index.html)
and works great with both `ordbetareg` and `glmmTMB`.

With both of these packages available, ordered beta regression has
**very strong** support in R. If you don’t have a preference for a
particular software package and want to use this model, I would
recommend R.

### Python

Ordered beta regression is implemented using the `scipy` package with
maximum likelihood estimation. At present, you’ll need to clone (i.e.,
download) this [Github
repository](https://github.com/saudiwin/ordbetareg_py) to install the
package as it is not yet available with `pip` or `conda forge`. The
package includes plotting functions specific to proportion outcomes and
reports coefficient values in the untransformed (logit) scale.
`marginaleffects` support for this package is coming soon.

### Stata

There is initial support for Stata as a set of .ado files that can be
downloaded and installed in the `ado` folder in your Stata machine. More
details are available on the [Github
repository](https://github.com/saudiwin/ordbetareg_stata). This package
supports using the `margins` command to convert coefficients to the
bounded outcome scale, although it does not support all Stata features
as of yet. There are plans to make this package available via SSC in the
near future.

### Julia

Ordered beta regression is available via the
[`SubjectiveScalesModels.jl`](https://github.com/DominiqueMakowski/SubjectiveScalesModels.jl/tree/825f3fc089e64361e72e1a336003ea78320dc496)
package, which offers support for a variety of regression models useful
in cognitive psychology and other fields with sliders/visual analog
scales. The package includes neat visualizations and is maintained by
[Dominique Makowksi](https://dominiquemakowski.github.io/).

### Stan

As is evident in the paper, the original model was written in Stan code.
While I no longer am maintaining the Stan code, [it is available in the
paper
repository](https://github.com/saudiwin/ordbetareg/blob/master/beta_logit.stan)
and works great. If you want to incorporate the likelihood or other Stan
code into your Stan model, feel free!
