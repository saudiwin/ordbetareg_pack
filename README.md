README
================

This is the repository for the `ordbetareg` package, which can be used
to fit the ordered beta regression model as specified in Kubinec (2022)
(see paper <https://osf.io/preprints/socarxiv/2sx6y/>). This package is
[on
CRAN](https://CRAN.R-project.org/package=ordbetareg),
but to use the development version, use this install command (requires
the `remotes` package to be installed first):

    remotes::install_github("saudiwin/ordbetareg_pack",build_vignettes=TRUE,dependencies = TRUE)

To learn more about beta regression modeling and where this package fits
in, you can read [my blog post on the
topic](https://www.robertkubinec.com/post/limited_dvs/). To learn how to
use the package, see the introduction vignette by using the
`vignette("package_introduction", package="ordbetareg")` command at the
R console or viewing the [vignette on
CRAN](https://CRAN.R-project.org/package=ordbetareg).

If you use the package, please cite it as:

Kubinec, Robert. “Ordered Beta Regression: A Parsimonious, Well-Fitting
Model for Continuous Data with Lower and Upper Bounds.” *Political
Analysis*. 2022. <https://doi.org/10.1017/pan.2022.20>.
