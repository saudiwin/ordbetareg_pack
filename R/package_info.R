utils::globalVariables(c('filter','resp'))

#' ordbetareg: A Model for Continuous Data with Lower and Upper Bounds
#'
#' The `ordbetareg` package is essentially a wrapper around `brms` that
#' enables the ordered beta regression model to be fit. This model has
#' advantages over other alternatives for continous data with upper
#' and lower bounds, such as survey sliders, indexes,
#' dose-response relationships,
#' and visual analog scales (among others). The package allows for all of the
#' many `brms` regression modeling functions to be used with the ordered
#' beta regression distribution.
#'
#' To learn more about how the package works, see the vignette by using
#' the command `browseVignettes(package='ordbetareg')`.
#'
#' For more info about the distribution, see
#' this paper: https://osf.io/preprints/socarxiv/2sx6y/
#'
#' To cite the package, please cite the following paper:
#'
#' Kubinec, Robert. "Ordered Beta Regression: A Parsimonious, Well-Fitting Model for Continuous Data with Lower and Upper Bounds." **Political Analysis**. 2022.
#' @docType package
#' @name ordbetareg
NULL


#' Pew American Trends Panel Wave 28
#'
#' A dataset with the non-missing responses for the 28th wave of the
#' Pew American Trends Panel survey.
#'
#' @format A data frame with 140 variables and 2,538 observations.
#' @source https://www.pewresearch.org/social-trends/dataset/american-trends-panel-wave-28/]
"pew"

#' Fitted Ordered Beta Regression Model
#'
#' A fitted ordered beta regression model to the mean of the thermometer
#' column from the pew data.
#' @format an `ordbetareg` object
"ord_fit_mean"

#' Fitted Ordered Beta Regression Model (Phi Regression)
#'
#' A fitted ordered beta regression model to the dispersion parameter
#' of the thermometer
#' column from the pew data.
#' @format an `ordbetareg` object
"ord_fit_phi"

#' Simulated Ordered Beta Regression Values
#'
#' The simulated draws used in the vignette for calculating statistical
#' power.
#' @format A dataframe
"sim_data"
