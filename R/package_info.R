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
