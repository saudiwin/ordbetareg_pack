utils::globalVariables(c('filter','resp','deparse0','SW',
                         'is_equal','all_vars','is.brmsfit',
                         'ulapply','combine_family_info',
                         'group_to_factor','family_names',
                         'get_dpar','softmax','component',
                         '.draw','rownum','pred',
                         'y','median_prop','y_agg',"y_plot",
                         'yplot','ymin','ymax','prop_ci_nice',
                         'restructure','prepare_predictions',
                         'collapse_comma',"is.bprepl","posterior_epred_trunc",
                         "is.bprepnl",
                         '.dpar_family','dpar_class',
                         'dpar_id','get_nlpar','is.mixfamily',
                         'conv_cats_dpars','seq_dim','p'))

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
#' @keywords internal
"_PACKAGE"


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

#' Fitted Ordered Beta Regression Model (Multivariate regression)
#'
#' A fitted ordered beta regression model with two responses,
#' one an ordered beta regression and the other a Gaussian/Normal
#' outcome. Useful for examining mediation analysis.
#' @format an `ordbetareg` object
"fit_multivariate"

#' Simulated Ordered Beta Regression Values
#'
#' The simulated draws used in the vignette for calculating statistical
#' power.
#' @format A dataframe
"sim_data"
