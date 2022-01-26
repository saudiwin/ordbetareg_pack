#' Normalize Outcome/Response to \\[0,1\\] Interval
#'
#' This function takes a continuous (double) column of data and converts it
#' to have 0 as the lower bound and 1 as the upper bound.
#'
#' Beta regression can only be done with a response that is continuous with a
#' lower bound of 0 and an upper bound of 1. However, it is
#' straightforward to transform any lower and upper-bounded continuous
#' variable to the \\[0,1\\] interval. This function does the transformation
#' and saves the original bounds as attributes so that the bounds can be
#' reverse-transformed.
#'
#' @param outcome Any non-character vector. Factors will be converted
#' to numeric via coercion.
#' @return A numeric vector with an upper bound of 1 and a lower bound of
#' 0. The original bounds are saved in the attributes "lower_bound" and
#' "upper_bound".
#' @examples
#' # set up arbitrary upper and lower-bounded vector
#' outcome <- runif(1000, min=-33, max=445)
#'
#' # normalize to \\[0,1\\]
#'
#' trans_outcome <- normalize(outcome=outcome)
#' summary(trans_outcome)
#'
#' # only works with numeric vectors and factors
#' \dontrun{
#'   normalize(outcome=c('a','b'))
#' }
#' @export
normalize <- function(outcome) {

  if(is.character(outcome)) {

    stop("Please do not pass a character vector to the function.\nThat really doesn't make any sense.")

  }

  if(is.factor(outcome)) {

    outcome <- as.numeric(outcome)

  }

  if(is.na(min(outcome, na.rm=T)) || is.infinite(min(outcome, na.rm=T))) {

    stop("The vector does not have enough non-missing data.")

  }

  trans_out <- (outcome - min(outcome, na.rm=T)) / (max(outcome, na.rm=T) - min(outcome, na.rm=T))

  attr(trans_out, "upper_bound") <- max(outcome, na.rm=T)
  attr(trans_out, "lower_bound") <- min(outcome, na.rm=T)

  return(trans_out)


}
