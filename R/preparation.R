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
normalize <- function(outcome,true_bounds=NULL) {

  if(is.character(outcome)) {

    stop("Please do not pass a character vector as a response/outcome.\nThat really doesn't make any sense.")

  }

  if(is.factor(outcome)) {

    print("Converting factor response variable to numeric.")

    outcome <- as.numeric(outcome)

  }

  if(is.na(min(outcome, na.rm=T)) || is.infinite(min(outcome, na.rm=T))) {

    stop("The vector does not have enough non-missing data.")

  }

  if(!is.null(true_bounds)) {

    min_out <- true_bounds[1]
    max_out <- true_bounds[2]
  } else {

    min_out <- min(outcome, na.rm=T)
    max_out <- max(outcome, na.rm=T)

    print(paste0("Normalizing using the observed bounds of ",min_out, " - ",
                 max_out,". If these are incorrect, please pass the bounds to use to the true_bounds parameter."))

  }

  trans_out <- (outcome - min_out) / (max_out - min_out)

  # handle values very close to 0

  attr(trans_out, "upper_bound") <- max_out
  attr(trans_out, "lower_bound") <- min_out

  return(trans_out)


}
