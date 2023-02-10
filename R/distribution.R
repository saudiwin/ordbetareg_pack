# functions that define the ordered beta statistical distribution



#' Probability Density Function for the Ordered Beta Distribution
#'
#'
#' This function will return the density of given variates of the
#' ordered beta distribution conditional on values for
#' the mean (`mu`), dispersion (`phi`) and cutpoints
#' governing the ratio of degenerate (discrete) to continuous
#' responses.
#'
#' @export
#' @param x Variates of the ordered beta distribution (should be in the \[0,1\] interval).
#' @param mu Value of the mean of the distribution.
#' Should be in the \(0,1\) interval (cannot be strictly equal to 0 or 1). If
#' length is greater than 1, should be of length x.
#' @param phi Value of the dispersion parameter. Should be strictly greater than 0. If
#' length is greater than 1, should be of length x.
#' @param cutpoints A vector of two numeric values for the cutpoints. Second value should
#' @param log where to return the log density
#' be strictly greater than the first value.
#' @examples
#'
#' # examine density (likelihood) of different possible values
#' # given fixed values for ordered beta parameters
#'
#' x <- seq(0, 1, by=0.01)
#'
#' x_dens <- dordbeta(x, mu = 0.3, phi=2, cutpoints=c(-2, 2))
#'
#' # Most likely value for x is approx 1
#' # Note discontinuity in density function between continuous/discrete values
#' # density function is a combined PMF/PDF, so not a real PDF
#' # can though be used for MLE
#'
#' plot(x_dens, x)
#'
#' # discrete values should be compared to each other:
#' # prob of discrete 0 > prob of discrete 1
#'
#' x_dens[x==0] > x_dens[x==1]
dordbeta <- function(x=.9,
                     mu=0.5,
                     phi=1,
                     cutpoints=c(-1,1),
                     log=FALSE) {

  if(!all(mu>0 & mu<1) ) {

    stop("Please pass a numeric value for mu that is between 0 and 1.")

  }

  if(!all(phi>0)) {

    stop("Please pass a numeric value for phi that is greater than 0.")

  }

  if(!(length(mu) %in% c(1,length(x)))) {

    stop("Please pass a vector for mu that is either length 1 or length N.")

  }

  if(!(length(phi) %in% c(1,length(x)))) {

    stop("Please pass a vector for phi that is either length 1 or length N.")

  }

  mu_ql <- qlogis(mu)

  if(length(mu_ql)==1) {
    mu_ql <- rep(mu_ql, length(x))
  }

  # probabilities for three possible categories (0, proportion, 1)
  low <- 1-plogis(mu_ql - cutpoints[1])
  middle <- plogis(mu_ql - cutpoints[1]) - plogis(mu_ql - cutpoints[2])
  high <- plogis(mu_ql - cutpoints[2])

  # we'll assume the same eta was used to generate outcomes

  out_beta <- dbeta(x = x,mu * phi, (1 - mu) * phi)

  # now determine which one we get for each observation
  outcomes <- sapply(1:length(x), function(i) {

      if(x[i]==0) {

        out <- low[i]

      } else if(x[i]==1) {

        out <- high[i]

      } else {

        out <- middle[i] * out_beta[i]

      }

    return(out)

  })

  if(log) {

    return(log(outcomes))

  } else {

    return(outcomes)

  }

}

#' Generate Ordered Beta Variates
#'
#'
#' This function will generate ordered beta random variates given
#' values for the mean (`mu`), dispersion (`phi`) and cutpoints
#' governing the ratio of degenerate (discrete) to continuous
#' responses.
#'
#' @export
#' @param n Number of variates to generate.
#' @param mu Value of the mean of the distribution.
#' Should be in the \(0,1\) interval (cannot be strictly equal to 0 or 1). If
#' length is greater than 1, should be of length `n`.
#' @param phi Value of the dispersion parameter. Should be strictly greater than 0. If
#' length is greater than 1, should be of length `n`.
#' @param cutpoints A vector of two numeric values for the cutpoints. Second value should
#' be strictly greater than the first value.
#' @examples
#'
#' # generate 100 random variates with an average of 0.7
#' # all will be in the closed interval \[0,1\]
#'
#' ordbeta_var <- rordbeta(n=100, mu=0.7, phi=2)
#'
#' # Will be approx mean = 0.7 with high positive skew
#'
#' summary(ordbeta_var)
rordbeta <- function(n=100,
                     mu=0.5,
                     phi=1,
                     cutpoints=c(-1,1)) {

  if(!all(mu>0 & mu<1) ) {

    stop("Please pass a numeric value for mu that is between 0 and 1.")

  }

  if(!all(phi>0)) {

    stop("Please pass a numeric value for phi that is greater than 0.")

  }

  if(!(length(mu) %in% c(1,n))) {

    stop("Please pass a vector for mu that is either length 1 or length N.")

  }

  if(!(length(phi) %in% c(1,n))) {

    stop("Please pass a vector for phi that is either length 1 or length N.")

  }

  mu_ql <- qlogis(mu)

  if(length(mu_ql)==1) {
    mu_ql <- rep(mu_ql, n)
  }

  # probabilities for three possible categories (0, proportion, 1)
  low <- 1-plogis(mu_ql - cutpoints[1])
  middle <- plogis(mu_ql - cutpoints[1]) - plogis(mu_ql - cutpoints[2])
  high <- plogis(mu_ql - cutpoints[2])

  # we'll assume the same eta was used to generate outcomes

  out_beta <- rbeta(n = n,mu * phi, (1 - mu) * phi)

  # now determine which one we get for each observation
  outcomes <- sapply(1:n, function(i) {

      sample(1:3,size=1,prob=c(low[i],middle[i],high[i]))

  })

  # now combine binary (0/1) with proportion (beta)

  final_out <- sapply(1:n,function(i) {
    if(outcomes[i]==1) {
      return(0)
    } else if(outcomes[i]==2) {
      return(out_beta[i])
    } else {
      return(1)
    }
  })


  return(final_out)

}
