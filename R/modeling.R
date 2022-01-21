# Main modeling functions for the package

#' Fit Ordered Beta Regression Model
#'
#' This function allows you to estimate an ordered beta regression model
#' via a formula syntax.
#'
#' This function is a wrapper around the [brms::brm] function, which is a
#' powerful Bayesian regresion modeling engine using Stan. To fully explore
#' the options available, including dynamic and hierarchical modeling, please
#' see the documentation for the `brm` function above. As the ordered beta
#' regression model is currently not available in `brms` natively, this modeling
#' function allows a `brms` model to be fit with the ordered beta regression
#' distribution.
#'
#' This function will also automatically normalize the outcome so that it
#' lies in the [0,1] interval, as required by beta regression. For furthur
#' information, see the documentation for the [normalize] function.
#'
#' @param formula Either an R formula in the form response/DV ~ var1 + var2
#'   etc. *or* formula object as created/called by the `brms`
#'   [brms::bf] function. *Please avoid using 0 or `Intercept` in the
#'   formula definition.
#' @param data An R data frame or tibble containing the variables in the formula
#' @param coef_prior_mean The mean of the Normal distribution prior on the
#'   regression coefficients (for predicting the mean of the response).
#'   Default is 0.
#' @param coef_prior_sd The SD of the Normal distribution prior on the
#'   regression coefficients (for predicting the mean of the response).
#'   Default is 5, which makes the prior weakly informative on the
#'   logit scale.
#' @param phi_prior The mean parameter of the exponential prior on
#'  phi, which determines the dispersion of the beta distribution. The
#'  default is .1, which equals a mean of 10 and is thus weakly
#'  informativce. If the response has very low variance (i.e. tightly)
#'  clusters around a specific value, then increasing this prior may be
#'  helpful. Checking the value of phi in the output of the model command
#'  will reveal if 10 is too small.
#' @param dirichlet_prior A vector of three integers
#'   corresponding to the prior parameters for the dirchlet distribution
#'   (alpha parameter) governing the location of the cutpoints between
#'   the components of the response (continuous vs. degenerate).
#'   The default is 1 which puts equal probability on
#'   degenerate versus continuous responses. Likely only needs to be
#'   changed in a repeated sampling situation to stabilize the cutpoint
#'   locations across samples.
#' @param phi_coef_prior_mean The mean of the Normal distribution prior on the
#'   regression coefficients for predicting phi, the dispersion parameter.
#'   Only useful if a linear model is being fit to phi.
#'   Default is 0.
#' @param phi_coef_prior_sd The SD of the Normal distribution prior on the
#'   regression coefficients for predicting phi, the dispersion parameter.
#'   Only useful if a linear model is being fit to phi.
#'   Default is 5, which makes the prior weakly informative on the
#'   logit scale.
#' @param extra_prior An additional prior, such as a prior for a specific
#'  regression coefficient, added to the model by passing one of the `brms`
#'  functions [brms::set_prior] or [brms::prior_string] with appropriate
#'  values.
#' @param ... All other arguments passed on to the `brm` function
#' @importFrom brms brm
#' @importForm brms bf
#' @export
ordbetareg <- function(formula=NULL,
                       data=NULL,
                       coef_prior_mean=0,
                       coef_prior_sd=5,
                       phi_prior=.1,
                       dirichlet_prior=c(1,1,1),
                       phi_coef_prior_mean=0,
                       phi_coef_prior_sd=5,
                       ...) {

  if(is.null(formula)) {

    stop("You must pass a formula to the formula argument.")

  }


  if('brmsformula' %in% class(formula)) {

    dv <- labels(terms(formula$formula))[1]

  } else if('mvbrmsformula' %in% class(formula)) {

    dv <- lapply(formula$forms, function(var) {

          labels(terms(var$formula))[1]

      })

  } else {

    dv <- labels(terms(formula))[1]

  }

  if(is.null(data)) {

    stop("You must pass a data frame or tibble to the data argument.")

  }

  # figure out where it is so we can edit it

  if(length(dv)==1) {

    dv_pos <- which(names(data)==dv)

    if(!(all(data[[dv_pos]] >= 0 & data[[dv_pos]] <= 1,na.rm=T))) {

      data[[dv_pos]] <- normalize(data[[dv_pos]])

    } else {

      attr(data[[dv_pos]], "upper_bound") <- 1
      attr(data[[dv_pos]], "lower_bound") <- 0

    }

  } else {

    # multivariate adjustment necessary

    dv_pos <- sapply(dv, function(d) {

      which(names(data)==d)

    })

    # check all outcomes

    data <- lapply(1:length(data), function(c) {

      if(c %in% dv_pos) {

        if(!(all(data[[c]] >= 0 & data[[c]] <= 1,na.rm=T))) {

          out_var <- normalize(data[[c]])

        } else {

          out_var <- data[[c]]

          attr(out_var, "upper_bound") <- 1
          attr(out_var, "lower_bound") <- 0

        }

      } else {

        out_var <- data[[c]]

      }

      return(out_var)

    })

  }


  # get ordered beta regression definition

  ordbeta_mod <- .load_ordbetareg(beta_prior=c(coef_prior_mean,
                                               coef_prior_sd),
                                  phireg_beta_prior = c(phi_coef_prior_mean,
                                                        phi_coef_prior_sd),
                                  dirichlet_prior=dirichlet_prior,
                                  phi_prior = phi_prior,
                                  extra_prior=extra_prior)


  brm(formula=formula, data=data,
      stanvars=ordbeta_mod$stan_vars,
      family=ordbeta_mod$family,
      priors=ordbeta_mod$priors,...)

}


#' Internal Function to Add Ordered Beta Regression Family
#'
#' Not exported.
#' @noRd
.load_ordbetareg <- function(beta_prior=NULL,
                             phireg_beta_prior=NULL,
                             dirichlet_prior=NULL,
                             phi_prior=NULL,
                             extra_prior=NULL) {

  # function called primarily for its side effects

  # custom family

  ord_beta_reg <- custom_family("ord_beta_reg",
                                dpars=c("mu","phi","cutzero","cutone"),
                                links=c("logit","log",NA,NA),
                                lb=c(NA,0,NA,NA),
                                type="real")

  # stan code for density of the model

  stan_funs <- "real ord_beta_reg_lpdf(real y, real mu, real phi, real cutzero, real cutone) {

    //auxiliary variables
    real mu_logit = logit(mu);
    vector[2] thresh;
    thresh[1] = cutzero;
    thresh[2] = cutzero + exp(cutone);

  if(y==0) {
      return log1m_inv_logit(mu_logit - thresh[1]);
    } else if(y==1) {
      return log_inv_logit(mu_logit  - thresh[2]);
    } else {
      return log(inv_logit(mu_logit  - thresh[1]) - inv_logit(mu_logit - thresh[2])) +
                beta_proportion_lpdf(y|mu,phi);
    }
  }"

  stanvars <- stanvar(scode=stan_funs,block="functions")

  # For pulling posterior predictions

  posterior_predict_ord_beta_reg <- function(i, draws, ...) {
    mu <- draws$dpars$mu[, i]
    phi <- draws$dpars$phi
    cutzero <- draws$dpars$cutzero
    cutone <- draws$dpars$cutone
    N <- length(phi)

    thresh1 <- cutzero
    thresh2 <- cutzero + exp(cutone)

    pr_y0 <- 1 - plogis(qlogis(mu) - thresh1)
    pr_y1 <- plogis(qlogis(mu) - thresh2)
    pr_cont <- plogis(qlogis(mu)-thresh1) - plogis(qlogis(mu) - thresh2)
    out_beta <- rbeta(n=N,mu*phi,(1-mu)*phi)

    # now determine which one we get for each observation
    outcomes <- sapply(1:N, function(i) {
      sample(1:3,size=1,prob=c(pr_y0[i],pr_cont[i],pr_y1[i]))
    })

    final_out <- sapply(1:length(outcomes),function(i) {
      if(outcomes[i]==1) {
        return(0)
      } else if(outcomes[i]==2) {
        return(out_beta[i])
      } else {
        return(1)
      }
    })

    final_out

  }

  # for calculating marginal effects/conditional expectations

  posterior_epred_ord_beta_reg<- function(draws) {

    cutzero <- draws$dpars$cutzero
    cutone <- draws$dpars$cutone

    mu <- draws$dpars$mu

    thresh1 <- cutzero
    thresh2 <- cutzero + exp(cutone)

    low <- 1 - plogis(qlogis(mu) - thresh1)
    middle <- plogis(qlogis(mu)-thresh1) - plogis(qlogis(mu) - thresh2)
    high <- plogis(qlogis(mu) - thresh2)

    low*0 + middle*mu + high
  }

  # for calcuating LOO and Bayes Factors

  log_lik_ord_beta_reg <- function(i, draws) {

    mu <- draws$dpars$mu[,i]
    phi <- draws$dpars$phi
    y <- draws$data$Y[i]
    cutzero <- draws$dpars$cutzero
    cutone <- draws$dpars$cutone

    thresh1 <- cutzero
    thresh2 <- cutzero + exp(cutone)

    if(y==0) {
      out <- log(1 - plogis(qlogis(mu) - thresh1))
    } else if(y==1) {
      out <- log(plogis(qlogis(mu) - thresh2))
    } else {
      out <- log(plogis(qlogis(mu)-thresh1) - plogis(qlogis(mu) - thresh2)) + dbeta(y,mu*phi,(1-mu)*phi,log=T)
    }

    out

  }

  ###### Code declaring induced dirichlet prior ####
  # code from Michael Betancourt/Staffan Betner
  # discussion here: https://discourse.mc-stan.org/t/dirichlet-prior-on-ordinal-regression-cutpoints-in-brms/20640
  dirichlet_prior <- "
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);

    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];

    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;

    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }

    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
"
  dirichlet_prior_stanvar <- stanvar(scode = dirichlet_prior, block = "functions")

  stanvar(scode = "ordered[2] thresh;
              thresh[1] = cutzero;
              thresh[2] = cutzero+exp(cutone);",
          block = "tparameters") -> # there might be a better way to specify this
    dirichlet_prior_ordbeta_stanvar

  stanvars <- stanvars + dirichlet_prior_stanvar + dirichlet_prior_ordbeta_stanvar

  # Feel free to add any other priors / change the priors on b,
  # which represent regression coefficients on the logit
  # scale

  priors <- set_prior("target += induced_dirichlet_lpdf(thresh | rep_vector(1, 3), 0)", check=FALSE) +
    set_prior(paste0("normal(",beta_prior[1],",",beta_prior[2],")"),class="b") +
    set_prior(paste0("exponential(",phi_prior,")"),class="phi")

  if(!is.null(extra_prior)) {

    priors <- priors + extra_prior

  }

  # priors <- set_prior("normal(0,5)",class="b") +
  #   prior(constant(0),class="b",coef="Intercept") +
  #   prior_string("target += normal_lpdf((cutzero + exp(cutone)) - cutzero|0,3) + cutone",check=F) +
  #   set_prior("exponential(.1)",class="phi")

  priors_phireg <- set_prior("normal(0,5)",class="b") +
    set_prior("target += induced_dirichlet_lpdf(thresh | rep_vector(1, 3), 0)", check=FALSE)

  return(list(priors=priors,
              priors_phireg=priors_phireg,
              stanvars=stanvars,
              log_lik=log_lik_ord_beta_reg,
              posterior_epred=posterior_epred_ord_beta_reg,
              stan_funs=stan_funs,
              family=ord_beta_reg))


}


#' Function to Simulate the Ordered Beta Regression Model
