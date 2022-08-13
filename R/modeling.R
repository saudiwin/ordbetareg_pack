# Main modeling functions for the package

#' Fit Ordered Beta Regression Model
#'
#' This function allows you to estimate an ordered beta regression model
#' via a formula syntax.
#'
#' This function is a wrapper around the [brms::brm] function, which is a
#' powerful Bayesian regression modeling engine using Stan. To fully explore
#' the options available, including dynamic and hierarchical modeling, please
#' see the documentation for the `brm` function above. As the ordered beta
#' regression model is currently not available in `brms` natively, this modeling
#' function allows a `brms` model to be fit with the ordered beta regression
#' distribution.
#'
#' For more information about the model, see the paper here: https://osf.io/preprints/socarxiv/2sx6y/.
#'
#' This function allows you to set priors on the dispersion parameter,
#' the cutpoints, and the regression coefficients (see below for options).
#' However, to add specific priors on individual covariates, you would need
#' to use the [brms::set_prior] function by specifying an individual covariate
#' (see function documentation) and passing the result of the function call
#' to the `extra_prior` argument.
#'
#' This function will also automatically normalize the outcome so that it
#' lies in the \\[0,1\\] interval, as required by beta regression. For furthur
#' information, see the documentation for the [normalize] function.
#'
#' @param formula Either an R formula in the form response/DV ~ var1 + var2
#'   etc. *or* formula object as created/called by the `brms`
#'   [brms::bf] function. *Please avoid using 0 or `Intercept` in the
#'   formula definition.
#' @param data An R data frame or tibble containing the variables in the formula
#' @param true_bounds If the true bounds of the outcome/response
#'   don't exist in the data, pass a length 2 numeric vector
#'   of the minimum and maximum bounds to properly normalize
#'   the outcome/response
#' @param use_brm_multiple (T/F) Whether the model should use
#'   [brms::brm_multiple] for multiple
#'   imputation over multiple dataframes passed
#'   as a list to the `data` argument
#' @param phi_reg T/F. Whether you are including a linear model predicting
#'   the dispersion  parameter, phi. Defaults to false.
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
#'  informative on the interval (0.4, 30). If the response has very low variance (i.e. tightly)
#'  clusters around a specific value, then decreasing this prior (and increasing the expected value)
#'  may be
#'  helpful. Checking the value of phi in the output of the model command
#'  will reveal if a value of 0.1 (mean of 10) is too small.
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
#'   exponential scale.
#' @param extra_prior An additional prior, such as a prior for a specific
#'  regression coefficient, added to the model by passing one of the `brms`
#'  functions [brms::set_prior] or [brms::prior_string] with appropriate
#'  values.
#' @param init This parameter is used to determine starting values for
#'   the Stan sampler to begin Markov Chain Monte Carlo sampling. It is
#'   set by default at 0 because the non-linear nature of beta regression
#'   means that it is possible to begin with extreme values depending on the
#'   scale of the covariates. Setting this to 0 helps the sampler find
#'   starting values. It does, on the other hand, limit the ability to detect
#'   convergence issues with Rhat statistics. If that is a concern, such as
#'   with an experimental feature of `brms`, set this to `"random"` to get
#'   more robust starting values (just be sure to scale the covariates so they are
#'   not too large in absolute size).
#' @param inits Alias of `init`. Overrides `init` when not NULL.
#' @param ... All other arguments passed on to the `brm` function
#' @return A `brms` object fitted with the ordered beta regression distribution.
#' @examples
#' # load survey data that comes with the package
#'
#' data("pew")
#'
#' # prepare data
#'
#' model_data <- select(pew,therm,
#'              education="F_EDUCCAT2_FINAL",
#'              region="F_CREGION_FINAL",
#'              income="F_INCOME_FINAL")
#'
#' # It takes a while to fit the models. Run the code
#' # below if you want to load a saved fitted model from the
#' # package, otherwise use the model-fitting code
#'
#' data("ord_fit_mean")
#'
#'   \donttest{
#'   # fit the actual model
#'
#'     ord_fit_mean <- ordbetareg(formula=therm ~ education + income +
#'     (1|region),
#'     data=model_data,
#'     cores=2,chains=2)
#'   }
#'
#' # access values of the coefficients
#'
#' summary(ord_fit_mean)
#' @importFrom brms brm brm_multiple
#' @importFrom brms bf
#' @export
#' 
ordbetareg <- function(formula=NULL,
                       data=NULL,
                       true_bounds=NULL,
                       phi_reg=FALSE,
                       use_brm_multiple=FALSE,
                       coef_prior_mean=0,
                       coef_prior_sd=5,
                       phi_prior=.1,
                       dirichlet_prior=c(1,1,1),
                       phi_coef_prior_mean=0,
                       phi_coef_prior_sd=5,
                       extra_prior=NULL,
                       init ="0",
                       inits = NULL
                       ...) {
  
  if(!is.null(inits)){
    init <- inits
  }
  
  if(is.null(formula)) {
    
    stop("You must pass a formula to the formula argument.")
    
  }
  
  
  if('brmsformula' %in% class(formula)) {
    
    dv <- all.vars(formula$formula)[1]
    
    formula$formula <- .update2.formula(formula$formula, "0 + Intercept ")
    
  } else if('mvbrmsformula' %in% class(formula)) {
    
    dv <- lapply(formula$forms, function(var) {
      
      all.vars(var$formula)[1]
      
    })
    
    formula$forms <- lapply(formula$forms, function(var) {
      
      if(is.null(var$family)) {
        
        var$formula <- .update2.formula(var$formula, "0 + Intercept ")
        
      }
      
      var
      
    })
    
  } else {
    
    dv <- all.vars(formula)[1]
    
    formula <- .update2.formula(formula, "0 + Intercept ")
    
  }
  
  if(is.null(data)) {
    
    stop("You must pass a data frame or tibble to the data argument.")
    
  }
  
  all_fam_types <- sapply(formula$forms, function(var) var$family$family)
  
  # figure out where it is so we can edit it
  
  if(use_brm_multiple) {
    
    if(is.data.frame(data))
      stop("To use brm_multiple with ordbetareg, please pass the multiple imputed datasets as a list to the data argument.\nMice objects are not currently supported.")
    
    if(length(dv)==1) {
      
      dv_pos <- which(names(data[[1]])==dv)
      
      
      data <- lapply(data, function(d) {
        
        if(!(all(d[[dv_pos]] >= 0 & d[[dv_pos]] <= 1,na.rm=T))) {
          
          d[[dv_pos]] <- normalize(d[[dv_pos]],true_bounds=true_bounds)
          
        } else {
          
          attr(d[[dv_pos]], "upper_bound") <- 1
          attr(d[[dv_pos]], "lower_bound") <- 0
          
        }
        
        d
        
      })
      
    } else {
      
      # multivariate adjustment necessary
      
      dv_pos <- sapply(dv, function(d) {
        
        which(names(data)==d)
        
      })
      
      data <- lapply(data, function(d) {
        
        d_prime <- lapply(1:length(d), function(c) {
          
          if(c %in% dv_pos) {
            
            if(is.null(formula$forms[[which(dv_pos==c)]]$family)) {
              
              if(!(all(d[[c]] >= 0 & d[[c]] <= 1,na.rm=T))) {
                
                out_var <- normalize(d[[c]],true_bounds=true_bounds)
                
              } else {
                
                out_var <- d[[c]]
                
                attr(out_var, "upper_bound") <- 1
                attr(out_var, "lower_bound") <- 0
                
              }
            } else {
              
              out_var <- d[[c]]
              
            }
            
          } else {
            
            out_var <- d[[c]]
            
          }
          
          return(out_var)
          
        })
        
        names(d_prime) <- names(d)
        
        d_prime
        
      })
      
      
      
      # check all outcomes
      
      
      
    }
    
    
  } else {
    
    if(length(dv)==1) {
      
      dv_pos <- which(names(data)==dv)
      
      if(!(all(data[[dv_pos]] >= 0 & data[[dv_pos]] <= 1,na.rm=T))) {
        
        data[[dv_pos]] <- normalize(data[[dv_pos]],true_bounds=true_bounds)
        
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
      
      d_prime <- lapply(1:length(data), function(c) {
        
        if(c %in% dv_pos) {
          
          if(is.null(formula$forms[[which(dv_pos==c)]]$family)) {
            
            if(!(all(data[[c]] >= 0 & data[[c]] <= 1,na.rm=T))) {
              
              out_var <- normalize(data[[c]],true_bounds=true_bounds)
              
            } else {
              
              out_var <- data[[c]]
              
              attr(out_var, "upper_bound") <- 1
              attr(out_var, "lower_bound") <- 0
              
            }
          } else {
            
            out_var <- data[[c]]
            
          }
          
        } else {
          
          out_var <- data[[c]]
          
        }
        
        return(out_var)
        
      })
      
      names(d_prime) <- names(data)
      
      data <- d_prime
      
    }
    
    
    
  }
  
  # get ordered beta regression definition
  
  
  sep_fam <- F
  suffix <- ""
  
  
  if('mvbrmsformula' %in% class(formula)) {
    # update formula objects with model families if they
    # are distinct families
    
    suffix <- paste0("_",names(all_fam_types))
    
    sep_fam <- (sum(sapply(all_fam_types, function(a) !is.null(a)))>0)
    
    if(sep_fam) {
      
      need_resp <- formula$resp[sapply(all_fam_types, is.null)]
      
    }
    
    #   if(any(!sapply(all_fam_types, is.null))) {
    #
    #     sep_fam <- T
    #     need_resp <- formula$resp[sapply(all_fam_types, is.null)]
    #     suffix <-
    #   } else {
    #     suffix <- ""
    #   }
    # } else {
    #
  }
  
  
  
  ordbeta_mod <- .load_ordbetareg(beta_prior=c(coef_prior_mean,
                                               coef_prior_sd),
                                  phireg_beta_prior = c(phi_coef_prior_mean,
                                                        phi_coef_prior_sd),
                                  dirichlet_prior_num=dirichlet_prior,
                                  phi_prior = phi_prior,
                                  extra_prior=extra_prior,
                                  suffix=suffix)
  
  if('mvbrmsformula' %in% class(formula)) {
    # update formula objects with model families if they
    # are distinct families
    
    if(sep_fam) {
      formula$forms <- lapply(formula$forms, function(var) {
        
        if(!is.null(var$family)) {
          
          return(var)
          
        } else {
          
          # add in all ordbetareg details here
          
          var <- bf(var$formula,
                    family=ordbeta_mod$family)
          
        }
        
      })
      
      # remove any ordbetareg priors for diff families
      
      all_fam_types <- ! (sapply(all_fam_types, is.null))
      
      ordbeta_mod$priors <- filter(ordbeta_mod$priors, !(grepl(x=class,
                                                               pattern="cutzero|cutone|phi") & all_fam_types[resp]))
      
      # get the outcome we need to change the priors
      
      #ordbeta_mod$priors$resp <- c(need_resp,need_resp,"",need_resp)
      
      
    }
    
  }
  
  
  
  
  
  
  
  if(use_brm_multiple) {
    
    if(phi_reg) {
      
      brm_multiple(formula=formula, data=data,
                   stanvars=ordbeta_mod$stanvars,
                   family=ordbeta_mod$family,
                   prior=ordbeta_mod$priors_phireg,
                   init=init,
                   ...)
      
    } else {
      
      if(sep_fam) {
        
        brm_multiple(formula=formula, data=data,
                     stanvars=ordbeta_mod$stanvars,
                     prior=ordbeta_mod$priors,
                     init=init,
                     ...)
        
      } else {
        
        brm_multiple(formula=formula, data=data,
                     stanvars=ordbeta_mod$stanvars,
                     family=ordbeta_mod$family,
                     prior=ordbeta_mod$priors,
                     init=init,
                     ...)
        
        
      }
      
      
      
    }
    
    
    
  } else {
    
    if(phi_reg) {
      
      brm(formula=formula, data=data,
          stanvars=ordbeta_mod$stanvars,
          family=ordbeta_mod$family,
          prior=ordbeta_mod$priors_phireg,
          init=init,
          ...)
      
    } else {
      
      if(sep_fam) {
        
        brm(formula=formula, data=data,
            stanvars=ordbeta_mod$stanvars,
            prior=ordbeta_mod$priors,
            init=init,
            ...)
        
        
      } else {
        
        brm(formula=formula, data=data,
            stanvars=ordbeta_mod$stanvars,
            family=ordbeta_mod$family,
            prior=ordbeta_mod$priors,
            init=init,
            ...)
        
        
      }
      
      
      
    }
    
  }
  
  
  
  
  
}


#' Internal Function to Add Ordered Beta Regression Family
#'
#' Not exported.
#' @importFrom brms custom_family stanvar set_prior
#' @noRd
.load_ordbetareg <- function(beta_prior=NULL,
                             phireg_beta_prior=NULL,
                             dirichlet_prior_num=NULL,
                             phi_prior=NULL,
                             extra_prior=NULL,
                             suffix="") {
  
  # function called primarily for its side effects
  
  # custom family
  
  ord_beta_reg <- custom_family("ord_beta_reg",
                                dpars=c("mu","phi","cutzero","cutone"),
                                links=c("identity","log",NA,NA),
                                lb=c(NA,0,NA,NA),
                                type="real")
  
  # stan code for density of the model
  
  stan_funs <- "

    real ord_beta_reg_lpdf(real y, real mu, real phi, real cutzero, real cutone) {

    vector[2] thresh;
    thresh[1] = cutzero;
    thresh[2] = cutzero + exp(cutone);

  if(y==0) {
      return log1m_inv_logit(mu - thresh[1]);
    } else if(y==1) {
      return log_inv_logit(mu  - thresh[2]);
    } else {
      return log_diff_exp(log_inv_logit(mu   - thresh[1]), log_inv_logit(mu - thresh[2])) +
                beta_lpdf(y|exp(log_inv_logit(mu) + log(phi)),exp(log1m_inv_logit(mu) + log(phi)));
    }
  }"
  
  stanvars <- stanvar(scode=stan_funs,block="functions")
  
  # For pulling posterior predictions
  
  posterior_predict_ord_beta_reg <- function(i, draws, ...) {
    mu <- draws$dpars$mu[, i]
    phi <- draws$dpars$phi
    cutzero <- draws$dpars$cutzero
    cutone <- draws$dpars$cutone
    N <- draws$ndraws
    
    thresh1 <- cutzero
    thresh2 <- cutzero + exp(cutone)
    
    pr_y0 <- 1 - plogis(mu - thresh1)
    pr_y1 <- plogis(mu - thresh2)
    pr_cont <- plogis(mu-thresh1) - plogis(mu - thresh2)
    out_beta <- rbeta(n=N,plogis(mu)*phi,(1-plogis(mu))*phi)
    
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
    
    low <- 1 - plogis(mu - thresh1)
    middle <- plogis(mu-thresh1) - plogis(mu - thresh2)
    high <- plogis(mu - thresh2)
    
    low*0 + middle*plogis(mu) + high
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
      out <- log(1 - plogis(mu - thresh1))
    } else if(y==1) {
      out <- log(plogis(mu - thresh2))
    } else {
      out <- log(plogis(mu-thresh1) - plogis(mu - thresh2)) + dbeta(y,plogis(mu)*phi,(1-plogis(mu))*phi,log=T)
    }
    
    out
    
  }
  
  ###### Code declaring induced dirichlet prior ####
  # code from Michael Betancourt/Staffan Betner
  # discussion here: https://discourse.mc-stan.org/t/dirichlet-prior-on-ordinal-regression-cutpoints-in-brms/20640
  dirichlet_prior <- "
  real induced_dirichlet_lpdf(real nocut, vector alpha, real phi, int cutnum, real cut1, real cut2) {
    int K = num_elements(alpha);
    vector[K-1] c = [cut1, cut1 + exp(cut2)]';
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);

    if(cutnum==1) {

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

    // divide in half for the two cutpoints

    // don't forget the ordered transformation

      return   dirichlet_lpdf(p | alpha)
           + log_determinant(J) + cut2;

    } else {

      return(0);

    }


  }

  real induced_dirichlet_rng(vector alpha, real phi, int cutnum, real cut1, real cut2) {

    int K = num_elements(alpha);
    vector[K] p;
    vector[K-1] cutpoints;

    // need to reverse the steps
    // first get the dirichlet probabilities conditional on alpha

    p = dirichlet_rng(alpha);

    // then do the *reverse* transformation to get cutpoints

    for(k in 1:(K-1)) {

       if(k==1) {

          cutpoints[k] = phi - logit(1 - p[k]);

       } else {

          cutpoints[k] = phi - logit(inv_logit(phi - cutpoints[k-1]) - p[k]);

       }

    }

    return  cutpoints[cutnum];
  }
"
  dirichlet_prior_stanvar <- stanvar(scode = dirichlet_prior, block = "functions")
  
  # stanvar(scode = "ordered[2] thresh;
  #             thresh[1] = cutzero;
  #             thresh[2] = cutzero+exp(cutone);",
  #         block = "tparameters") -> # there might be a better way to specify this
  #   dirichlet_prior_ordbeta_stanvar
  
  stanvars <- stanvars + dirichlet_prior_stanvar
  
  # Feel free to add any other priors / change the priors on b,
  # which represent regression coefficients on the logit
  # scale
  
  if(length(suffix)>1) {
    
    # multiple outcomes
    
    cutzero <- paste0("cutzero",suffix)
    cutone <- paste0("cutone",suffix)
    
    priors <- set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                            collapse=","),"]', 0, 1,", cutzero[1],",", cutone[1],")"),
                        class="cutzero",resp=substring(suffix[1],2)) +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 2,", cutzero[1],",", cutone[1],")"),
                class="cutone",resp=substring(suffix[1],2)) +
      set_prior(paste0("normal(",beta_prior[1],",",beta_prior[2],")"),class="b",resp=substring(suffix[1],2)) +
      set_prior(paste0("exponential(",phi_prior,")"),class="phi",resp=substring(suffix[1],2)) +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 1,", cutzero[2],",", cutone[2],")"),
                class="cutzero",resp=substring(suffix[2],2)) +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 2,", cutzero[2],",", cutone[2],")"),
                class="cutone",resp=substring(suffix[2],2)) +
      set_prior(paste0("normal(",beta_prior[1],",",beta_prior[2],")"),class="b",resp=substring(suffix[2],2)) +
      set_prior(paste0("exponential(",phi_prior,")"),class="phi",resp=substring(suffix[2],2))
    
    priors_phireg <- set_prior("normal(0,5)",class="b") +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 1, cutzero, cutone)"),
                class="cutzero") +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 2, cutzero, cutone)"),
                class="cutone")
    
  } else {
    
    
    cutzero <- paste0("cutzero",suffix)
    cutone <- paste0("cutone",suffix)
    
    priors <- set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                            collapse=","),"]', 0, 1,", cutzero,",", cutone,")"),
                        class="cutzero") +
      set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                    collapse=","),"]', 0, 2,", cutzero,",", cutone,")"),
                class="cutone") +
      set_prior(paste0("normal(",beta_prior[1],",",beta_prior[2],")"),class="b") +
      set_prior(paste0("exponential(",phi_prior,")"),class="phi")
    
    
  }
  
  priors_phireg <- set_prior("normal(0,5)",class="b") +
    set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                  collapse=","),"]', 0, 1, cutzero, cutone)"),
              class="cutzero") +
    set_prior(paste0("induced_dirichlet([",paste0(dirichlet_prior_num,
                                                  collapse=","),"]', 0, 2, cutzero, cutone)"),
              class="cutone")
  
  
  
  if(!is.null(extra_prior)) {
    
    priors <- priors + extra_prior
    
  }
  
  # priors <- set_prior("normal(0,5)",class="b") +
  #   prior(constant(0),class="b",coef="Intercept") +
  #   prior_string("target += normal_lpdf((cutzero + exp(cutone)) - cutzero|0,3) + cutone",check=F) +
  #   set_prior("exponential(.1)",class="phi")
  
  return(list(priors=priors,
              priors_phireg=priors_phireg,
              stanvars=stanvars,
              log_lik=log_lik_ord_beta_reg,
              posterior_epred=posterior_epred_ord_beta_reg,
              stan_funs=stan_funs,
              family=ord_beta_reg))
  
  
}

#' Helper function to add 0 + Intercept to function call\
#' @importFrom stats as.formula dbeta plogis qlogis quantile rbeta rnorm runif update var
#' @noRd
.update2.formula <- function (old, new, ...)
{
  
  # treat it like a string, break it apart on the formula sign
  
  string_form <- as.character(old)
  
  new_form <- paste0(string_form[2],string_form[1], new,
                     " + ",
                     string_form[3])
  
  as.formula(new_form)
  
}


#' Power Calculation via Simulation of the Ordered Beta Regression Model
#'
#' This function allows you to calculate power curves (or anything else)
#' via simulating the ordered beta regression model.
#'
#' This function implements the simulation found in Kubinec (2022). This
#' simulation allows you to vary the sample size, number & type of predictors,
#' values of the predictors (or treatment values), correlation between
#' predictors, and the power to target. The function returns a data frame
#' with one row per simulation draw and covariate `k`.
#' @param N The sample size for the simulation. Include a vector of integers
#'   to examine power/results for multiple sample sizes.
#' @param k The number of covariates/predictors.
#' @param iter The number of simulations to run. For power calculation,
#'   should be at least 500 (yes, this will take some time).
#' @param cores The number of cores to use to parallelize the simulation.
#' @param rho The correlation (between -1 and 1) of the predictors `k`.
#' @param phi Value of the dispersion parameter in the beta distribution.
#' @param cutpoints Value of the two cutpoints for the ordered model.
#'   By default are the values -1 and +1 (these are interpreted in the
#'   logit scale and so should not be too large). The farther apart,
#'   the fewer degenerate (0 or 1) responses there will be in the distribution.
#' @param beta_coef If not null, a vector of length `k` of the true
#'   predictor coefficients/treatment values to use for the simulation.
#'   Otherwise, coefficients are drawn from a random uniform distribution
#'   from -1 to 1 for each predictor.
#' @param beta_type Can be either `continuous` or `binary`. Use the latter
#'   for conventional treatments with two values.
#' @param treat_assign If `beta_type` is set to `binary`,
#'   you can use this parameter to set the proportion
#'   of `N` assigned to treatment. By default,
#'   the parameter is set to 0.5 for
#'   equal/balanced treatment control groups.
#' @param return_data Whether to return the simulated dqta as a list
#'   in the `data` column of the returned data frame.
#' @param seed The seed to use to make the results reproducible. Set
#'   automatically to a date-time stamp.
#' @param ... Any other arguments are passed on to the [brms::brm] function
#'   to control modeling options.
#' @return a tibble data frame with columns of simulated and estimated values and
#'   rows for each simulation iteration X coefficient combination. I.e.,
#'   if there are five predictors, and 1,000 iterations, the resulting data frame
#'   will have 1,000 rows. If there are multiple values for `N`,
#'   then each value
#'   of `N` will have its own set of iterations, making the final size of the
#'   data a multiple of the number of sample sizes to iterate over. The
#'   data frame will have the following columns:
#'   1.
#' @examples
#' # This function takes a while to run as it has
#' # to fit an ordered beta regression to each
#' # draw. The package comes with a saved
#' # simulation dataset you can inspect to see what the
#' # result looks like
#'
#' data("sim_data")
#'
#' # will take a while to run this
#' \donttest{
#'     sim_data <- sim_ordbeta(N=c(250,750),
#'     k=1,
#'     beta_coef = .5,
#'     iter=5,cores=2,
#'     beta_type="binary",
#'     treat_assign=0.3)
#'
#' }
#'
#' # to get the power values by N, simply summarize/group
#' # by N with functions from the R package dplyr
#'
#' sim_data %>%
#'   group_by(N) %>%
#'   summarize(mean_power=mean(power))
#'
#'
#' @importFrom dplyr bind_rows mutate tibble slice as_tibble arrange group_by pull summarize %>% bind_cols
#' @importFrom tidyr unchop
#' @importFrom faux rnorm_multi
#' @importFrom brms posterior_epred
#' @export
sim_ordbeta <- function(N=1000,k=5,
                        iter=1000,
                        cores=1,
                        rho=.5,
                        phi=1,
                        cutpoints=c(-1,1),
                        beta_coef=NULL,
                        beta_type='continuous',
                        treat_assign=0.5,
                        return_data=FALSE,
                        seed=as.numeric(Sys.time()),
                        ...) {
  
  # silly R package things
  
  marg_eff <- marg_eff_est <- high_marg <- low_marg <- high <- low <- x_col <- sum_stat <-  NULL
  
  set.seed(seed)
  
  if(is.null(beta_coef)) {
    
    beta_coef <- runif(k, min=-1, max=1)
    
  }
  
  if(length(beta_coef) != k) {
    
    stop("If passing in fixed beta_coef, please pass a vector of length k.")
    
  }
  
  simul_data <- lapply(1:length(N), function(n) {
    
    tibble(N=rep(N[n],iter),
           k=k,
           rho=rho,
           phi=phi,
           cutpoints1=cutpoints[1],
           cutpoints2=cutpoints[2],
           X_beta=list(beta_coef))
    
    
  }) %>% bind_rows
  
  
  print(paste0("Iterating for ", nrow(simul_data), " simulations."))
  
  # marginal effects calc
  
  eps <- 1e-7
  
  setstep <- function(x) {
    x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
  }
  
  # fit a template model
  
  print("Compiling model for re-use.")
  
  X_temp <- matrix(rep(1, k), ncol=k)
  colnames(X_temp) <- paste0(rep("Var",k),1:k)
  X_temp <- as_tibble(X_temp)
  X_temp$outcome <- .5
  
  template_mod <- ordbetareg(formula = as.formula(paste0("outcome ~ ",
                                                         paste0(rep("Var",k),1:k,
                                                                collapse=" + "))),
                             data=X_temp,
                             chains=0,cores=1,iter=100,
                             seed=seed,
                             silent=0,refresh=0,...)
  
  
  all_simul_data <- parallel::mclapply(1:nrow(simul_data), function(i) {
    
    this_data <- slice(simul_data,i)
    
    # Draw from ordered beta regression ---------------------------------------
    
    N <- this_data$N
    
    X_beta <- this_data$X_beta[[1]]
    
    # need to create X
    
    if(k>1) {
      
      X <- rnorm_multi(n=N,vars=this_data$k,r=this_data$rho,as.matrix=T)
      
    } else {
      
      X <- matrix(rnorm(n=N),ncol=1)
      
    }
    
    
    # convert to binary sampling if binary requested
    
    if(beta_type=="binary") {
      
      treat_val <- quantile(c(X),prob=treat_assign)
      
      X <- apply(X, c(1,2), function(i) {
        
        as.numeric(i<treat_val)
        
      })
      
    }
    
    
    eta <- X%*%matrix(X_beta)
    
    # ancillary parameter of beta distribution
    phi <- this_data$phi
    
    # predictor for ordered model
    mu1 <- eta[,1]
    # predictor for beta regression
    mu2 <- eta[,1]
    
    cutpoints <- c(this_data$cutpoints1,this_data$cutpoints2)
    
    # probabilities for three possible categories (0, proportion, 1)
    low <- 1-plogis(mu2 - cutpoints[1])
    middle <- plogis(mu2-cutpoints[1]) - plogis(mu2-cutpoints[2])
    high <- plogis(mu2 - cutpoints[2])
    
    # we'll assume the same eta was used to generate outcomes
    
    out_beta <- rbeta(N,plogis(mu1) * phi, (1 - plogis(mu1)) * phi)
    
    print(paste0("Now on iteration ",i))
    
    # now determine which one we get for each observation
    outcomes <- sapply(1:N, function(i) {
      sample(1:3,size=1,prob=c(low[i],middle[i],high[i]))
    })
    
    # now combine binary (0/1) with proportion (beta)
    
    final_out <- sapply(1:length(outcomes),function(i) {
      if(outcomes[i]==1) {
        return(0)
      } else if(outcomes[i]==2) {
        return(out_beta[i])
      } else {
        return(1)
      }
    })
    
    # check for floating point errors
    
    final_out <- ifelse(final_out>(1 - 1e-10) & final_out<1,final_out - 1e-10,
                        final_out)
    
    final_out <- ifelse(final_out<(0 + 1e-10) & final_out>0,final_out + 1e-10,
                        final_out)
    
    # calculate "true" marginal effects
    
    # loop over k
    
    X_low <- X
    X_high <- X
    
    marg_eff <- sapply(1:this_data$k, function(tk)  {
      
      X_low[,tk] <- X[,tk] - setstep(X[,tk])
      X_high[,tk] <- X[,tk] + setstep(X[,tk])
      
      y0 <- .predict_ordbeta(cutpoints=cutpoints,
                             X=X_low,
                             X_beta=X_beta)
      
      y1 <- .predict_ordbeta(cutpoints=cutpoints,
                             X=X_high,
                             X_beta=X_beta)
      
      mean((y1-y0)/(X_high[,tk]-X_low[,tk]))
      
      
    })
    
    
    
    
    # Fit models --------------------------------------------------------------
    
    # now fit ordinal beta
    
    # fit all models
    
    X_brms <- X
    colnames(X_brms) <- paste0(rep("Var",ncol(X)),1:ncol(X))
    X_brms <- as_tibble(X_brms)
    X_brms <- mutate(X_brms, outcome=final_out)
    
    # we only need to update this once at a time
    
    fit_model <- update(template_mod,
                        chains=1,iter=1000,
                        newdata=X_brms,recompile=F,
                        refresh=0)
    
    if('try-error' %in% class(fit_model)) {
      
      print(paste0("Estimation failed for row ",i,"\n"))
      
      this_data$status <- "estimation_failure"
      
      return(this_data)
      
    }
    
    yrep_ord <- try(posterior_epred(fit_model,draws=100))
    
    if('try-error' %in% class(yrep_ord)) {
      
      print(paste0("Posterior prediction failed for row ",i,"\n"))
      
      this_data$status <- "posterior_prediction_failure"
      
      return(this_data)
      
    }
    
    this_data$status <- "success"
    
    
    # Calculate estimands -----------------------------------------------------
    
    # now return the full data frame
    
    X_beta_ord <- as.matrix(fit_model,variable=paste0("b_Var",1:this_data$k))
    
    # calculate rmse
    
    rmse_ord <- sqrt(mean(apply(yrep_ord,1,function(c) { (c - final_out)^2 })))
    
    # calculate marginal effects
    
    cutpoints_est <- as.matrix(fit_model, variable=c('cutzero','cutone'))
    
    margin_ord <- lapply(1:ncol(X), function(tk) {
      
      X_low[,tk] <- X[,tk] - setstep(X[,tk])
      X_high[,tk] <- X[,tk] + setstep(X[,tk])
      
      tibble(marg_eff=sapply(1:nrow(X_beta_ord), function(i) {
        y0 <- .predict_ordbeta(cutpoints=cutpoints_est[i,],
                               X=X_low,
                               X_beta=c(X_beta_ord[i,]))
        
        y1 <- .predict_ordbeta(cutpoints=cutpoints_est[i,],
                               X=X_high,
                               X_beta=c(X_beta_ord[i,]))
        
        marg_eff <- (y1-y0)/(X_high[,tk]-X_low[,tk])
        
        mean(marg_eff)
      }),
      x_col=tk)
      
    }) %>% bind_rows
    
    this_data$marg_eff <- list(marg_eff)
    
    # Combine estimates -------------------------------------------------------
    
    sum_marg <- function(d,func,...) {
      
      ret_vec <-  arrange(d,x_col) %>%
        group_by(x_col) %>%
        summarize(sum_stat=func(marg_eff,...)) %>%
        pull(sum_stat)
      
      list(ret_vec)
      
    }
    
    
    out_d <- bind_cols(this_data,
                       tibble(model="Ordinal Beta Regression",
                              med_est=list(apply(X_beta_ord,2,mean)),
                              high=list(apply(X_beta_ord,2,quantile,.95)),
                              low=list(apply(X_beta_ord,2,quantile,.05)),
                              var_calc=list(apply(X_beta_ord,2,var)),
                              rmse=rmse_ord,
                              marg_eff_est=sum_marg(margin_ord,mean),
                              high_marg=sum_marg(margin_ord,quantile,.95),
                              low_marg=sum_marg(margin_ord,quantile,.05),
                              var_marg=sum_marg(margin_ord,var)))
    
    
    if(return_data) {
      
      out_d$data <- list(X_brms)
      
    }
    
    return(out_d)
    
  },mc.cores=cores) %>%
    bind_rows %>%
    unchop(c("med_est","X_beta","marg_eff","high","low","var_calc",
             'marg_eff_est','high_marg','low_marg','var_marg'))
  
  
  all_simul_data %>%
    mutate(s_err=sign(marg_eff)!=sign(marg_eff_est),
           m_err=abs(marg_eff_est)/abs(marg_eff),
           coverage=ifelse(marg_eff>0,marg_eff<high_marg & marg_eff>low_marg,
                           marg_eff<high_marg & marg_eff>low_marg),
           power=as.numeric(ifelse(sign(marg_eff)==sign(high) & sign(marg_eff)==sign(low),
                                   1,
                                   0)),
           treat_assign=treat_assign,
           beta_type=beta_type,
           cutpoints=list(cutpoints),
           seed=seed)
  
}


#' Internal function for calculating predicting ordered beta for simulation
#' @noRd
.predict_ordbeta <- function(cutpoints=NULL,X=NULL,X_beta=NULL,
                             combined_out=T) {
  
  # we'll assume the same eta was used to generate outcomes
  eta <- X%*%matrix(X_beta)[,1]
  
  # probabilities for three possible categories (0, proportion, 1)
  low <- 1-plogis(eta - cutpoints[1])
  middle <- plogis(eta-cutpoints[1]) - plogis(eta-cutpoints[2])
  high <- plogis(eta - cutpoints[2])
  
  # check for whether combined outcome or single outcome
  
  if(combined_out) {
    low*0 + middle*plogis(eta) + high*1
  } else {
    list(pr_zero=low,
         pr_proportion=middle,
         proportion_value=plogis(eta),
         pr_one=high)
  }
  
}
