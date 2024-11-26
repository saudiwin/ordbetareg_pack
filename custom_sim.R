# custom sim with rordbeta

library(ordbetareg)
library(glmmTMB)
library(tidyverse)

# we'll use glmmTMB because it's much faster and we're only going to use it
# for point estimation

# parameters to loop over and get power curves for

treat_effects <- seq(0,1, by=0.1)
num_subjects <- seq(10,20,by=1)
num_measures_per_day <- seq(1,5,by=1)
num_days <- 10
# measure of how much variance there is within-subject over time
# equivalent to sd of random intercepts
sd_within_subject <- .2

# need to simulate dispersion (phi)
# higher values = more concentration (sort of like lower variance)

phi <- 1

# p-value threshold to use

p_threshold <- 0.05

# make a big grid

all_params <- expand.grid(treat_effects,num_subjects,
                          num_measures_per_day,
                          num_days,
                          sd_within_subject,
                          phi)

names(all_params) <- c("treat_effects","num_subjects",
                       "num_measures_per_day","num_days","sd_within_subject",
                       "phi")

# number of sims per parameter combination

sims <- 100

# parallel processing to make it run faster
# set to just vanilla lapply if on Windows

over_params <- parallel::mclapply(1:nrow(all_params), function(i) {

  this_param <- slice(all_params, i)

  # now loop over iterations

  over_sims <- lapply(1:sims, function (s) {

    # first draw random intercept for each respondent

    var_int <- rnorm(n=this_param$num_subjects,
                     mean=0,sd=this_param$sd_within_subject)

    # make vector with one intercept per participant
    # multiple measures = multiple obs per respondent

    var_int_all <- rep(var_int,each=this_param$num_measures_per_day*this_param$num_days)

    # simulate outcome with function rordbeta given parameters

    # first simulate predictor as random uniform

    covar <- runif(n=this_param$num_subjects * this_param$num_measures_per_day * this_param$num_days,
                   min = 0,
                   max=1)

    # note: no overall intercept, we assume it is zero but that could be changed

    # logit function for linear model

    linear_model <- plogis(this_param$treat_effects * covar +
                             var_int_all)

    out_ordbeta <- rordbeta(n=this_param$num_subjects * this_param$num_measures_per_day * this_param$num_days,
                            mu = linear_model,
                            phi=this_param$phi)

    # fit a glmmTMB model for speed

    to_model <- tibble(out_ordbeta=out_ordbeta,
                       covar=covar,
                       subject=rep(1:this_param$num_subjects,
                                   each=this_param$num_measures_per_day * this_param$num_days),
                       days=rep(1:this_param$num_days, times=this_param$num_measures_per_day*this_param$num_subjects))

    # fit varying intercepts model if > 1 obs per respondent

    if(this_param$num_measures_per_day>1 || this_param$num_days>1) {

      glmtmb_fit <- glmmTMB(out_ordbeta ~ covar + (1|subject),data = to_model,
                            family=ordbeta)

    } else {

      glmtmb_fit <- glmmTMB(out_ordbeta ~ covar,data = to_model,
                            family=ordbeta)

    }

    # now get estimates and see if coef is significant and what the bias is

    glm_sum <- summary(glmtmb_fit)

    # simulation results

    sim_res <- tibble(param_vals=i,
                      iteration=s,
                      treat_est=glm_sum$coefficients$cond['covar','Estimate'],
                      treat_pvalue=glm_sum$coefficients$cond['covar','Pr(>|z|)'],
                      treat_err=glm_sum$coefficients$cond['covar','Std. Error'],
                      true_treat=this_param$treat_effects,
                      true_phi=this_param$phi,
                      num_subjects=this_param$num_subjects,
                      num_measures_per_day=this_param$num_measures_per_day,
                      sd_within_subject=this_param$sd_within_subject,
                      p_threshold=p_threshold)

    # power = p < threshold

    sim_res <- mutate(sim_res,
                      treat_sig=as.numeric(treat_pvalue > p_threshold))

    return(sim_res)



  }) %>% bind_rows

  return(over_sims)

},mc.cores=parallel::detectCores()) %>% bind_rows

# estimate power by param combination

power_est <- group_by(over_params, param_val,true_treat,
                      true_phi,num_subjects, num_measures_per_day,
                      sd_within_subject,p_threshold) %>%
  summarize(power_est=mean(treat_sig))

# plot power curve for varying N
# with different lines for different true treatment effects
# facet wrap by num_measures_per_day

library(ggplot2)

power_est %>%
  ggplot(aes(y=power_est,
             x=num_subjects)) +
  geom_line(aes(colour=true_treat)) +
  facet_wrap(~num_measures_per_day)
