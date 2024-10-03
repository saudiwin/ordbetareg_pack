# calculate Pr(Y=1) and Pr(Y=0) (i.e. the discrete categories of the scale)
# for different groups in the data
# useful for treatment/control comparisons

library(ordbetareg)
library(dplyr)
library(tidyr)
library(posterior)
library(tidybayes)
library(ggplot2)

data("pew")

# prepare data

model_data <- select(pew,therm,
                     education="F_EDUCCAT2_FINAL",
                     region="F_CREGION_FINAL",
                     income="F_INCOME_FINAL")

# we will look at probabilities of y=1 or y=0 for different levels of education

table(model_data$education)

ord_fit_mean <- ordbetareg(formula=therm ~ education + income +
                             (1|region),
                           data=model_data,
                           cores=2,chains=2)

# function that will return tibble with Pr(y=1) and Pr(y=0)

return_pr_01 <- function(object) {

  # get linear prediction

  ord_pred <- add_linpred_draws(newdata=object$data,
                                object=object,ndraws=200)

  # need top and bottom cutpoints

  cut_draws <- spread_draws(object,cutzero,cutone)

  # merge with linear prediction

  ord_pred <- left_join(ord_pred,
                        distinct(select(cut_draws,.draw,
                               cutzero,cutone)), by=c(".draw"))

  # calculate pr_0 and pr_1

  ord_pred <- group_by(ord_pred, .iteration) %>%
                     mutate(pr_y_1=plogis(.linpred - cutone),
                            pr_y_0=1 - plogis(.linpred - cutzero))

  return(ungroup(ord_pred))

}

# what we want is the difference in Pr(Y=1) between two groups in the data:
# education = Less than high school and education = Postgraduate
# collapse to this quantity for each posterior draw (group by .draw)

ord_pry1y0 <- return_pr_01(ord_fit_mean)

change_groups_y1 <- group_by(ord_pry1y0, .draw) %>%
  summarize(group_diff=mean(pr_y_1[education=="Postgraduate"]) - mean(pr_y_1[education=="Less than high school"]))

# look at the estimate of the difference in terms of posterior draws:

hist(change_groups_y1$group_diff)
summary(change_groups_y1$group_diff)

# difference is mostly positive (more 1s in the postgraduate group relative to low-education group)

 # test with 5% - 95% interval

quantile(change_groups_y1$group_diff, prob=c(.05, .95))

# this interval is only positive values

# also plot with ggplot & tidybayes

change_groups_y1 %>%
  ggplot(aes(x=group_diff)) +
  stat_halfeye()

# now the same for Pr(y=0)

change_groups_y0 <- group_by(ord_pry1y0, .draw) %>%
  summarize(group_diff=mean(pr_y_0[education=="Postgraduate"]) - mean(pr_y_0[education=="Less than high school"]))

# look at the estimate of the difference in terms of posterior draws:

hist(change_groups_y0$group_diff)
summary(change_groups_y0$group_diff)

# difference is more negative (more 0s in the postgraduate group relative to low-education group)

# test with 5% - 95% interval

quantile(change_groups_y0$group_diff, prob=c(.05, .95))

# this interval includes only negative values

# also plot with ggplot & tidybayes

change_groups_y0 %>%
  ggplot(aes(x=group_diff)) +
  stat_halfeye()

# combined plot

bind_rows(list(`Pr(Y=1)`=change_groups_y1,
               `Pr(Y=0)`=change_groups_y0),
          .id="Outcome Type") %>%
  ggplot(aes(x=group_diff)) +
  stat_halfeye(aes(fill=`Outcome Type`),
               alpha=0.5) +
  labs(x="Difference in Probability of Outcome from Low to High Education",
       y="")





