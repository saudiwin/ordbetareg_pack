library(tidyverse)
library(dagitty)
library(ggdag)

dag_kruger <- dagitty("dag{ A -> P;
                         C -> P;
                         P -> S;
                         A -> S;
                         P -> Ahat;
                         A -> Phat;
                         S -> Phat;
                         S -> Ahat;
                         S -> B
                      }")

theme_set(theme_minimal())

thematic::auto_config_set(
  thematic::auto_config(
    bg = "#D8DEE9",
    fg = "#2E3440",
    accent = "#5E81AC"
  )
)
thematic::thematic_on()

in_range <- function(x, lower, upper) {
  case_when(
    x > upper ~ upper,
    x < lower ~ lower,
    TRUE      ~ x
  )
}

n = 10000
uncertainty = 0.1

# change this so that perceived - ability, but reported is censored.

over_unc <- lapply(seq(0.01,3, by=0.01), function(u) {

  sim <- crossing(
    participant = 1:n
  ) %>%
    mutate(
      ability = runif(n(), 0, 1),
      ability2 = 1.89*ability,
      test_score=ability + rnorm(n(), 0, sd=u),
      perceived_ab=ability + rnorm(n(), 0, sd=u),
      # confidence = rnorm(n(), 0, sd=uncertainty),
      # perceived_accurate=ability + confidence,
      perceived_ab=plogis(ability),
      perceived_test=plogis(test_score),
      # actual = ability %>% in_range(0, 1),
      # bias = perceived - actual,
      # real_bias = perceived_accurate - ability
    )

  tibble(unc=u,
         est=lm(perceived_test ~ perceived_ab, data=sim)$coefficients[2])


}) %>% bind_rows

plot(over_unc$unc, over_unc$est)

sim <- crossing(
  participant = 1:n
) %>%
  mutate(
    ability = rbeta(n(),1,1),
    ability2 = 1.89*ability,
    step1 = qlogis(ability),
    err=rnorm(n(),sd=1),
    perceived_ab=,
    # confidence = rnorm(n(), 0, sd=uncertainty),
    # perceived_accurate=ability + confidence,
    # actual = ability %>% in_range(0, 1),
    # bias = perceived - actual,
    # real_bias = perceived_accurate - ability
  )

this_sim <- sim_ordbeta(N=1000, k=1, beta_coef=3,return_data=T,iter=1,cutpoints = c(-2,2))

m1 <- ordbetareg(outcome ~ Var1 + I(Var1^2) + I(Var1^3), data=this_sim$data[[1]],backend="cmdstanr",
                 threads=threading(16),chains=1)
summary(m1)
summary(lm(outcome ~ Var1 + I(Var1^2) + I(Var1^3), data=this_sim$data[[1]]))

ggplot(this_sim$data[[1]], aes(y=outcome, x=Var1))  + geom_point() + stat_smooth(method="lm") +
  geom_abline(slope=1, intercept=0, linetype=2, colour="red")
ggplot(this_sim$data[[1]], aes(y=outcome, x=Var1))  + geom_point() + stat_smooth() +
  geom_abline(slope=1, intercept=0, linetype=2, colour="red")
ggplot(this_sim$data[[1]], aes(y=outcome, x=Var1))  + geom_point() + stat_smooth(method = "lm",formula = y ~ poly(x,2)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="red")



sim_quantiles <- sim %>%
  summarise(
    quantile_breaks = quantile(actual),
    quantile = names(quantile_breaks),
    percieved_quantile_breaks = quantile(perceived),
    percieved_quantile = names(percieved_quantile_breaks)
  )

# sim <- sim %>%
#   mutate(
#     actual_q = cut(actual,
#                    breaks = sim_quantiles$quantile_breaks,
#                    labels = sim_quantiles$quantile[-1],
#                    include.lowest = TRUE),
#     perceived_q = cut(perceived,
#                       breaks = sim_quantiles$percieved_quantile_breaks,
#                       labels = sim_quantiles$percieved_quantile[-1],
#                       include.lowest = TRUE)
#   )
#
# sim_q_summary <- sim %>%
#   group_by(actual_q) %>%
#   summarise(across(c(actual, perceived),
#                    mean)) %>%
#   mutate(quartile_mean = actual) %>%
#   pivot_longer(c(actual, perceived))

# model this with ordered beta regression

library(ordbetareg)
library(brms)

sim_perc <- mutate(sim, perc=cut(actual,
                                 breaks=quantile(actual, prob=seq(0,1,by=0.01)),
                                 labels=names(quantile(actual, prob=seq(0,1,by=0.01)))[-1],
                                  include.lowest=TRUE),
                   perc=as.numeric(stringr::str_remove(perc,"%")),
                   bias=actual-perceived)

est_mod <- ordbetareg(formula=perceived ~ actual + I(actual^2) +
                        I(actual^3),data=sim_perc,
                      backend="cmdstanr",chains=1,
                      threads=threading(10))

# est_mod <- ordbetareg(formula=perceived_actual ~ actual,data=sim_perc,
#                       backend="cmdstanr",chains=1,
#                       threads=threading(10))

summary(est_mod)

check_mod <- lm(bias ~ ability*confidence, data=sim_perc)

summary(check_mod)

summary(lm(perceived ~ actual + I(actual^2) + I(actual^3), data=sim_perc))

summary(lm(confidence ~ ability + perceived, data=sim_perc))

summary(lm(confidence ~ ability + perceived_accurate, data=sim_perc))

summary(lm(perceived_accurate ~ ability*confidence, data=sim_perc))

summary(lm(perceived ~ actual + ability, data=sim_perc))
