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
uncertainty = 0.5

sim <- crossing(
  participant = 1:n
) %>%
  mutate(
    ability = rnorm(n(), 0.5, sd = 0.20),
    confidence = rnorm(n(), 0, sd=uncertainty),
    perceived_accurate=ability + confidence,
    perceived=in_range(perceived_accurate,0, 1),
    actual = ability %>% in_range(0, 1),
    bias = perceived - actual,
    real_bias = perceived_accurate - ability
  )


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
