library(palmerpenguins)
library(tidyverse)
library(brms)
library(ordbetareg)
library(mice)

penguins_mice <- penguins |>
  mutate(body_mass_g = ordbetareg::normalize(body_mass_g))


obr <- ordbetareg(formula = bf(body_mass_g ~ sex + species + island,
                               cutone ~ sex,
                               cutzero ~ sex),
                              iter = 500,
                              backend = "cmdstanr",
                              silent = 0,
                              cores = 8,
                              chains = 2,
                  manual_prior = set_prior("normal(0,1)") +
                    set_prior("normal(0,1)", class="b",dpar="cutone") +
                    set_prior("normal(0,1)", class="b",dpar="cutzero"),
                              data = penguins_mice)

# see how it works by looking at predictions by sex

plot_heiss(obr, grouping_fac = "sex")

# ok not a huge difference, haha, but you get the idea
