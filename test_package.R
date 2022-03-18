data("pew")

# prepare data

model_data <- select(pew,therm,
                     education="F_EDUCCAT2_FINAL",
                     region="F_CREGION_FINAL",
                     income="F_INCOME_FINAL",
                     THERMO_THERMBA_W28) %>%
  mutate(THERMO_THERMBA_W28=as.numeric(THERMO_THERMBA_W28))


# fit the actual model

ord_fit_one_mod <- ordbetareg(formula=bf(THERMO_THERMBA_W28~ education),
                               data=model_data,
                               cores=2,chains=2,sample_prior="only",
                               backend="cmdstanr")

ord_fit_one_mod_newbounds <- ordbetareg(formula=bf(THERMO_THERMBA_W28~ education),
                              data=model_data,
                              true_bounds=c(0,100),
                              cores=2,chains=2,sample_prior="only",
                              backend="cmdstanr")

ord_fit_two_mods <- ordbetareg(formula=bf(THERMO_THERMBA_W28~ education) +
                                 bf(therm ~ region),
                               data=model_data,
                               cores=2,chains=2,sample_prior="only",
                               backend="cmdstanr")

ord_fit_two_mods_dfam <- ordbetareg(formula=bf(THERMO_THERMBA_W28~ education, family="gaussian") +
                             bf(therm ~ region),
                           data=model_data,use_brm_multiple = T,
                           cores=2,chains=2,sample_prior="only",
                           backend="cmdstanr")
ord_fit_two_mods_mi <- ordbetareg(formula=bf(THERMO_THERMBA_W28~ education) +
                                 bf(therm ~ region),
                               data=list(model_data,model_data),
                               use_brm_multiple = T,
                               cores=2,chains=2,sample_prior="only",
                               backend="cmdstanr")

ord_fit_two_mods_dfam_mi <- ordbetareg(formula=bf(THERMO_THERMBA_W28~ education, family="gaussian") +
                                      bf(therm ~ region),
                                    data=list(model_data,model_data),use_brm_multiple = T,
                                    cores=2,chains=2,sample_prior="only",
                                    backend="cmdstanr")

# check with cmdstnar

require(cmdstanr)

dirich <- cmdstan_model("~/test.stan",quiet=F,compile = T,force_recompile=T)

dirich_fit <- dirich$sample(data=list(K=3))
