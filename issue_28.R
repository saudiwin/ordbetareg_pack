library(ordbetareg)
library(tidyverse)

# Read the data
# you have data but repo is private
fluency_data <- fluency_data <- read.csv("~/beta_regression_tutorial/data/miko_data.csv") %>%
  # Rename the columns for better readability
  rename(
    "Participant" = "ResponseID", # Rename "ResponseID" to "Participant"
    "Fluency" = "Condition", # Rename "Condition" to "Fluency"
    "Time" = "name", # Rename "name" to "Time"
    "Accuracy" = "value" # Rename "value" to "Accuracy"
  ) %>%
  # Transform the data
  mutate(
    Accuracy = Accuracy * 10 / 100, # Convert Accuracy values to proportions
    Fluency = ifelse(Fluency == 1, "Fluent", "Disflueny"), # rename levels
    Fluency_dummy = ifelse(Fluency == "Fluent", 0, 1), # Recode
    # Fluency: 1 becomes 0, others become 1
    Fluency_dummy = as.factor(Fluency_dummy) # turn fluency cond to dummy code
  ) %>%
  filter(Time == "FreeCt1AVG") %>% # only choose first time point

  # Drop the column "X" and "time" from the dataframe
  select(-X, -Time) %>%
  # move columns around

  relocate(Accuracy, .after = last_col())

m.phi <- ordbetareg(bf(Accuracy ~ Fluency_dummy,
                       phi ~ Fluency_dummy),
                    data=fluency_data,
                    backend = "cmdstanr",
                    file = "model_beta_ordbeta_phi",
                    iter = 200,
                    cores=2,chains=2,
                    phi_reg='both') # log phi

plot_heiss(m.phi, grouping_fac="Fluency_dummy", recode_group_labels = c("Fluency","Disfluency"))

