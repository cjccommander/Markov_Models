#### Model Comparison ####
# *** update to reflect new model code

# Model 1: Null Model
# Model 2: Region Model
# Model 3: Otter Occupation Time (Categories) Model
# Model 4: Region + OOT_Cat

library(cmdstanr)
library(posterior)
library(tidyverse)
library(parallel)
library(rstan)
library(reshape2)
library(ggplot2)
library(Matrix)
rstan::rstan_options(javascript=FALSE)

# Empirical Data
community_data = read.csv("~/Markov/Markov_Models/R/CommunityStates_WCVIx2_CCBC_HG_combined_REVISED.csv") %>% 
  arrange(site, year)
# select desired columns
community_data = community_data %>% 
  select(region, location, site, month, year, year.of.otter.occupancy,
         otter.occupation.time, state, comm_state)

# We need to create a column of time differentials (time passed from the prior
# time point, i.e., previous survey year)
community_data = community_data %>% 
  group_by(site) %>% 
  mutate(time_diff = year - lag(year))
# in time differential column, replace NA's with 0
community_data[, 10][is.na(community_data[, 10])] = 0

# remove the obs. with the "empty" state category
community_data = community_data %>% filter(comm_state != "empty")

# number of observations to simulate
n_obs = nrow(community_data)
# number of states to simulate
n_states = 6
# time difference among observations (the initial is zero)
time_diff = community_data$time_diff
# number of chains to run in the Stan model
cores = detectCores()
n_chains = min(10,cores-3)
# number of iterations to run in the Stan model per chain
n_iter = 250

# We need to create a matrix of observed empirical states, each row being a single 
# observation (i.e., site-year) with a 1 or 0 for each state column.
# In the community_data df, create a column for each community state and 
# assign a 0 or 1 for each row's observed state:
community_data = community_data  %>% 
  mutate(state_1 = ifelse(state == 1, 1, 0)) %>% # urchins
  mutate(state_2 = ifelse(state == 2, 1, 0)) %>% # understory
  mutate(state_3 = ifelse(state == 3, 1, 0)) %>% # epibenthic_canopy
  mutate(state_4 = ifelse(state == 4, 1, 0)) %>% # macro_canopy
  mutate(state_5 = ifelse(state == 5, 1, 0)) %>% # mixed_algae
  mutate(state_6 = ifelse(state == 6, 1, 0))     # coexistence

# matrix of empirical data (observed states) for below
obs_emp = as.matrix(community_data[, c(11:16)])

# replace otter occupation time NAs with 0
community_data = community_data %>% 
  mutate(otter.occupation.time = ifelse(is.na(otter.occupation.time), as.numeric(0), otter.occupation.time))

# Calculate otter occupation time as categories (e.g., none, low, med, high)
# use prefix to re-level factor labels
community_data = community_data %>% 
  mutate(OOT_Cat = ifelse(otter.occupation.time == 0, "A_None", 
                          ifelse(otter.occupation.time > 0 & otter.occupation.time <= 2, "B_Low", 
                                 ifelse(otter.occupation.time > 2 & otter.occupation.time <= 8, "C_Med", 
                                        ifelse(otter.occupation.time > 8, "D_High", "")))))

## Model 1: Null Model ##

# Design matrix
X = model.matrix(~ 1, data = community_data) # just the intercept (a single column of 1's)

# create Stan data
stan_data_model1 <- list(
  # number of observations (single integer)
  n_obs = n_obs, 
  # number of states (single integer)
  n_states = n_states,
  #
  n_conts = length(1), # no contrasts since it's just the intercept
  # difference in time among points (vector n_obs long)
  time_diff = time_diff,
  # the design matrix
  X = X,
  # 
  conts = 1,
  # number of parameters to generate each q matrix entry 
  n_pars = ncol(X), # matches the design matrix
  # if we want a Dirichlet-multinomial, set to 1 (we don't)
  dm = 0, 
  # the data (an integer matrix of n_states columns and n_obs rows)
  y_obs = obs_emp
)

## Model 2: Region Model ##

# Design matrix
X = model.matrix(~ region - 1, data = community_data)

# Contrasts we're interested in:
# (i.e., rows in the design matrix (X) that correspond to the region)
contrast_rows = c(1, 18, 5)
# 3 contrasts:
# 1 = CCBC
# 18 = HG
# 5 = WCVI

# create Stan data
stan_data_model2 <- list(
  # number of observations (single integer)
  n_obs = n_obs, 
  # number of states (single integer)
  n_states = n_states,
  #
  n_conts = length(contrast_rows),
  # difference in time among points (vector n_obs long)
  time_diff = time_diff,
  # the design matrix
  X = X,
  # 
  conts = contrast_rows,
  # number of parameters to generate each q matrix entry 
  n_pars = ncol(X), # matches the design matrix
  # if we want a Dirichlet-multinomial, set to 1 (we don't)
  dm = 0, 
  # the data (an integer matrix of n_states columns and n_obs rows)
  y_obs = obs_emp
)

## Model 3: Otter Occupation Time (categories) Model ##

# Design matrix
X = model.matrix(~ OOT_Cat - 1, data = community_data)

# Contrasts we're interested in:
# (i.e., rows in the design matrix (X) that correspond to the otter occupation time category)
contrast_rows = c(5, 1, 4, 9)
# 4 contrasts:
# 5 = None (no otters)
# 1 = Low
# 4 = Med
# 9 = High

# create Stan data
stan_data_model3 <- list(
  # number of observations (single integer)
  n_obs = n_obs, 
  # number of states (single integer)
  n_states = n_states,
  #
  n_conts = length(contrast_rows),
  # difference in time among points (vector n_obs long)
  time_diff = time_diff,
  # the design matrix
  X = X,
  # 
  conts = contrast_rows,
  # number of parameters to generate each q matrix entry 
  n_pars = ncol(X), # matches the design matrix
  # if we want a Dirichlet-multinomial, set to 1 (we don't)
  dm = 0, 
  # the data (an integer matrix of n_states columns and n_obs rows)
  y_obs = obs_emp
)

## Model 4: Full Model (Region + OOT) ##

# Design matrix
X = model.matrix(~ region + OOT_Cat - 1, data = community_data)

# Contrasts we're interested in:
# (i.e., rows in the design matrix (X) that correspond to the region and otter occupation time category)
contrast_rows = c(207, 1, 4, 20, 18, 5, 31, 8, 9)
# 9 contrasts:
# 207 = CCBC_None; 1 = CCBC_Low; 4 = CCBC_Med; 20 = CCBC_High
# 18 = HG_None
# 5 = WCVI_None; 31 = WCVI_Low; 8 = WCVI_Med; 9 = WCVI_High

# create Stan data
stan_data_model4 <- list(
  # number of observations (single integer)
  n_obs = n_obs, 
  # number of states (single integer)
  n_states = n_states,
  #
  n_conts = length(contrast_rows),
  # difference in time among points (vector n_obs long)
  time_diff = time_diff,
  # the design matrix
  X = X,
  # 
  conts = contrast_rows,
  # number of parameters to generate each q matrix entry 
  n_pars = ncol(X), # matches the design matrix
  # if we want a Dirichlet-multinomial, set to 1 (we don't)
  dm = 0, 
  # the data (an integer matrix of n_states columns and n_obs rows)
  y_obs = obs_emp
)

library(bridgesampling)
# we need to use RStan here (instead of CmdStanR)
stan_model_Rstan = stan_model("~/Markov/Markov_Models/Stan/Markov_intensities_v2_model_comparison.stan")
posterior_model1 <- rstan::sampling(stan_model_Rstan, data = stan_data_model1, iter = 20000, warmup = 1000, chains = 4, cores = 4)
posterior_model2 <- rstan::sampling(stan_model_Rstan, data = stan_data_model2, iter = 20000, warmup = 1000, chains = 4, cores = 4)
posterior_model3 <- rstan::sampling(stan_model_Rstan, data = stan_data_model3, iter = 20000, warmup = 1000, chains = 4, cores = 4)
posterior_model4 <- rstan::sampling(stan_model_Rstan, data = stan_data_model4, iter = 20000, warmup = 1000, chains = 4, cores = 4)

# save models
saveRDS(posterior_model1, "posterior_model1_FINAL.rds")
saveRDS(posterior_model2, "posterior_model2_FINAL.rds")
saveRDS(posterior_model3, "posterior_model3_FINAL.rds")
saveRDS(posterior_model4, "posterior_model4_FINAL.rds")

bridge_model1 <- bridge_sampler(posterior_model1)
bridge_model2 <- bridge_sampler(posterior_model2)
bridge_model3 <- bridge_sampler(posterior_model3)
bridge_model4 <- bridge_sampler(posterior_model4)

# save bridge models
saveRDS(bridge_model1, "bridge_model1.rds")
saveRDS(bridge_model2, "bridge_model2.rds")
saveRDS(bridge_model3, "bridge_model3.rds")
saveRDS(bridge_model4, "bridge_model4.rds")

bridge_model1 # bridge_model1$logml # log marginal likelihood
bridge_model2
bridge_model3
bridge_model4

# Bayes Factor: examines relative evidence of one model over another 
# (i.e., the relative probability of the observed data under each model)
bayes_factor(bridge_model1, bridge_model2)

# posterior probability of support
post_prob(bridge_model1, bridge_model2, bridge_model3, bridge_model4)

# compute error measures for estimated marginal likelihood [*** ask Dan about error measures]
error_measures(bridge_model1)
error_measures(bridge_model2)
error_measures(bridge_model3)
error_measures(bridge_model4)
# re2: approximate relative mean squared error for marginal likelihood estimate
# cv: coefficient of variation for marginal likelihood estimate (assumes that bridge estimate is unbiased)
# percentage: percentage error of marginal likelihood estimate
