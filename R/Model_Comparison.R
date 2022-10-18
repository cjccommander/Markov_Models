#### Model Comparison ####

# Model 1: Null Model
# Model 2: Region Model
# Model 3: Otter Occupation Time (Categories) Model
# Model 4: Region + OOT_Cat
# Model 5: Region + OOT + Region:OOT_Cat # run this and compare (Bayes Factor)

# Model 1: null model
X = matrix(1, nrow = n_obs, ncol = 1) # Null model (just a single column of 1's)
stan_data_model1 <- list(
  # number of observations (single integer)
  n_obs = n_obs, 
  # number of states (single integer)
  n_states = n_states,  
  # difference in time among points (vector n_obs long)
  time_diff = time_diff,
  # the design matrix
  X = X,
  # number of parameters to generate each q matrix entry 
  n_pars = ncol(X), # matches the design matrix
  # if we want a Dirichlet-multinomial, set to 1 (we don't)
  dm = 0, 
  # the data (an integer matrix of n_states columns and n_obs rows)
  y_obs = obs_emp
)

# Model 2: region model
X = model.matrix(~ region - 1, data = community_data) # Region model
stan_data_model2 <- list(
  # number of observations (single integer)
  n_obs = n_obs, 
  # number of states (single integer)
  n_states = n_states,  
  # difference in time among points (vector n_obs long)
  time_diff = time_diff,
  # the design matrix
  X = X,
  # number of parameters to generate each q matrix entry 
  n_pars = ncol(X), # matches the design matrix
  # if we want a Dirichlet-multinomial, set to 1 (we don't)
  dm = 0, 
  # the data (an integer matrix of n_states columns and n_obs rows)
  y_obs = obs_emp
)

# Model 3: otter occupation time categories
X = model.matrix(~ OOT_Cat - 1, data = community_data) # OOT categories model
stan_data_model3 <- list(
  # number of observations (single integer)
  n_obs = n_obs, 
  # number of states (single integer)
  n_states = n_states,  
  # difference in time among points (vector n_obs long)
  time_diff = time_diff,
  # the design matrix
  X = X,
  # number of parameters to generate each q matrix entry 
  n_pars = ncol(X), # matches the design matrix
  # if we want a Dirichlet-multinomial, set to 1 (we don't)
  dm = 0, 
  # the data (an integer matrix of n_states columns and n_obs rows)
  y_obs = obs_emp
)

# Model 4: full model (region + otter occupation time)
X = model.matrix(~ region + OOT_Cat - 1, data = community_data) # full model
stan_data_model4 <- list(
  # number of observations (single integer)
  n_obs = n_obs, 
  # number of states (single integer)
  n_states = n_states,  
  # difference in time among points (vector n_obs long)
  time_diff = time_diff,
  # the design matrix
  X = X,
  # number of parameters to generate each q matrix entry 
  n_pars = ncol(X), # matches the design matrix
  # if we want a Dirichlet-multinomial, set to 1 (we don't)
  dm = 0, 
  # the data (an integer matrix of n_states columns and n_obs rows)
  y_obs = obs_emp
)

library(bridgesampling)
# we need to use RStan here (instead of CmdStanR)
stan_model_Rstan = stan_model("~/Markov/Markov_intensities_alt.stan")
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
bayes_factor(bridge_model3, bridge_model4)

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
