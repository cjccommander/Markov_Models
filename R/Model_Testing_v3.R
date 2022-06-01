#### Bayesian Markov Model v3 ####

library(cmdstanr)
library(posterior)
library(dplyr)
library(parallel)
library(rstan)
library(reshape2)
library(ggplot2)
library(Matrix)
rstan::rstan_options(javascript=FALSE)

# Empirical Data
community_data = read.csv("CommunityStates_WCVIx2_CCBC_HG_combined.csv") %>% 
  arrange(site, year)
# select desired columns
community_data = community_data %>% 
  select(region, location, site, month, year, year.of.otter.occupancy,
         otter.occupation.time, state, comm_state)
# We need to create a column of time differentials (time passed from the prior
# time point, i.e. previous survey year)
community_data = community_data %>% 
  group_by(site) %>% 
  mutate(time_diff = year - lag(year))
# replace NAs with 0
community_data[, 10][is.na(community_data[, 10])] = 0

# number of observations to simulate
n_obs = 584
# number of states to simulate
n_states = 6
# time difference among observations (the initial is zero)
time_diff = community_data$time_diff
# number of chains to run in the Stan model
cores = detectCores()
n_chains = min(10,cores-3)
# number of iterations to run in the Stan model per chain
n_iter = 250

#### *** What's need from here to below??? *** ####

# *** [Replace this with true observed intensities?] ***
# simulate "true" intensities
set.seed(123)
intensities <- exp(rnorm(n_states * (n_states - 1),0,1) * 1 - 0.5*log(n_states-1))
intensities <- matrix(intensities, nrow= n_states)

# build Q matrix
q_mat = matrix(0, ncol = n_states, nrow = n_states)
for (i in 1 : n_states){
  for (j in 1 : n_states){
    # if parameter is before the diagonal, place in order
    if (j < i){
      q_mat[i, j] = intensities[i, j];
    } 
    # if parameter is after the diagonal, place one ahead
    if (j > i){
      q_mat[i, j] = intensities[i, j - 1];
    }
    # if diagonal, temporarily hold with zero
    if (j == i){
      q_mat[i, j] = 0;
    }
  } 
  #set diagonal to sum of the off diagonals for the row
  q_mat[i, i] = - sum(q_mat[i, 1 : n_states])
}
p_mat_true = as.matrix(expm(q_mat))

# unit test: ensure transition probabilities after unit time sum to one
all.equal(rowSums(diag(rep(1, n_states)) %*% expm(q_mat)), rep(1, n_states))
true_df <- melt(as.matrix(expm(q_mat)))
names(true_df) = c("from","to","prob")

# simulate one realization for the time domain defined (n_obs)
set.seed(123)
out <- matrix(0, ncol = n_states, nrow = n_obs)
out[1, ] <- rmultinom(1, size = 1, prob= rep(1 / n_states, n_states)) 
for (i in 2 : n_obs){
  out[i, ] = rmultinom(1, 1, as.vector(out[i - 1, ] %*% 
                                         expm(time_diff[i] * q_mat)))
}

#### *** What's need from here to above??? *** ####

# We need to create a matrix of observed empirical states (each row being a single 
# observation, i.e. site-year, with a 1 or 0 for each state column.
# In the community_data df, create a column for each community state and 
# assign a 0 or 1 for each row's observed state:
community_data = community_data  %>% 
  mutate(state_1 = ifelse(state == 1, 1, 0)) %>% 
  mutate(state_2 = ifelse(state == 2, 1, 0)) %>% 
  mutate(state_3 = ifelse(state == 3, 1, 0)) %>% 
  mutate(state_4 = ifelse(state == 4, 1, 0)) %>% 
  mutate(state_5 = ifelse(state == 5, 1, 0)) %>% 
  mutate(state_6 = ifelse(state == 6, 1, 0))

# matrix of empirical data (observed states) for below
obs_emp = as.matrix(community_data[, c(11:16)])

# Simulate parameters
indxs = matrix(0,nrow = n_states, ncol= n_states-1)
for(i in 1:n_states){
  tmp = seq(1,n_states); tmp = tmp[-i]
  indxs[i,] = tmp
}
rm(tmp)
# create Stan data
stan_data <- list(
  # number of observations (single integer)
  n_obs = n_obs, 
  # number of states (single integer)
  n_states = n_states,  
  # difference in time among points (vector n_obs long)
  time_diff = time_diff,
  # number of parameters to generate each q matrix entry 
  n_pars = 1, 
  # if we want a Dirichlet-multinomial, set to 1 (we don't)
  dm = 0, 
  # the data (an integer matrix of n_states columns and n_obs rows)
  y_obs = obs_emp,
  # the design matrix (her just a single column of 1s)
  X = matrix(1, nrow = n_obs, ncol= 1),
  # index of column numbers per row excluding diagonal
  indxs = indxs
)

# Compile model
stan_model <- cmdstan_model("../Stan/Markov_intensities_v2.stan",
                            cpp_options = list(opencl = TRUE))

stan_model

# Run the model
system.time(
  suppressMessages(                # Suppress messages/warnings (if desired)
    suppressWarnings (
      posterior_draws <- stan_model$sample(
        data = stan_data,
        iter_warmup = n_iter,
        iter_sampling = n_iter,
        chains = n_chains,
        refresh = ceiling(n_iter/5),
        parallel_chains = getOption("mc.cores", n_chains)
      )
    )
  )
)
