#### Bayesian Markov Model v4 ####

# *** Future Work: Develop R function to project comm. state changes under different scenarios 
# (e.g., varying otter occupation time and varying otter density [need Tim for this])

## Have a plot of just self-transitions (4 panel plot (OOTCat))
## Maybe make plot with inverse of from this state/to this state (To this state as panel, from this state as x axis)

# urchin model comm. states: barren, the 4 kelp states, a sparse kelp state (???)

### *** TO DO LIST *** ###
# (1) Model with urchins as independent drivers [for this model, re-define community states 
# without urchins ("barren" instead of urchin state). Start with urchin categories, use criteria from OG community states]
#     (1a) re-define community states: ??? []
#     (1b) create column of urchin density categories: ??? [previously it was >2 or <2]
# (2) Length of transient period: "The half-life to equilibrium evaluates the 
# time of convergence to the steady-state and the length of the transient period.
# the convergence RATE to the equilibrium distribution can be measured 
# using the damping ratio (Hill et al., 2004)" [Brice et al. 2020]: 
# rho = lambda_A1 / lambda_A2 where lambda_A1 and lambda_A2 are the largest 
# and second largest eigenvalues of A (transition prob. matrix)
# The convergence time (length of the transient period) can be approximated 
# using the half-life to equilibrium: t_0.5 = log(2)/log(rho)
# (3) Entropy: use Brice et al. 2020 formula

# Maybe one more model?
# Region + OOT_Cat + Physical Variable(s) (e.g., exposure) ???

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
# select desired columns, including community state columns
# * Select comm. state columns depending on density threshold:
# "comm_state" and "state" = Threshold 1: 2 urchins, 0.5 algae
# "comm_state_Threshold2" and "state_Threshold2" = 2.5 urchins, 1 algae
# "comm_state_Threshold3" and "state_Threshold3" = 1 urchin, 1 algae
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
# * replace "comm_state" with alternative (e.g., "comm_state_Threshold2") as needed

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
# * replace "state" with alternative (e.g., "state_Threshold2") as needed
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

# calculate mean otter occupation time and std dev
OOT_mean = mean(community_data$otter.occupation.time)
OOT_StdDev = sd(community_data$otter.occupation.time)

# scale otter occupation time
community_data = community_data %>% 
  mutate(OOT_scaled = (otter.occupation.time - OOT_mean)/(2*OOT_StdDev))

# Calculate otter occupation time as categories (e.g., none, low, med, high)
# use prefix to re-level factor labels
community_data = community_data %>% 
  mutate(OOT_Cat = ifelse(otter.occupation.time == 0, "A_None", 
                          ifelse(otter.occupation.time > 0 & otter.occupation.time <= 2, "B_Low", 
                                 ifelse(otter.occupation.time > 2 & otter.occupation.time <= 8, "C_Med", 
                                        ifelse(otter.occupation.time > 8, "D_High", "")))))

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
stan_data <- list(
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

# Compile model
stan_model <- cmdstan_model("~/Markov/Markov_intensities_v2.stan",
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

# Save posterior draws
#saveRDS(posterior_draws, "posterior_draws_RegionOOT.rds")

# Extract posterior draws
stanfit <- rstan::read_stan_csv(posterior_draws$output_files())
betas <- rstan::extract(stanfit)$beta

# Save betas
#saveRDS(betas, "betas_RegionOOT.rds")

## Model Outputs ##

# Extract Qmats (# estimated transition intensities)
qmats <- rstan::extract(stanfit)$q_mat_gen
#saveRDS(qmats, "Qmats_RegionOOT.rds")

# Load estimated Qmats
Qmats_RegionOOT = readRDS("~/Markov/Qmats_RegionOOT.rds")
qmats = Qmats_RegionOOT
dim(qmats) # 2500 iterations of intensities/q's for each region--OOT contrast (9) and each state transition (6 x 6)
qmat_df <- melt(qmats) # make into df
names(qmat_df) <- c("iter","contrast","from","to","Qval")
# First, get mean, median, and upper & lower 90% and 95% quantiles
Qmat = qmat_df %>%
  group_by(from, to, contrast) %>%
  summarise(mean = mean(Qval),
            median = median(Qval),
            pct90_lower = quantile(Qval, probs = 0.05),
            pct90_upper = quantile(Qval, probs = 0.95),
            pct95_lower = quantile(Qval, probs = 0.025),
            pct95_upper = quantile(Qval, probs = 0.975))

# filter data and create df's for each contrast
# * Note: Contrasts are ordered per the design matrix above
Qmat_CCBC_None = Qmat %>% filter(contrast == "1") 
Qmat_CCBC_Low = Qmat %>% filter(contrast == "2")
Qmat_CCBC_Med = Qmat %>% filter(contrast == "3")
Qmat_CCBC_High = Qmat %>% filter(contrast == "4")
Qmat_HG_None = Qmat %>% filter(contrast == "5")
Qmat_WCVI_None = Qmat %>% filter(contrast == "6") 
Qmat_WCVI_Low = Qmat %>% filter(contrast == "7")
Qmat_WCVI_Med = Qmat %>% filter(contrast == "8")
Qmat_WCVI_High = Qmat %>% filter(contrast == "9")

# Save Qmats for each contrast
#saveRDS(Qmat_CCBC_None, "Qmat_CCBC_None.rds")
#saveRDS(Qmat_CCBC_Low, "Qmat_CCBC_Low.rds")
#saveRDS(Qmat_CCBC_Med, "Qmat_CCBC_Med.rds")
#saveRDS(Qmat_CCBC_High, "Qmat_CCBC_High.rds")
#saveRDS(Qmat_HG_None, "Qmat_HG_None.rds")
#saveRDS(Qmat_WCVI_None, "Qmat_WCVI_None.rds")
#saveRDS(Qmat_WCVI_Low, "Qmat_WCVI_Low.rds")
#saveRDS(Qmat_WCVI_Med, "Qmat_WCVI_Med.rds")
#saveRDS(Qmat_WCVI_High, "Qmat_WCVI_High.rds")

# Create mean Q matrices for each contrast
# use pivot_wider() to transform df to matrix
# CCBC None
Qmat_CCBC_None = Qmat_CCBC_None %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# CCBC Low
Qmat_CCBC_Low = Qmat_CCBC_Low %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# CCBC Med
Qmat_CCBC_Med = Qmat_CCBC_Med %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# CCBC High
Qmat_CCBC_High = Qmat_CCBC_High %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)

# HG None
Qmat_HG_None = Qmat_HG_None %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)

# WCVI None
Qmat_WCVI_None = Qmat_WCVI_None %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# WCVI Low
Qmat_WCVI_Low = Qmat_WCVI_Low %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# WCVI Med
Qmat_WCVI_Med = Qmat_WCVI_Med %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# WCVI High
Qmat_WCVI_High = Qmat_WCVI_High %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)

## Turnover (sojourn) time: time spent in one state before transitioning to a different state
# Turnover_r = -1 / q_rr where q_rr is the rth entry on the diagonal of the estimated Q matrix
# (i.e., -1 divided by each self transition (diagonal entries))
# Since we have 2500 iterations of q estimates, we need to calculate turnover time 
# for each iteration of self-transitions (q_rr), then average
# Use original qmat_df (all iterations), add turnover time column (-1/Qval)
# summarise, create "turnover time matrices", extract the diagonals (self-transitions)
qmat_df$turnover = -1/qmat_df$Qval

# Get turnover mean, median, and upper & lower 90% and 95% quantiles
stderror = function(x) sd(x)/sqrt(length(x)) # std. error function
Q_turnover = qmat_df %>%
  group_by(from, to, contrast) %>%
  summarise(mean = mean(turnover),
            median = median(turnover),
            sd = sd(turnover),
            stderror = stderror(turnover),
            pct90_lower = quantile(turnover, probs = 0.05),
            pct90_upper = quantile(turnover, probs = 0.95),
            pct95_lower = quantile(turnover, probs = 0.025),
            pct95_upper = quantile(turnover, probs = 0.975))

# filter data and create df's for each contrast
Turnover_CCBC_None = Q_turnover %>% filter(contrast == "1") 
Turnover_CCBC_Low = Q_turnover %>% filter(contrast == "2")
Turnover_CCBC_Med = Q_turnover %>% filter(contrast == "3")
Turnover_CCBC_High = Q_turnover %>% filter(contrast == "4")
Turnover_HG_None = Q_turnover %>% filter(contrast == "5")
Turnover_WCVI_None = Q_turnover %>% filter(contrast == "6") 
Turnover_WCVI_Low = Q_turnover %>% filter(contrast == "7")
Turnover_WCVI_Med = Q_turnover %>% filter(contrast == "8")
Turnover_WCVI_High = Q_turnover %>% filter(contrast == "9")

# select just self-transitions --> rows = c(1,8,15,22,29,36)
# add region_OOT and comm_state
Turnover_CCBC_None = Turnover_CCBC_None[c(1,8,15,22,29,36), 4:11]
Turnover_CCBC_None$region_OOT = "CCBC_None"
Turnover_CCBC_None$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")

Turnover_CCBC_Low = Turnover_CCBC_Low[c(1,8,15,22,29,36), 4:11]
Turnover_CCBC_Low$region_OOT = "CCBC_Low"
Turnover_CCBC_Low$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")

Turnover_CCBC_Med = Turnover_CCBC_Med[c(1,8,15,22,29,36), 4:11]
Turnover_CCBC_Med$region_OOT = "CCBC_Med"
Turnover_CCBC_Med$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")

Turnover_CCBC_High = Turnover_CCBC_High[c(1,8,15,22,29,36), 4:11]
Turnover_CCBC_High$region_OOT = "CCBC_High"
Turnover_CCBC_High$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")

Turnover_HG_None = Turnover_HG_None[c(1,8,15,22,29,36), 4:11]
Turnover_HG_None$region_OOT = "HG_None"
Turnover_HG_None$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")

Turnover_WCVI_None = Turnover_WCVI_None[c(1,8,15,22,29,36), 4:11]
Turnover_WCVI_None$region_OOT = "WCVI_None"
Turnover_WCVI_None$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")

Turnover_WCVI_Low = Turnover_WCVI_Low[c(1,8,15,22,29,36), 4:11]
Turnover_WCVI_Low$region_OOT = "WCVI_Low"
Turnover_WCVI_Low$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")

Turnover_WCVI_Med = Turnover_WCVI_Med[c(1,8,15,22,29,36), 4:11]
Turnover_WCVI_Med$region_OOT = "WCVI_Med"
Turnover_WCVI_Med$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")

Turnover_WCVI_High = Turnover_WCVI_High[c(1,8,15,22,29,36), 4:11]
Turnover_WCVI_High$region_OOT = "WCVI_High"
Turnover_WCVI_High$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")

# combine dfs
Turnover_all = rbind(Turnover_CCBC_None, Turnover_CCBC_Low, Turnover_CCBC_Med, Turnover_CCBC_High, 
                     Turnover_HG_None, 
                     Turnover_WCVI_None, Turnover_WCVI_Low, Turnover_WCVI_Med, Turnover_WCVI_High)

Turnover_all = Turnover_all %>% relocate(region_OOT, .before = mean) %>% 
  relocate(comm_state, .after = region_OOT) %>% rename(mean_turnover = mean) %>% 
  rename(median_turnover = median)

# Plot mean turnover time by contrast (x = comm. state, y = mean turnover time, color = region_OOT)
# ggplot options for visual appeal:
gg_options <- function() theme_bw()+theme(
  panel.grid=element_blank(), # removes ugly grid lines
  plot.background= element_blank(), # removes ugly plot background
  strip.background=element_blank(), # removes ugly strip background
  panel.background= element_blank(), # removes ugly panel background
  legend.title=element_blank(), # not a fan of legend titles if not needed
  legend.background= element_blank(), # removes ugly panel background
  legend.box.background=element_blank(), # removes legend box background
  legend.key= element_blank()) #
# Plot
# re-order comm_state for x-axis and region_OOT for legend [***run plot before this to see if this is necessary]
Turnover_all$comm_state = factor(Turnover_all$comm_state, levels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence"))
Turnover_all$region_OOT = factor(Turnover_all$region_OOT, levels = c("CCBC_None", "CCBC_Low", "CCBC_Med", 
                                                                   "CCBC_High", "HG_None", "WCVI_None", 
                                                                   "WCVI_Low", "WCVI_Med", "WCVI_High"))
# plot with +/- 1 standard error bars (uncertainty around the mean)
ggplot(Turnover_all, aes(comm_state, mean_turnover, fill = region_OOT)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = mean_turnover-stderror, ymax = mean_turnover+stderror), width = 0.2, position = position_dodge(0.9)) + 
  xlab("Community State")+ 
  ylab("Turnover Time") + 
  gg_options()


## Pmats
# Extract Pmats (# estimated transition probabilities)
pmats <- rstan::extract(stanfit)$p_mat_gen
#saveRDS(pmats, "Pmats_RegionOOT.rds")

# Load Pmats
Pmats_RegionOOT = readRDS("~/Markov/Pmats_RegionOOT.rds")
pmats = Pmats_RegionOOT
dim(pmats)
pmat_df <- melt(pmats)
names(pmat_df) <- c("iter","contrast","from","to","Prob")

# First, get mean, median, sd, stderror, and upper & lower 90% and 95% quantiles
Pmat = pmat_df %>%
  group_by(from, to, contrast) %>%
  summarise(mean = mean(Prob),
            median = median(Prob),
            sd = sd(Prob),
            stderror = stderror(Prob),
            pct90_lower = quantile(Prob, probs = 0.05),
            pct90_upper = quantile(Prob, probs = 0.95),
            pct95_lower = quantile(Prob, probs = 0.025),
            pct95_upper = quantile(Prob, probs = 0.975))

# filter data and create df's for each contrast
# * Note: Contrasts are ordered per the design matrix above
Pmat_CCBC_None = Pmat %>% filter(contrast == "1") 
Pmat_CCBC_Low = Pmat %>% filter(contrast == "2")
Pmat_CCBC_Med = Pmat %>% filter(contrast == "3")
Pmat_CCBC_High = Pmat %>% filter(contrast == "4")
Pmat_HG_None = Pmat %>% filter(contrast == "5")
Pmat_WCVI_None = Pmat %>% filter(contrast == "6") 
Pmat_WCVI_Low = Pmat %>% filter(contrast == "7")
Pmat_WCVI_Med = Pmat %>% filter(contrast == "8")
Pmat_WCVI_High = Pmat %>% filter(contrast == "9")

# Save Pmats for each contrast
#saveRDS(Pmat_CCBC_None, "Pmat_CCBC_None.rds")
#saveRDS(Pmat_CCBC_Low, "Pmat_CCBC_Low.rds")
#saveRDS(Pmat_CCBC_Med, "Pmat_CCBC_Med.rds")
#saveRDS(Pmat_CCBC_High, "Pmat_CCBC_High.rds")
#saveRDS(Pmat_HG_None, "Pmat_HG_None.rds")
#saveRDS(Pmat_WCVI_None, "Pmat_WCVI_None.rds")
#saveRDS(Pmat_WCVI_Low, "Pmat_WCVI_Low.rds")
#saveRDS(Pmat_WCVI_Med, "Pmat_WCVI_Med.rds")
#saveRDS(Pmat_WCVI_High, "Pmat_WCVI_High.rds")

# Create mean P matrices for each contrast
# use pivot_wider() to transform df to matrix
# CCBC None
Pmat_CCBC_None = Pmat_CCBC_None %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# CCBC Low
Pmat_CCBC_Low = Pmat_CCBC_Low %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# CCBC Med
Pmat_CCBC_Med = Pmat_CCBC_Med %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# CCBC High
Pmat_CCBC_High = Pmat_CCBC_High %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)

# HG None
Pmat_HG_None = Pmat_HG_None %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)

# WCVI None
Pmat_WCVI_None = Pmat_WCVI_None %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# WCVI Low
Pmat_WCVI_Low = Pmat_WCVI_Low %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# WCVI Med
Pmat_WCVI_Med = Pmat_WCVI_Med %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)
# WCVI High
Pmat_WCVI_High = Pmat_WCVI_High %>% 
  select(from, to, mean) %>% 
  pivot_wider(names_from = to, values_from = mean)


# Steady-State Distribution
# Take the left eigenvector of the transition probability matrix, 
# re-scale so the elements sum to 1 (Norris, 1997). Left eigenvectors = right eigenvectors 
# of the transposed matrix. Use eigen() with transposed matrix (i.e., Pmat dataframe); 
# left eigenvector is the 1st column.

# For each contrast and each iteration, create a P matrix, calculate left eigenvector, store
eigen_array = array(0, dim = c(2500, 9, 6))
mat = matrix(data = 0.0, nrow = 6, ncol = 6)
for(i in 1:2500){
  for(j in 1:9) {
    mat = subset(pmat_df, iter == i & contrast == j)
    mat = t(matrix(mat$Prob, ncol = 6, nrow = 6)) # Left eigenvectors = right eigenvectors 
    # of the transposed matrix. Use eigen() with transposed matrix; 
    # left eigenvector is the 1st column.
    eigen_array[i,j,] = eigen(mat)$vectors[,1]
  }
}

# Average eigenvectors across iterations and by contrasts
eigen_means = as.data.frame(apply(eigen_array, c(2,3), mean)) # ***There should be 9 eigenvectors
eigen_means = eigen_means %>%
  rename(Urchins = V1,
         Understory = V2,
         Epibenthic = V3,
         Macro = V4,
         Mixed = V5,
         Coexistence = V6)
# check to see if apply() is working the way we want
#contrast1 = as.data.frame(eigen_array[,1,])
#colMeans(contrast1) # compare to eigen_means
# Now re-scale each eigenvector so they sum to 1 [use absolute values to change negative signs]
eigen_means = as.data.frame(t(apply(eigen_means[1:6], 1, function(x) abs(x)/sum(abs(x)))))
eigen_means$Region_OOT = c("CCBC_None", "CCBC_Low", "CCBC_Med", "CCBC_High", 
                           "HG_None", 
                           "WCVI_None", "WCVI_Low", "WCVI_Med", "WCVI_High")
eigen_means = eigen_means %>% relocate(Region_OOT, .before = Urchins)

# Save
#write.csv(eigen_means, "SteadyStateDistribution.csv")

# use pivot to get df in order for plotting
eigen_means = eigen_means %>% pivot_longer(names_to = "State", values_to = "Steady", Urchins:Coexistence)
eigen_means
# Modify Region_OOTCat factor levels to get the order we want (not alphabetical order)
# [***run plot before this to see if this is necessary]
eigen_means$Region_OOT = factor(eigen_means$Region_OOT, 
                          levels = c("CCBC_None", "CCBC_Low", "CCBC_Med", "CCBC_High", 
                                     "HG_None", 
                                     "WCVI_None", "WCVI_Low", "WCVI_Med", "WCVI_High"))
# Modify community state factor levels to get the order we want (not alphabetical order)
eigen_means$State = factor(eigen_means$State, 
                         levels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence"))

ggplot(eigen_means, aes(State, Steady, fill = Region_OOT)) + 
  geom_col(position = "dodge") + 
  xlab("Community State")+ 
  ylab("Steady-State Distribution (%)") + 
  gg_options()

### Pmat Plotting ###

# *** Plot of just self-transitions (4 panel plot, 1 for each OOTCat, lines denote regions)
# Make 4 figures, use ggarrange to make 4 panel plot
# Panel 1 (no otters) = contrasts 1, 5, 6
# Panel 2 (low otters) = contrasts 2, 7
# Panel 3 (med otters) = contrasts 3, 8
# Panel 4 (high otters) = contrasts 4, 9

# filter for just self transitions
Pmat_selfs = Pmat # copy Pmat df
Pmat_selfs = Pmat_selfs %>% 
  filter(from == "1" & to == "1" | from == "2" & to == "2" | 
           from == "3" & to == "3" | from == "4" & to == "4" | 
           from == "5" & to == "5" | from == "6" & to == "6")

# add region
Pmat_selfs = Pmat_selfs %>% 
  mutate(region = ifelse(contrast == 1 | contrast == 2 | contrast == 3 | contrast == 4, "CCBC", 
                                 ifelse(contrast == 5, "HG", 
                                        ifelse(contrast == 6 | contrast == 7 | contrast == 8 | contrast == 9, "WCVI", ""))))

# add comm_state
Pmat_selfs = Pmat_selfs %>% 
  mutate(comm_state = ifelse(from == 1, "Urchins", 
                         ifelse(from == 2, "Understory",
                                ifelse(from == 3, "Epibenthic",
                                       ifelse(from == 4, "Macro",
                                              ifelse(from == 5, "Mixed",
                                                     ifelse(from == 6, "Coexistence", "")))))))

# re-order comm_state for x-axis and region_OOT for legend [***run plot before this to see if this is necessary]
Pmat_selfs$comm_state = factor(Pmat_selfs$comm_state, levels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence"))
Pmat_selfs$region = factor(Pmat_selfs$region, levels = c("CCBC", "HG", "WCVI"))

# Panel 1: filter Pmat_selfs df for contrasts with no otters (contrasts 1, 5, 6)
Pmat_OOT_None = Pmat_selfs %>% filter(contrast == "1" | contrast == "5" | contrast == "6")
Panel1 = ggplot(Pmat_OOT_None, aes(comm_state, mean, color = region)) + 
  geom_point(position = position_dodge(0.9)) + 
  scale_color_manual("Region", values = c("CCBC" = "cyan", "HG" = "purple", "WCVI" = "orange")) +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), size = 1, width = 0, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.9)) +
  xlab("Community State")+ 
  ylab("Predicted Self-Transition Probability") + 
  ggtitle("No Otters") +
  gg_options()

# Panel 2: filter Pmat_selfs df for contrasts with low otters (contrasts 2, 7)
Pmat_OOT_Low = Pmat_selfs %>% filter(contrast == "2" | contrast == "7")
Panel2 = ggplot(Pmat_OOT_Low, aes(comm_state, mean, color = region)) + 
  geom_point(position = position_dodge(0.5)) + 
  scale_color_manual("Region", values = c("CCBC" = "cyan", "WCVI" = "orange")) +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), size = 1, width = 0, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.5)) +
  xlab("Community State")+ 
  ylab("Predicted Self-Transition Probability") + 
  ggtitle("Low Otters") +
  gg_options()

# Panel 3: filter Pmat_selfs df for contrasts with med otters (contrasts 3, 8)
Pmat_OOT_Med = Pmat_selfs %>% filter(contrast == "3" | contrast == "8")
Panel3 = ggplot(Pmat_OOT_Med, aes(comm_state, mean, color = region)) + 
  geom_point(position = position_dodge(0.5)) + 
  scale_color_manual("Region", values = c("CCBC" = "cyan", "WCVI" = "orange")) +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), size = 1, width = 0, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.5)) +
  xlab("Community State")+ 
  ylab("Predicted Self-Transition Probability") + 
  ggtitle("Med Otters") +
  gg_options()

# Panel 4: filter Pmat_selfs df for contrasts with high otters (contrasts 4, 9)
Pmat_OOT_High = Pmat_selfs %>% filter(contrast == "4" | contrast == "9")
Panel4 = ggplot(Pmat_OOT_High, aes(comm_state, mean, color = region)) + 
  geom_point(position = position_dodge(0.5)) + 
  scale_color_manual("Region", values = c("CCBC" = "cyan", "WCVI" = "orange")) +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), size = 1, width = 0, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.5)) +
  xlab("Community State")+ 
  ylab("Predicted Self-Transition Probability") + 
  ggtitle("High Otters") +
  gg_options()

library("ggpubr")
self_transitions = ggarrange(Panel1, Panel2, Panel3, Panel4,
                               ncol = 2, nrow = 2)
self_transitions
# save as PDF
#ggsave(".pdf", self_transitions, width=11, height=8.5)

# Just No and High Otters predicted self-transitions
self_transitions2 = ggarrange(Panel1, Panel4,
                             ncol = 1, nrow = 2)
self_transitions2
# save as PDF
#ggsave(".pdf", self_transitions2, width=11, height=8.5)

# *** Plot of state transitions (6 panels, just none and high OOT categories)
# add region to Pmat df
Pmat = Pmat %>% 
  mutate(region = ifelse(contrast == 1 | contrast == 2 | contrast == 3 | contrast == 4, "CCBC", 
                         ifelse(contrast == 5, "HG", 
                                ifelse(contrast == 6 | contrast == 7 | contrast == 8 | contrast == 9, "WCVI", ""))))

# filter for just No Otters contrasts (contrasts 1, 5, 6)
Pmat_None = Pmat %>% filter(contrast == "1" | contrast == "5" | contrast == "6")
StateTransitions_NoOtters = ggplot(aes(y = mean, x = factor(to), color = factor(region)), data = Pmat_None)+
  geom_point(position = position_dodge(0.9)) + 
  scale_color_manual(name = "Region", values = c("CCBC" = "cyan", "HG" = "purple", "WCVI" = "orange")) +
  geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.9)) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("Transition To This State")+
  ylab("Predicted Transition Probability") + 
  ggtitle("No Otters") +
  gg_options()
StateTransitions_NoOtters
# save as PDF
#ggsave("State Transitions_No Otters.pdf", StateTransitions_NoOtters, width=11, height=8.5)

# filter for just High Otters contrasts (contrasts 4, 9)
Pmat_High = Pmat %>% filter(contrast == "4" | contrast == "9")
StateTransitions_HighOtters = ggplot(aes(y = mean, x = factor(to), color = factor(region)), data = Pmat_High)+
  geom_point(position = position_dodge(0.5)) + 
  scale_color_manual(name = "Region", values = c("CCBC" = "cyan", "WCVI" = "orange")) +
  geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.5)) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("Transition To This State")+
  ylab("Predicted Transition Probability") + 
  ggtitle("High Otters") +
  gg_options()
StateTransitions_HighOtters
# save as PDF
#ggsave("State Transitions_High Otters.pdf", StateTransitions_HighOtters, width=11, height=8.5)


########### Old plot code ###########

## Plot Pmat by OOT_Cat (HG_None, WCVI_None, CCBC_None)
# filter Pmat df for OOT_Cat contrasts: OOTCat = "None"
Pmat_OOT_None = pmat_df %>% filter(contrast == "1" | contrast == "5" | contrast == "6")

ggplot(aes(y = Prob, x = factor(to), color = factor(contrast)), data = Pmat_OOT_None)+
  geom_pointrange(
    stat = "summary",
    fun = mean, position = position_dodge(width = 0.5)
  ) + 
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.025)},
    fun.max = function(z) {quantile(z,0.975)},
    size= 1, position = position_dodge(width = 0.5)
  ) + 
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.05)},
    fun.max = function(z) {quantile(z,0.95)},
    size= 2, position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("To This State")+
  ylab("Predicted Transition Probability")+
  scale_colour_discrete(name = "Region_OOT", labels = c("CCBC_None", "HG_None", "WCVI_None")) +
  gg_options()

# filter Pmat df for OOT_Cat contrasts: OOTCat = "Low"
Pmat_OOT_Low = pmat_df %>% filter(contrast == "2" | contrast == "7")

ggplot(aes(y = Prob, x = factor(to), color = factor(contrast)), data = Pmat_OOT_Low)+
  geom_pointrange(
    stat = "summary",
    fun = mean, position = position_dodge(width = 0.5)
  )+
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.025)},
    fun.max = function(z) {quantile(z,0.975)},
    size= 1, position = position_dodge(width = 0.5)
  ) + 
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.05)},
    fun.max = function(z) {quantile(z,0.95)},
    size= 2, position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("To This State")+
  ylab("Predicted Transition Probability")+
  scale_colour_discrete(name = "Region_OOT", labels = c("CCBC_Low", "WCVI_Low")) +
  gg_options()


# filter Pmat df for OOT_Cat contrasts: OOTCat = "Med"
Pmat_OOT_Med = pmat_df %>% filter(contrast == "3" | contrast == "8")

ggplot(aes(y = Prob, x = factor(to), color = factor(contrast)), data = Pmat_OOT_Med)+
  geom_pointrange(
    stat = "summary",
    fun = mean, position = position_dodge(width = 0.5)
  )+
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.025)},
    fun.max = function(z) {quantile(z,0.975)},
    size= 1, position = position_dodge(width = 0.5)
  ) + 
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.05)},
    fun.max = function(z) {quantile(z,0.95)},
    size= 2, position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("To This State")+
  ylab("Predicted Transition Probability")+
  scale_colour_discrete(name = "Region_OOT", labels = c("CCBC_Med", "WCVI_Med")) +
  gg_options()

# filter Pmat df for OOT_Cat contrasts: OOTCat = "High"
Pmat_OOT_High = pmat_df %>% filter(contrast == "4" | contrast == "9")

ggplot(aes(y = Prob, x = factor(to), color = factor(contrast)), data = Pmat_OOT_High)+
  geom_pointrange(
    stat = "summary",
    fun = mean, position = position_dodge(width = 0.5)
  )+
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.025)},
    fun.max = function(z) {quantile(z,0.975)},
    size= 1, position = position_dodge(width = 0.5)
  ) + 
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.05)},
    fun.max = function(z) {quantile(z,0.95)},
    size= 2, position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("To This State")+
  ylab("Predicted Transition Probability")+
  scale_colour_discrete(name = "Region_OOT", labels = c("CCBC_High", "WCVI_High")) +
  gg_options()


## Plot Pmat by region (we can't fit all 9 contrasts on 1 figure)
# filter Pmat df for just WCVI contrasts (contrast = 6-9)
Pmat_WCVI = pmat_df %>% filter(contrast == "6" | contrast == "7" | contrast == "8" | contrast == "9")

ggplot(aes(y = Prob, x = factor(to), color = factor(contrast)), data = Pmat_WCVI)+
  geom_pointrange(
    stat = "summary",
    fun = mean, position = position_dodge(width = 0.5)
  )+
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.025)},
    fun.max = function(z) {quantile(z,0.975)},
    size= 1, position = position_dodge(width = 0.5)
  ) + 
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.05)},
    fun.max = function(z) {quantile(z,0.95)},
    size= 2, position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("To This State")+
  ylab("Predicted Transition Probability")+
  scale_colour_discrete(name = "Region_OOT", labels = c("WCVI_None", "WCVI_Low", "WCVI_Med", "WCVI_High")) +
  gg_options()

# filter Pmat df for just CCBC contrasts (contrast = 1-4)
Pmat_CCBC = pmat_df %>% filter(contrast == "1" | contrast == "2" | contrast == "3" | contrast == "4")

ggplot(aes(y = Prob, x = factor(to), color = factor(contrast)), data = Pmat_CCBC)+
  geom_pointrange(
    stat = "summary",
    fun = mean, position = position_dodge(width = 0.5)
  )+
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.025)},
    fun.max = function(z) {quantile(z,0.975)},
    size= 1, position = position_dodge(width = 0.5)
  ) + 
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.05)},
    fun.max = function(z) {quantile(z,0.95)},
    size= 2, position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("To This State")+
  ylab("Predicted Transition Probability")+
  scale_colour_discrete(name = "Region_OOT", labels = c("CCBC_None", "CCBC_Low", "CCBC_Med", "CCBC_High")) +
  gg_options()

# filter Pmat df for just HG contrast (contrast = 5)
Pmat_HG = pmat_df %>% filter(contrast == "5")

ggplot(aes(y = Prob, x = factor(to), color = factor(contrast)), data = Pmat_HG)+
  geom_pointrange(
    stat = "summary",
    fun = mean, position = position_dodge(width = 0.5)
  )+
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.025)},
    fun.max = function(z) {quantile(z,0.975)},
    size= 1, position = position_dodge(width = 0.5)
  ) + 
  geom_linerange(
    stat = "summary",
    fun.min = function(z) {quantile(z,0.05)},
    fun.max = function(z) {quantile(z,0.95)},
    size= 2, position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("To This State")+
  ylab("Predicted Transition Probability")+
  scale_colour_discrete(name = "Region_OOT", labels = c("HG_None")) +
  gg_options()
