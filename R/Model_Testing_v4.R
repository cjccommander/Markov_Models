#### Bayesian Markov Model v4 ####

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

### NOTES ###

# "estimated baseline hazards (i.e., instantaneous
# risk of moving from one state to another when all covariates are set to 0 
# (i.e., the means of standardized covariates)" - From Brice et al. 2020
# These show the background rate of community changes for the study area 
# (i.e., the greatest transition intensities are the dominant state transitions)

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
community_data = read.csv("~/Markov/CommunityStates_WCVIx2_CCBC_HG_combined.csv") %>% 
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
# in time differential column, replace NA's with 0
community_data[, 10][is.na(community_data[, 10])] = 0

# remove the 3 obs. with the "empty" state category
community_data = community_data %>% filter(comm_state != "empty")

# number of observations to simulate
n_obs = 581
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

# Which contrasts I'm interested in:
# [rows in design matrix X that correspond to region and OOTCat]
contrast_rows = c(207, 1, 4, 20, 18, 5, 31, 8, 9) # 207 = CCBC, none; 1 = CCBC, low; etc.
# 18 = HG, none (no otters)
# 5 = WCVI, none; 31 = WCVI, low; etc.

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

## Qmats [*** loop over Qmats, get mean and CI for each Qmat, then take average (same procedure as steady state)]
# Extract Qmats (# estimated transition intensities/hazards)
qmats <- rstan::extract(stanfit)$q_mat_gen
#saveRDS(qmats, "Qmats_RegionOOT.rds")
# load Qmats
Qmats_RegionOOT = readRDS("~/Markov/Qmats_RegionOOT.rds")
qmats = Qmats_RegionOOT
dim(qmats) # 2500 iterations for each region--OOT contrast (9) and each state transition (6 x 6)
qmat_df <- melt(qmats)
names(qmat_df) <- c("iter","contrast","from","to","Qval")
# Create estimated transition intensity/hazard matrices for each contrast
# First, get mean, median, and upper & lower 90% and 95% intervals
Qmat = qmat_df %>%
  group_by(from, to, contrast) %>%
  summarise(mean = mean(Qval),
            median = median(Qval),
            pct90_lower = quantile(Qval, probs = 0.05),
            pct90_upper = quantile(Qval, probs = 0.95),
            pct95_lower = quantile(Qval, probs = 0.025), # maybe show 2 error bars: 90% and 95%
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


## Turnover/sojourn time:
# Turnover_r = -1 / Qrr where qrr is the rth entry on the diagonal of the estimated Q matrix
# -1 divided by each self transition (diagonal entries)
Sojourn_CCBC_None = -1/Qmat_CCBC_None
Sojourn_CCBC_Low = -1/Qmat_CCBC_Low
Sojourn_CCBC_Med = -1/Qmat_CCBC_Med
Sojourn_CCBC_High = -1/Qmat_CCBC_High
Sojourn_HG_None = -1/Qmat_HG_None
Sojourn_WCVI_None = -1/Qmat_WCVI_None
Sojourn_WCVI_Low = -1/Qmat_WCVI_Low
Sojourn_WCVI_Med = -1/Qmat_WCVI_Med
Sojourn_WCVI_High = -1/Qmat_WCVI_High # did WCVI Med and WCVI High get switched? Plot results and look. 
# WCVI Med results look like what you'd expect from WCVI High, and vice versa.

# extract diagonal (turnover time) from each matrix and create df
Sojourn_CCBC_None = as.data.frame(diag(as.matrix(Sojourn_CCBC_None[, -1])))
Sojourn_CCBC_None$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")
Sojourn_CCBC_None$region_OOT = "CCBC_None"
names(Sojourn_CCBC_None)[1] = "turnover_time"

Sojourn_CCBC_Low = as.data.frame(diag(as.matrix(Sojourn_CCBC_Low[, -1])))
Sojourn_CCBC_Low$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")
Sojourn_CCBC_Low$region_OOT = "CCBC_Low"
names(Sojourn_CCBC_Low)[1] = "turnover_time"

Sojourn_CCBC_Med = as.data.frame(diag(as.matrix(Sojourn_CCBC_Med[, -1])))
Sojourn_CCBC_Med$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")
Sojourn_CCBC_Med$region_OOT = "CCBC_Med"
names(Sojourn_CCBC_Med)[1] = "turnover_time"

Sojourn_CCBC_High = as.data.frame(diag(as.matrix(Sojourn_CCBC_High[, -1])))
Sojourn_CCBC_High$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")
Sojourn_CCBC_High$region_OOT = "CCBC_High"
names(Sojourn_CCBC_High)[1] = "turnover_time"

Sojourn_HG_None = as.data.frame(diag(as.matrix(Sojourn_HG_None[, -1])))
Sojourn_HG_None$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")
Sojourn_HG_None$region_OOT = "HG_None"
names(Sojourn_HG_None)[1] = "turnover_time"

Sojourn_WCVI_None = as.data.frame(diag(as.matrix(Sojourn_WCVI_None[, -1])))
Sojourn_WCVI_None$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")
Sojourn_WCVI_None$region_OOT = "WCVI_None"
names(Sojourn_WCVI_None)[1] = "turnover_time"

Sojourn_WCVI_Low = as.data.frame(diag(as.matrix(Sojourn_WCVI_Low[, -1])))
Sojourn_WCVI_Low$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")
Sojourn_WCVI_Low$region_OOT = "WCVI_Low"
names(Sojourn_WCVI_Low)[1] = "turnover_time"

Sojourn_WCVI_Med = as.data.frame(diag(as.matrix(Sojourn_WCVI_Med[, -1])))
Sojourn_WCVI_Med$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")
Sojourn_WCVI_Med$region_OOT = "WCVI_Med"
names(Sojourn_WCVI_Med)[1] = "turnover_time"

Sojourn_WCVI_High = as.data.frame(diag(as.matrix(Sojourn_WCVI_High[, -1])))
Sojourn_WCVI_High$comm_state = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")
Sojourn_WCVI_High$region_OOT = "WCVI_High"
names(Sojourn_WCVI_High)[1] = "turnover_time"

# combine dfs
Sojourn_all = rbind(Sojourn_CCBC_None, Sojourn_CCBC_Low, Sojourn_CCBC_Med, Sojourn_CCBC_High, 
                    Sojourn_HG_None, 
                    Sojourn_WCVI_None, Sojourn_WCVI_Low, Sojourn_WCVI_Med, Sojourn_WCVI_High)

# Plot turnover time by contrast (x = comm. state, y = turnover time, color = region_OOT)
# need 3 columns: contrast, comm. state, turnover time
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
# re-order comm_state for x-axis and region_OOT for legend
Sojourn_all$comm_state = factor(Sojourn_all$comm_state, levels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence"))
Sojourn_all$region_OOT = factor(Sojourn_all$region_OOT, levels = c("CCBC_None", "CCBC_Low", "CCBC_Med", 
                                                                   "CCBC_High", "HG_None", "WCVI_None", 
                                                                   "WCVI_Low", "WCVI_Med", "WCVI_High"))
ggplot(Sojourn_all, aes(comm_state, turnover_time, fill = region_OOT)) + 
  geom_col(position = "dodge") + 
  xlab("Community State")+ 
  ylab("Turnover Time") + 
  gg_options()


## Pmats
# Extract Pmats
pmats <- rstan::extract(stanfit)$p_mat_gen
#saveRDS(pmats, "Pmats_RegionOOT.rds")
# load Pmats
Pmats_RegionOOT = readRDS("~/Markov/Pmats_RegionOOT.rds")
pmats = Pmats_RegionOOT
dim(pmats)
pmat_df <- melt(pmats)
names(pmat_df) <- c("iter","contrast","from","to","Prob")

## Plot Pmat by OOT_Cat (HG_None, WCVI_None, CCBC_None)
# filter Pmat df for OOT_Cat contrasts: OOTCat = "None"
Pmat_OOT_None = pmat_df %>% filter(contrast == "1" | contrast == "5" | contrast == "6")

ggplot(aes(y = Prob, x = factor(to), color = factor(contrast)), data = Pmat_OOT_None)+
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


## # Steady-state distribution:
# Take the left eigenvector of the transition probability matrix, 
# re-scale so the elements sum to 1 (Norris, 1997). Left eigenvectors = right eigenvectors 
# of the transposed matrix. Use eigen() with transposed matrix (i.e., Pmat dataframe); 
# left eigenvector is the 1st column.

# [*** We'll need to use the raw Pmat file for generating left eigenvectors for each iteration (filtered by category).
# Then average those vectors (within each category) to get the steady-state distribution]
# For each OOTCat, loop over each of the 250 iterations of Pmat, 
# create a matrix, calculate left eigenvector, store, average

