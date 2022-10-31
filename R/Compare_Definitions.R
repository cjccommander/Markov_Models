#### Testing Alternative Density Thresholds ####

# Test if model results are robust to different community state definitions. 
# Compare community states based on alternative density thresholds.

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
# time point, i.e. previous survey year)
community_data = community_data %>% 
  group_by(site) %>% 
  mutate(time_diff = year - lag(year))
# in time differential column, replace NA's with 0
community_data[, 10][is.na(community_data[, 10])] = 0

# remove obs. with the "empty" state category
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

# Which contrasts I'm interested in:
# [rows in design matrix X that correspond to region and OOTCat]
contrast_rows = c(207, 1, 4, 20, 18, 5, 31, 8, 9) # 207 = CCBC_none; 1 = CCBC_low; etc.
# 18 = HG_none (no otters)
# 5 = WCVI_none; 31 = WCVI_low; etc.

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
stan_model <- cmdstan_model("~/Markov/Markov_Models/Stan/Markov_intensities_v2.stan",
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
#saveRDS(posterior_draws, "posterior_draws_Threshold1.rds")

# Extract posterior draws
stanfit <- rstan::read_stan_csv(posterior_draws$output_files())
# Extract betas
betas <- rstan::extract(stanfit)$beta

# Save betas
#saveRDS(betas, "betas_Threshold1.rds")

# Extract Qmats (# estimated transition intensities)
qmats <- rstan::extract(stanfit)$q_mat_gen
#saveRDS(qmats, "Qmats_Threshold1.rds")

# Extract Pmats (# estimated transition probabilities)
pmats <- rstan::extract(stanfit)$p_mat_gen
#saveRDS(pmats, "Pmats_Threshold1.rds")

## Model Outputs ##

# Load estimated Qmats
qmats = readRDS("~/Markov/Markov_Models/Qmats_Threshold1.rds")
dim(qmats) # 2500 iterations of intensities/q's for each state transition (6 x 6)
qmat_df <- melt(qmats) # make into df
names(qmat_df) <- c("iter","contrast","from","to","Qval")
# First, get mean, median, and upper & lower 90% and 95% intervals
stderror = function(x) sd(x)/sqrt(length(x)) # std. error function
Qmat = qmat_df %>%
  group_by(from, to, contrast) %>%
  summarise(mean = mean(Qval),
            median = median(Qval),
            sd = sd(Qval),
            stderror = stderror(Qval),
            pct90_lower = quantile(Qval, probs = 0.05),
            pct90_upper = quantile(Qval, probs = 0.95),
            pct95_lower = quantile(Qval, probs = 0.025),
            pct95_upper = quantile(Qval, probs = 0.975))

# After reading (line 163) and processing Qmats_Threshold1, Qmats_Threshold2, Qmats_Threshold3:
# *** Create df of q's from different thresholds, filter for just WCVI and no & high otters (contrasts 6 and 9),
# create figure of q values compared across the 3 thresholds (2 6-panel figures with 3 colors for thresholds)
# The q values will differ but the dominant transitions (highest q's) should be consistent
Qs_Threshold1 = Qmat
Qs_Threshold2 = Qmat
Qs_Threshold3 = Qmat

Qs_Threshold1_WCVI_NoneHigh = Qs_Threshold1 %>% filter(contrast == "6" | contrast == "9")
Qs_Threshold2_WCVI_NoneHigh = Qs_Threshold2 %>% filter(contrast == "6" | contrast == "9")
Qs_Threshold3_WCVI_NoneHigh = Qs_Threshold3 %>% filter(contrast == "6" | contrast == "9")

# add threshold designation
Qs_Threshold1_WCVI_NoneHigh$threshold = "1"
Qs_Threshold2_WCVI_NoneHigh$threshold = "2"
Qs_Threshold3_WCVI_NoneHigh$threshold = "3"

# add region
Qs_Threshold1_WCVI_NoneHigh$region = "WCVI"
Qs_Threshold2_WCVI_NoneHigh$region = "WCVI"
Qs_Threshold3_WCVI_NoneHigh$region = "WCVI"

# combine dfs
Qs_AllThresholds_WCVI_NoneHigh = rbind(Qs_Threshold1_WCVI_NoneHigh, 
                                       Qs_Threshold2_WCVI_NoneHigh, 
                                       Qs_Threshold3_WCVI_NoneHigh)

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

## WCVI No Otters plot of Q's across thresholds
# filter df for WCVI No Otters (contrast 6)
Qs_AllThresholds_WCVI_None = Qs_AllThresholds_WCVI_NoneHigh %>% filter(contrast == "6")
ggplot(aes(y = mean, x = factor(to), color = factor(threshold)), data = Qs_AllThresholds_WCVI_None)+
  geom_point(position = position_dodge(0.9)) + 
  scale_color_manual(name = "Threshold", values = c("1" = "cyan", "2" = "purple", "3" = "orange")) +
  #geom_errorbar(aes(ymin = mean-stderror, ymax = mean+stderror), size = 1, width = 0, position = position_dodge(0.9)) +
  #geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.9)) +
  #geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.9)) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("Transition To This State")+
  ylab("Transition Intensity (Q)") + 
  ggtitle("WCVI No Otters") +
  gg_options()

## WCVI High Otters plot of Q's across thresholds
# filter df for WCVI High Otters (contrast 9)
Qs_AllThresholds_WCVI_High = Qs_AllThresholds_WCVI_NoneHigh %>% filter(contrast == "9")
ggplot(aes(y = mean, x = factor(to), color = factor(threshold)), data = Qs_AllThresholds_WCVI_High)+
  geom_point(position = position_dodge(0.9)) + 
  scale_color_manual(name = "Threshold", values = c("1" = "cyan", "2" = "purple", "3" = "orange")) +
  #geom_errorbar(aes(ymin = mean-stderror, ymax = mean+stderror), size = 1, width = 0, position = position_dodge(0.9)) +
  #geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.9)) +
  #geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.9)) +
  facet_wrap(~from)+ # *** [ facet_wrap doesn't allow label re-naming??? ] ***
  scale_x_discrete(labels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence")) +
  xlab("Transition To This State")+
  ylab("Transition Intensity (Q)") + 
  ggtitle("WCVI High Otters") +
  gg_options() #+ coord_cartesian(ylim = c(-40, 30))


## Pmats

# Load Pmats
pmats = readRDS("~/Markov/Markov_Models/Pmats_Threshold1.rds")
dim(pmats)
pmat_df <- melt(pmats)
names(pmat_df) <- c("iter","contrast","from","to","Prob")

# First, get mean, median, and upper & lower 90% and 95% intervals
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

# *** Plot of just self-transitions (2 panel plot, no and high otter categories, lines denote regions)
# Make 2 figures, use ggarrange to make 2 panel plot
# Panel 1 (no otters) = contrasts 1, 5, 6
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

# re-order comm_state for x-axis and region_OOT for legend
Pmat_selfs$comm_state = factor(Pmat_selfs$comm_state, levels = c("Urchins","Understory","Epibenthic","Macro","Mixed","Coexistence"))
Pmat_selfs$region = factor(Pmat_selfs$region, levels = c("CCBC", "HG", "WCVI"))

# *** Just CCBC & WCVI since Panel 2 is just CCBC & WCVI?
# Panel 1: filter Pmat_selfs df for contrasts with no otters (contrasts 1, 5, 6)
Pmat_OOT_None = Pmat_selfs %>% filter(contrast == "1" | contrast == "5" | contrast == "6")
Panel1 = ggplot(Pmat_OOT_None, aes(comm_state, mean, color = region)) + 
  geom_point(position = position_dodge(0.9)) + 
  scale_color_manual(name = "Region", values = c("CCBC" = "cyan", "HG" = "purple", "WCVI" = "orange")) +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), size = 1, width = 0, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.9)) +
  xlab("Community State")+ 
  ylab("Predicted Self-Transition Probability") + 
  ggtitle("No Otters") +
  gg_options()

# Panel 2: filter Pmat_selfs df for contrasts with high otters (contrasts 4, 9)
Pmat_OOT_High = Pmat_selfs %>% filter(contrast == "4" | contrast == "9")
Panel2 = ggplot(Pmat_OOT_High, aes(comm_state, mean, color = region)) + 
  geom_point(position = position_dodge(0.5)) + 
  scale_color_manual(name = "Region", values = c("CCBC" = "cyan", "WCVI" = "orange")) + 
  labs(color = "Your Title Here") +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), size = 1, width = 0, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pct90_lower, ymax = pct90_upper), size = 1, width = 0, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = pct95_lower, ymax = pct95_upper), size = 0.5, width = 0, position = position_dodge(0.5)) +
  xlab("Community State")+ 
  ylab("Predicted Self-Transition Probability") + 
  ggtitle("High Otters") +
  gg_options()

library("ggpubr")
self_transitions = ggarrange(Panel1, Panel2,
                             ncol = 1, nrow = 2)
self_transitions

# save as PDF
#ggsave("Threshold1.pdf", self_transitions, width=11, height=8.5)
