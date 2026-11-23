library(dplyr)
library(readxl)
library(tidyverse)
library(rstan)

# # Filepaths to rw, climate data
# climData_path <- "/home/ailunf/final_model_data/climaticData.xlsx"
# aRW_path <- "/home/ailunf/final_model_data/Arboringwidthsangiosperms.xlsx"
# gRW_path <- "/home/ailunf/final_model_data/Arboringwidthsgymnosperms.xlsx"
# GDD_path <- "/home/ailunf/final_model_data/GDD_data.xlsx"

climData_path <- "C:/Temporal Ecology Lab/arboretaclimsens/final_model_data/climaticData.xlsx"
aRW_path <- "C:/Temporal Ecology Lab/arboretaclimsens/final_model_data/Arboringwidthsangiosperms.xlsx"
gRW_path <- "C:/Temporal Ecology Lab/arboretaclimsens/final_model_data/Arboringwidthsgymnosperms.xlsx"
GDD_path <- "C:/Temporal Ecology Lab/arboretaclimsens/final_model_data/GDD_data.xlsx"

# Load all data into workspace
prov_clims <- read_excel(climData_path)
a_RW <- read_excel(aRW_path)
g_RW <- read_excel(gRW_path)
GDD_data <- read_excel(GDD_path)
all_cores <- bind_rows(a_RW, g_RW)

# Pre-computed climate metrics for the arboretum (Boston)
arbPre_1950_1979 <- 11.4586002082825
arbPre_1980_2009 <- 12.0016001985868
arbTmp_1950_1979 <- 9.52888905708161
arbTmp_1980_2009 <- 10.0055557200478

GDD_data$GDD <- GDD_data$GDD / 100
GDD_base <- mean(GDD_data$GDD)

# Combine ring widths with provenance data
sorted <- left_join(prov_clims, all_cores, by = "TreeID")
sorted$dT_1950_1979 <- abs(sorted$tmp_1950_1979 - arbTmp_1950_1979)
sorted$dT_1980_2009 <- abs(sorted$tmp_1980_2009 - arbTmp_1980_2009)
sorted$dP_1950_1979 <- abs(sorted$pre_1950_1979 - arbPre_1950_1979)
sorted$dP_1980_2009 <- abs(sorted$pre_1980_2009 - arbPre_1980_2009)

# Subset the select years of data
yr_range <- as.character(seq(2014, 1975, -1))
entries <- c("TreeID", "SciName", "dT_1980_2009", "dP_1980_2009", yr_range)
sample <- sorted %>% select(all_of(entries))
sample_sub <- sample[1:100,]

# Assign a species # for each tree
K <- length(unique(sample_sub$SciName))
sample_sub$sp_tree <- as.numeric(as.factor(sample_sub$SciName))

long_data <- sample_sub %>%
  pivot_longer(
    cols = all_of(yr_range),
    names_to = "Year",                              
    values_to = "rw"                               
  )

# Combine with GDD Data and remove all NA rows
long_data$Year <- as.numeric(long_data$Year)
long_data <- left_join(long_data, GDD_data, by = "Year")
long_data <- na.omit(long_data)

# Set up model data_list
N <- nrow(long_data)
N_tree <- nrow(sample_sub)
N_all_year_obs <- length(yr_range)

rw_years <- subset(sample_sub, select=yr_range)
N_obs_tree <- rowSums(!is.na(rw_years))

end <- cumsum(N_obs_tree)
start <- end - N_obs_tree + 1

# model data_list initialized here
data_list <- list(
  N_obs = N,
  K = K,
  N_tree = N_tree,
  N_all_year_obs = N_all_year_obs, 
  GDD_base = GDD_base,
  
  delta_tmp = long_data$dT_1980_2009,
  delta_pre = long_data$dP_1980_2009,
  GDD = long_data$GDD,
  sp = long_data$sp_tree,
  sp_tree = sample_sub$sp_tree,
  
  start = as.integer(start), 
  end = as.integer(end),
  N_obs_tree = N_obs_tree,
  all_year_obs = as.numeric(yr_range),
  
  year_obs = long_data$Year,
  log_rw_obs = log(long_data$rw, 10)
)

# Avoid recompiling unchanged model
rstan_options(auto_write = TRUE)

setwd("C:/Temporal Ecology Lab/arboretaclimsens/Model")
gp_model <- stan_model('GP_multispecies.stan')
fit1 <- sampling(gp_model, data = data_list, chains = 4, cores = 4, save_warmup = FALSE)

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

diagnostics <- util$extract_hmc_diagnostics(fit1)
util$check_all_hmc_diagnostics(diagnostics)
samples <- util$extract_expectand_vals(fit1)
params <- c("alpha_tilde", "m_alpha", "s_alpha", "beta_tilde_GDD", "m_beta_GDD", "s_beta_GDD", "beta_tmp",
            "beta_pre", "sigma", "log_tilde_rho", "m_log_rho", "s_log_rho", "log_tilde_gamma", "m_log_gamma", "s_log_gamma");
base_samples <- util$filter_expectands(samples, params, check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)


# i <- 9
# idxs <- start[i]:end[i]
# names <- paste0('log_rw_pred[', idxs, ']')
# util$plot_conn_pushforward_quantiles(samples, names, data_list$year_obs[idxs], display_ylim = c(-1,2))
# points(x = data_list$year_obs[idxs], y = data_list$log_rw_obs[idxs], col = "white", cex = 2, pch = 20)
# points(x = data_list$year_obs[idxs], y = data_list$log_rw_obs[idxs], pch = 20)
# 
# summary(fit1, pars=c("rho"))
# 
# util$plot_disc_pushforward_quantiles(samples,paste0('rho[',1:data_list$K, ']'),
#                                  xticklabs =  1:data_list$K,
#                                  ylab = bquote(rho[species]),
#                                  xlab = 'Species')
# 
# util$plot_disc_pushforward_quantiles(samples,paste0('gamma[',1:data_list$K, ']'),
#                                      xticklabs =  1:data_list$K,
#                                      ylab = bquote(gamma[species]),
#                                      xlab = 'Species')
# 
# util$plot_disc_pushforward_quantiles(samples,paste0('alpha[',1:data_list$K, ']'),
#                                      xticklabs =  1:data_list$K,
#                                      ylab = bquote(alpha[species]),
#                                      xlab = 'Species')
# 
# util$plot_disc_pushforward_quantiles(samples,paste0('beta_GDD[',1:data_list$K, ']'),
#                                      xticklabs =  1:data_list$K,
#                                      ylab = bquote(beta_GDD[species]),
#                                      xlab = 'Species')



