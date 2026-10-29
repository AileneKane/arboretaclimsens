rm(list=ls())
# setwd("C:/Temporal Ecology Lab/arboretaclimsens/Model")
setwd("/home/victor/projects/arboretaclimsens")

library(dplyr)
library(rstan)

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

# Simulate some fake climate and tree ring data
set.seed(82)

# Initial data simulation settings 
N_all_year_obs <- 20 # Number of Years
N_tree <- 10 # Number of trees (individuals) 
all_year_obs <- seq(2014, 2014 - N_all_year_obs + 1, -1) # Range of years observed
N_obs_tree <- sample(10:N_all_year_obs, size = N_tree, replace=TRUE) # num years observed per tree
N_obs <- sum(N_obs_tree) # total num of observations

# GDD data simulation
GDD_vals <- rnorm(N_all_year_obs, mean = 26, sd = 2)
GDD_base <- 26

# Climate data simulation
pre_diff <- abs(rnorm(N_tree, mean = 2.55, sd = 2.44))
tmp_diff <- abs(rnorm(N_tree, mean = 2.8, sd = 2.4))

# Test parameters
alpha <- 0.35
beta_GDD <- 0.01
beta_tmp <- -0.005
beta_pre <- -0.005

gamma <- 0.6
rho <- 2.9
sigma <- 0.4

# Data vectors
delta_pre <- rep(pre_diff, N_obs_tree)
delta_tmp <- rep(tmp_diff, N_obs_tree)

year_obs <- c()
end <- cumsum(N_obs_tree)
start <- end - N_obs_tree + 1

for (n in 1:N_tree){
  year_obs[(length(year_obs)+1):(length(year_obs)+N_obs_tree[n])] <- seq(2014, 2014 - N_obs_tree[n]+1, -1)
}

GDD_data <- data.frame(
  year_obs = all_year_obs,
  GDD = GDD_vals
)

all_mean_data <- data.frame(
  year_obs = year_obs,
  delta_tmp = delta_tmp,
  delta_pre = delta_pre
)

all_mean_data <- left_join(all_mean_data, GDD_data, by= "year_obs")

all_mean_data$mu_log_rw <- rep(alpha, N_obs) + beta_GDD*(all_mean_data$GDD - GDD_base) + 
  beta_tmp*all_mean_data$delta_tmp + beta_pre*all_mean_data$delta_pre

# Simulate the GP now
cov_func <- function(delta_x, gamma, rho){
  return(gamma^2*exp(-1/2*(delta_x/rho)^2))
}

grid <- data.frame(x1 = seq(2014, 2014-N_all_year_obs+1, -1), x2 = seq(2014,2014-N_all_year_obs+1, -1))
dmat <- as.matrix(dist(grid, diag = T, upper = T))

cov <- apply(dmat, 1:2, cov_func, gamma = gamma, rho = rho)
f <- c()

for (i in 1:N_tree){
  f[start[i]:end[i]] <- MASS::mvrnorm(1, all_mean_data$mu_log_rw[start[i]:end[i]], cov[1:N_obs_tree[i],1:N_obs_tree[i]])
}

log_rw_sim <- rnorm(N_obs, mean = f, sd = sigma)

data_list <- list(
  N_obs = N_obs,
  N_tree = N_tree,
  N_all_year_obs = N_all_year_obs, 
  GDD_base = GDD_base,
  
  delta_tmp = all_mean_data$delta_tmp,
  delta_pre = all_mean_data$delta_pre,
  GDD = all_mean_data$GDD,
  start = as.integer(start), 
  end = as.integer(end),
  
  N_obs_tree = N_obs_tree,
  all_year_obs = all_year_obs,
  
  year_obs = all_mean_data$year_obs,
  log_rw_obs = log_rw_sim
)

i <- 5
idxs <- start[i]:end[i]
plot(x =  data_list$year_obs[idxs], y =  f[idxs], type = 'l', ylim = c(-1,3))
points(x = data_list$year_obs[idxs], y = data_list$log_rw_obs[idxs], pch = 20)


gp_model <- stan_model('stan/GP_fit.stan')
fit1 <- sampling(gp_model, data = data_list,
                chains = 4, cores = 4)

diagnostics <- util$extract_hmc_diagnostics(fit1)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit1)
util$check_all_expectand_diagnostics(samples)

i <- 10
idxs <- start[i]:end[i]
names <- paste0('f[', idxs, ']')
util$plot_conn_pushforward_quantiles(samples, names, data_list$year_obs[idxs], display_ylim = c(-1,2))
lines(x = data_list$year_obs[idxs], y = f[idxs], col = "black", lwd = 2)
points(x = data_list$year_obs[idxs], y = data_list$log_rw_obs[idxs], pch = 20)

names <- paste0('log_rw_pred[', idxs, ']')
util$plot_conn_pushforward_quantiles(samples, names, data_list$year_obs[idxs], display_ylim = c(-1,2))
lines(x = data_list$year_obs[idxs], y = f[start[i]:end[i]], col = "white", lwd = 4)
lines(x = data_list$year_obs[idxs], y = f[start[i]:end[i]], col = "black", lwd = 2)
points(x = data_list$year_obs[idxs], y = data_list$log_rw_obs[idxs], col = "white", cex = 2, pch = 20)
points(x = data_list$year_obs[idxs], y = data_list$log_rw_obs[idxs], pch = 20)

util$plot_expectand_pushforward(samples[['rho']], 30, 'rho')
util$plot_expectand_pushforward(samples[['gamma']], 30, 'gamma')
util$plot_expectand_pushforward(samples[['beta_GDD']], 30, 'beta_GDD')
