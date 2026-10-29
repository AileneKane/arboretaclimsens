functions {
  vector gp_pred_rng(real[] x2, vector y1, real[] x1, vector mu,
                     real alpha, real rho, real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K = gp_exp_quad_cov(x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1 - mu);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = gp_exp_quad_cov(x1, x2, alpha, rho);
      vector[N2] f2_mu = mu + (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 = gp_exp_quad_cov(x2, alpha, rho) - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}

data {
  int<lower=1> N_obs; # Number of total observations
  int<lower=1> K; # Number of species
  int<lower=1> N_tree; # Number of trees (individuals) 
  int<lower=1> N_all_year_obs; # Number of Years 
  real<lower=0> GDD_base; 
  
  vector[N_obs] delta_tmp;
  vector[N_obs] delta_pre;
  vector[N_obs] GDD;
  array[N_obs] int<lower=1, upper=K> sp; # Species listed for every obs
  array[N_tree] int<lower=1, upper=K> sp_tree; # Species for every tree
  array[N_tree] int<lower=1> start;
  array[N_tree] int<lower=1> end;
  
  array[N_tree] int<lower=1, upper=N_obs> N_obs_tree; # Number of Observations per Tree
  array[N_all_year_obs] int all_year_obs; # List of all Years Used
  
  # Year and rw observations
  real year_obs[N_obs]; 
  vector[N_obs] log_rw_obs;
}

parameters {
  // Main parameters for linear regression
  vector[K] alpha;
  real m_alpha;
  real<lower=0> s_alpha;
  
  vector[K] beta_tilde_GDD;
  real m_beta_GDD;
  real<lower=0> s_beta_GDD;
  
  real beta_tmp;
  real beta_pre;
  real<lower=0> sigma;
  
  // Hyperparameters for GP
  vector[K] log_rho;
  //vector[K] log_tilde_rho;
  real m_log_rho;
  real<lower=0> s_log_rho;
  
  vector[K] log_gamma;
  //vector[K] log_tilde_gamma;
  real m_log_gamma;
  real<lower=0> s_log_gamma;
}

transformed parameters{
  vector[K] beta_GDD;
  beta_GDD = m_beta_GDD + s_beta_GDD * beta_tilde_GDD;
  
  vector[K] rho;
  rho = exp(log_rho);
  vector[K] gamma;
  gamma = exp(log_gamma);
  
  vector[N_obs] y_clim;
  for (n in 1:N_obs){
    y_clim[n] = alpha[sp[n]] + beta_GDD[sp[n]]*(GDD[n]-GDD_base) + beta_tmp*delta_tmp[n] + beta_pre*delta_pre[n];
  }
}

model {
  // Parameters with partial pooling
  m_alpha ~ normal(0.5, 0.1);
  s_alpha ~ normal(0, 0.15);
  m_beta_GDD ~ normal(0, 0.05);
  s_beta_GDD ~ normal(0, 0.12);
  
  alpha ~ normal(m_alpha, s_alpha);
  beta_tilde_GDD ~ normal(0, 1);
  
  // Parameters without partial pooling
  beta_tmp ~ normal(-0.1, 0.3);
  beta_pre ~ normal(-0.1, 0.3);
  sigma ~ normal(0, 0.5);
  
  // Hyperparameters with partial pooling
  m_log_rho ~ normal(1, 0.2);
  s_log_rho ~ normal(0, 0.3);
  m_log_gamma ~ normal(-1, 0.2);
  s_log_gamma ~ normal(0, 0.3);
  
  //log_tilde_rho ~ normal(0, 1);
  //log_tilde_gamma ~ normal(0, 1)
  log_rho ~ normal(m_log_rho, s_log_rho);
  log_gamma ~ normal(m_log_gamma, s_log_gamma);
  
  // Define an array of cholesky decomposed covariance matrices
  array[K] matrix[N_all_year_obs, N_all_year_obs] L_cov;
  
  // Calculate covariance matrix for every species
  for (species in 1:K) {
    matrix[N_all_year_obs, N_all_year_obs] cov = gp_exp_quad_cov(all_year_obs, exp(log_gamma[species]), 
    exp(log_rho[species])) + diag_matrix(rep_vector(square(sigma), N_all_year_obs));
    L_cov[species] = cholesky_decompose(cov);
  }
  
  # Subset L_cov for each tree
    for (i in 1:N_tree){
      matrix[N_obs_tree[i], N_obs_tree[i]] L_cov_tree = block(L_cov[sp_tree[i]], 1, 1, N_obs_tree[i], N_obs_tree[i]);
      log_rw_obs[start[i]:end[i]] ~ multi_normal_cholesky(y_clim[start[i]:end[i]], L_cov_tree);
    }
}

generated quantities {
  vector[N_obs] f;
  array[N_obs]real log_rw_pred;
  
  for (i in 1:N_tree){
    f[start[i]:end[i]] = gp_pred_rng(year_obs[start[i]:end[i]], log_rw_obs[start[i]:end[i]], 
    year_obs[start[i]:end[i]], y_clim[start[i]:end[i]], gamma[sp_tree[i]], rho[sp_tree[i]], sigma, 1e-10);
  
    log_rw_pred[start[i]:end[i]] = normal_rng(f[start[i]:end[i]], sigma);
  }
}
