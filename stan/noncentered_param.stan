data {
  int<lower=0> N;
  int<lower=0> K;
  real<lower=0> GDD_base;
  
  array[N] int<lower=1, upper=K> sp_id;
  vector[N] delta_tmp;
  vector[N] delta_pre;
  vector[N] GDD;
  
  vector[N] log_rw;
}

parameters {
  vector[K] alpha;
  real m_alpha;
  real<lower=0> s_alpha;
  
  vector[K] beta_tilde_GDD;
  real m_beta_GDD;
  real<lower=0> s_beta_GDD;
  
  real beta_tmp;
  real beta_pre;
  real<lower=0> sigma;
}

transformed parameters{
  vector[K] beta_GDD;
  beta_GDD = m_beta_GDD + s_beta_GDD * beta_tilde_GDD;
  
  vector[N] y;
  for (n in 1:N){
    y[n] = alpha[sp_id[n]] + beta_GDD[sp_id[n]]*(GDD[n]-GDD_base) + beta_tmp*delta_tmp[n] + beta_pre*delta_pre[n];
  }
}

model {
  m_alpha ~ normal(0.5, 0.1);
  s_alpha ~ normal(0, 0.15);
  m_beta_GDD ~ normal(0, 0.05);
  s_beta_GDD ~ normal(0, 0.12);
  
  alpha ~ normal(m_alpha, s_alpha);
  beta_tilde_GDD ~ normal(0, 1);
  
  beta_tmp ~ normal(-0.1, 0.3);
  beta_pre ~ normal(-0.1, 0.3);
  sigma ~ normal(0, 0.5);
  
  log_rw ~ normal(y, sigma);
}

generated quantities {
  vector[N] log_rwSim;
  for (i in 1:N) {
    log_rwSim[i] = normal_rng(y[i], sigma);
  }
}




