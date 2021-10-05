
data {
  
  // dimension parameters
  int<lower=0> N; // number of observations
  int<lower=0> K; // number of species
  int<lower=0> p_beta; // number of abundance parameters
  int<lower=0> p_phi; // number of catchability parameters
  
  // data
  int<lower=0> Y[N, K]; // total caught
  int<lower=0> y_vec[N*K]; // vector of Y
  matrix[N, K] effort; // effort
  matrix[N, p_beta] X; // design matrix abundance
  matrix[N, p_phi] Z; // design matrix catchability
  
}


parameters {
  
  vector[K] beta_0; // intercept
  // real beta_0; // intercept
  matrix[p_beta, K] beta; // abundance parameters
  matrix[p_phi, K] phi; // catchability parameters
  vector[K] omega; // species random effect
  
  corr_matrix[K] Sigma_species; // species correlation
  // cholesky_factor_corr[K] Lcorr;
  vector<lower=0>[K] tau; // species scale
  
}


model {
  
  // transform variables
  matrix[N, K] B0; // intercept
  matrix[N, K] theta_star; // unbounded catchability
  matrix[N, K] theta; // bounded catchability
  matrix[N, K] gamma; // abundance
  matrix[N, K] lambda; // mean intensity function
  matrix[N, K] OMEGA; // larger omega
  vector[N*K] lambda_vec; // lambda vector
  vector[K] zeros = rep_vector(0, K);
  
  // define intensity variables
  OMEGA = rep_matrix(omega, N)'; // map omega to full design matrix
  B0 = rep_matrix(beta_0, N)';
  
  gamma = B0 + X * beta + OMEGA;
  
  theta_star = Z * phi;
  theta = theta_star - mean(theta_star);
  
  lambda = log(effort) + gamma + theta;
  
  
  // priors
  beta_0 ~ normal(0, 10); // beto 0 prior
  to_vector(beta) ~ normal(0, 10); // beta prior
  to_vector(phi) ~ normal(0, 10); // phi prior
  
  tau ~ cauchy(0, 2.5); // scale prior
  Sigma_species ~ lkj_corr(1); // Sigma_species prior
  // Lcorr ~ lkj_corr_cholesky(1); // Sigma_species prior
  omega ~ multi_normal(zeros, quad_form_diag(Sigma_species, tau)); // omega prior
  // omega ~ multi_normal(zeros, diag_pre_multiply(tau, Lcorr)); // omega prior
 

  // log likelihood
  lambda_vec = to_vector(lambda);
  y_vec ~ poisson_log(lambda_vec);
  
}
