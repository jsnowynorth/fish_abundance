
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
  // matrix[N, p_phi] Z; // design matrix catchability
  
}


parameters {
  
  vector[K] beta_0; // intercept
  // real beta_0; // intercept
  matrix[p_beta, K] beta; // abundance parameters
  // matrix[p_phi, K] phi; // catchability parameters
  
}


model {
  
  // transform variables
  matrix[N, K] B0; // intercept
  // matrix[N, K] theta_star; // unbounded catchability
  // matrix[N, K] theta; // bounded catchability
  matrix[N, K] gamma; // abundance
  matrix[N, K] lambda; // mean intensity function
  vector[N*K] lambda_vec; // lambda vector
  
  // define intensity variables
  B0 = rep_matrix(beta_0, N)';
  
  gamma = B0 + X * beta;
  
  // theta_star = Z * phi;
  // theta = theta_star - mean(theta_star);
  
  // lambda = log(effort) + gamma + theta;
  lambda = log(effort) + gamma;
  
  
  // priors
  beta_0 ~ normal(0, 1); // beto 0 prior
  to_vector(beta) ~ normal(0, 1); // beta prior
  // to_vector(phi) ~ normal(0, 1); // phi prior

  

  // log likelihood
  lambda_vec = to_vector(lambda);
  y_vec ~ poisson_log(lambda_vec);
  
}

