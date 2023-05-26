
functions{
  
  // #include ssm.stan
  matrix kronecker_prod(matrix A, matrix B) {
    matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
    int m;
    int n;
    int p;
    int q;
    m = rows(A);
    n = cols(A);
    p = rows(B);
    q = cols(B);
    for (i in 1:m) {
      for (j in 1:n) {
        int row_start;
        int row_end;
        int col_start;
        int col_end;
        row_start = (i - 1) * p + 1;
        row_end = (i - 1) * p + p;
        col_start = (j - 1) * q + 1;
        col_end = (j - 1) * q + q;
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
    return C;
  }
  
  matrix to_matrix_colwise(vector v, int m, int n) {
    matrix[m, n] res;
    for (j in 1:n) {
      for (i in 1:m) {
        res[i, j] = v[(j - 1) * m + m];
      }
    }
    return res;
  }
  
  matrix to_matrix_rowwise(vector v, int m, int n) {
    matrix[m, n] res;
    for (i in 1:n) {
      for (j in 1:m) {
        res[i, j] = v[(i - 1) * n + n];
      }
    }
    return res;
  }
  
  // A is input matrix of smaller dimension
  // N is number of output rows
  // each is number of times to repeat each row of A, integer vector of different numbers
  matrix extend_matrix(matrix A, int N, int[] each){
    int n = rows(A);
    int K = cols(A);
    matrix[N,K] C;
    int count = 1;

    for(i in 1:n){
      for(j in 1:each[i]){
        C[count,] = A[i,];
        count += 1;
      }
    }
    return C;
  }
  
  
  
}

data {
  
  // dimension parameters
  int<lower=0> N; // number of observations
  int<lower=0> Nstar; // number of observations
  int<lower=0> K; // number of species
  int<lower=0> p_beta; // number of abundance parameters
  int<lower=0> p_phi; // number of catchability parameters
  
  // data
  int<lower=0> Y[N, K]; // total caught
  int<lower=0> y_vec[N*K]; // vector of Y
  matrix[N, K] effort; // effort
  matrix[N, p_beta] X; // design matrix abundance
  matrix[N, p_phi] Z; // design matrix catchability
  matrix[Nstar, p_phi] Zstar; // design matrix catchability
  
  int n_lakes; // number of lakes
  matrix[N, n_lakes] each_lake; // mapping for spatial random effect
  
}

transformed data{

  // row_vector[n_lakes] mean_ones = rep_row_vector(1, n_lakes);
  // vector[K * n_lakes] A =  rep_vector(1, n_lakes * K);
  // matrix[n_lakes, n_lakes] Qpre = diag_matrix(rep_vector(1, n_lakes));
  
  vector[K] A =  rep_vector(1, K);
  // matrix[K, K] Sig = diag_matrix(A);

}

parameters {
  
  // row_vector[K] beta_0; // intercept
  // row_vector[K-1] beta_0_eff; // intercept
  // matrix[p_beta, K] beta; // abundance parameters
  vector[p_beta*K - 1] beta_eff; // catchability parameters
  matrix[p_phi, K] phi; // catchability parameters
  // vector[p_phi*K - 1] phi_eff; // catchability parameters
  cholesky_factor_corr[K] Sigma_species; // species correlation
  vector<lower=0>[K] tau; // species scale
  matrix[n_lakes, K] z; // random effect sample - the matt
  
}

transformed parameters{

  // force mean zero
  // matrix[n_lakes, K] omega = z * diag_pre_multiply(tau, Sigma_species); // species random effect
  
  // conditional sample mean zero
  matrix[K, K] Sig = diag_pre_multiply(tau, Sigma_species) * diag_pre_multiply(tau, Sigma_species)';
  matrix[n_lakes, K] omega = z - rep_matrix((inverse(Sig) * A * (1/(sum(inverse(Sig)*n_lakes))) * sum(z))', n_lakes);
  
  // row_vector[K] beta_0 = append_col(beta_0_eff, 0); // catchability parameters
  vector[p_beta*K] beta_tmp = append_row(0, beta_eff); // catchability parameters
  matrix[p_beta, K] beta = to_matrix(beta_tmp, p_beta, K); // catchability parameters

  // phi = to_matrix(phi_tmp, p_phi, K);

  // omega_star = z * quad_form_diag(Sigma_species, tau);


}

model {
  
  // define model parameters

  // matrix[N, K] B0; // global intercept
  real mux; // mean for log normal
  real sx; // sd for log normal
  matrix[N, K] theta_star; // unbounded effort scaling
  matrix[Nstar, K] theta_star_star; // unbounded effort scaling
  matrix[N, K] theta; // bounded effort scaling
  matrix[N, K] gamma; // relative abundance
  matrix[N, K] lambda; // mean intensity function
  matrix[N, K] OMEGA; // lake-wise random effect
  vector[N*K] lambda_vec; // lambda vector
  // matrix[n_lakes, K] omega; // species random effect
  // matrix[n_lakes * K, n_lakes * K] Q = kronecker_prod(inverse(diag_pre_multiply(tau, Sigma_species)), Qpre); // species random effect

  // Global intercept
  // B0 = rep_matrix(beta_0, N);
  
  // Random effects - mean zero by species
  // omega_mean = (mean_ones * omega_star)./n_lakes; // calculate mean of omega
  // omega = omega_star - rep_matrix(omega_mean, n_lakes);
  // OMEGA = each_lake * omega;
  
  // omega = each_lake * omega_star;
  // OMEGA = each_lake * (omega - mean(omega)); // force mean zero
  OMEGA = each_lake * omega; // sample mean zero
  
  // omega = z - to_matrix(Q * A * (1/(A' * Q * A)) * sum(to_vector(z)), n_lakes, K);
  // OMEGA = each_lake * omega;
  
  // relative abundance
  // gamma = B0 + X * beta + OMEGA;
  gamma = X * beta + OMEGA;

  // effort scaling
  // theta = Z * phi;
  theta_star = Z * phi;
  theta_star_star = Zstar * phi;
  mux = mean(exp(theta_star_star))^2;
  sx = variance(exp(theta_star_star));
  theta = theta_star - log(mux/sqrt(mux+sx)) - log(1 + sx/mux)/2;
  
  // log lambda
  lambda = log(effort) + gamma + theta;
  // lambda = log(effort) + gamma;
  lambda_vec = to_vector(lambda);
  
  // priors - https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  // beta_0 ~ normal(0, 10); // beto 0 prior
  // beta_0_eff ~ std_normal(); // beto 0 prior
  // to_vector(beta) ~ std_normal(); // beta prior
  // to_vector(phi) ~ std_normal(); // phi prior
  // to_vector(beta) ~ normal(0, 10); // beta prior
  to_vector(beta_eff) ~ normal(0, 10); // beta prior
  to_vector(phi) ~ normal(0, 10); // phi prior
  // phi_eff ~ std_normal(); // phi prior
  to_vector(z) ~ std_normal(); // omega prior
  tau ~ cauchy(0, 2.5); // scale prior
  Sigma_species ~ lkj_corr_cholesky(0.5); // Sigma_species prior, parameter determines strength of correlation, https://distribution-explorer.github.io/multivariate_continuous/lkj.html
  
  
  // likelihood
  y_vec ~ poisson_log(lambda_vec);
  
}

// generated quantities{

  
  // matrix[n_lakes, K] omega; // species random effect
  
  // return random effects mean zero by species
  // omega = z - to_matrix(Q * A' * (1/(A * Q * A')) * sum(to_vector(z)), n_lakes, K);

// }


