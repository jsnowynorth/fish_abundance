
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
  
  // matrix extend_matrix(matrix A, int N, int[] each){
  //   int n = rows(A);
  //   int K = cols(A);
  //   matrix[N,K] C;
  //   // int count = 1;
  //   
  //   for(i in 1:n){
  //     C = append_row(C, rep_matrix(A[i,], each[i]));
  //     // for(j in 1:each[i]){
  //     //   C[count,] = A[i,];
  //     //   count += 1;
  //     // }
  //   }
  //   return C;
  // }
  

}

data {
  
  // dimension parameters
  int<lower=0> N; // number of observations
  int<lower=0> K; // number of species
  int<lower=0> p_beta; // number of abundance parameters
  int<lower=0> p_phi; // number of catchability parameters
  
  // data
  int<lower=0> Y[N, K]; // total caught
  matrix[N, K] effort; // effort
  matrix[N, p_beta] X; // design matrix abundance
  matrix[N, p_phi] Z; // design matrix catchability
  
  int n_lakes; // number of lakes
  int each_lake[n_lakes]; // mapping for spatial random effect
  matrix[n_lakes, n_lakes] Sigma_spatial; // Spatial correlation matrix
  
  
}

transformed data{
  
  matrix[n_lakes, n_lakes] Sigma_spatial_chol = cholesky_decompose(Sigma_spatial); // Spatial correlation matrix
  
}

parameters {
  
  // vector[K] beta_0; // intercept
  
  matrix[p_beta, K] beta; // abundance parameters
  matrix[p_phi, K] phi; // catchability parameters
  matrix[n_lakes, K] omega; // species random effect
  corr_matrix[K] Sigma_species; // species correlation
  vector<lower=0>[K] tau; // species scale
  
}

model {
  
  // transform variables
  matrix[N, K] theta_star; // unbounded catchability
  matrix[N, K] theta; // bounded catchability
  matrix[N, K] gamma; // abundance
  matrix[N, K] lambda; // mean intensity function
  matrix[N, K] OMEGA; // lake-wise random effect
  vector[n_lakes*K] zeros = rep_vector(0, n_lakes*K);
  
  // matrix[n_lakes*K, n_lakes*K] L_K; // cholesky decomp
  matrix[n_lakes*K, n_lakes*K] Sigma; // sigma species

  OMEGA = extend_matrix(omega, N, each_lake);
  // int count = 1;
  // for(i in 1:n_lakes){
  //   for(j in 1:each_lake[i]){
  //     OMEGA[count,] = omega[i,];
  //     count += 1;
  //   }
  // }
  
  // define intensity variables
  gamma = X * beta + OMEGA;
  
  theta_star = Z * phi;
  theta = theta_star - mean(theta_star);
  
  lambda = log(effort) + gamma + theta;
  
  
  Sigma = kronecker_prod(cholesky_decompose(quad_form_diag(Sigma_species, tau)), Sigma_spatial_chol);
  
  // priors
  // beta_0 ~ normal(0, 10); // beto 0 prior
  to_vector(beta) ~ normal(0, 1); // beta prior
  to_vector(phi) ~ normal(0, 1); // phi prior
  
  tau ~ cauchy(0, 2.5); // scale prior
  Sigma_species ~ lkj_corr(1); // Sigma_species prior, parameter determines strength of correlation, https://distribution-explorer.github.io/multivariate_continuous/lkj.html
  
  to_vector(omega) ~ multi_normal_cholesky(zeros, Sigma); // omega prior
  

  // log likelihood
  for(k in 1:K){
    Y[,k] ~ poisson_log(lambda[,k]);
  }
  
}




