
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
  int<lower=0> K; // number of species
  int<lower=0> p_beta; // number of abundance parameters
  int<lower=0> p_phi; // number of catchability parameters
  
  // data
  int<lower=0> Y[N, K]; // total caught
  int<lower=0> y_vec[N*K]; // vector of Y
  matrix[N, K] effort; // effort
  matrix[N, p_beta] X; // design matrix abundance
  matrix[N, p_phi] Z; // design matrix catchability
  
  int n_lakes; // number of lakes
  matrix[N, n_lakes] each_lake; // mapping for spatial random effect
  // matrix[n_lakes, n_lakes] Sigma_spatial; // Spatial correlation matrix
  // matrix[n_lakes, n_lakes] P; // Spatial projection matrix
  
  
}

// transformed data{
// 
//   matrix[n_lakes, n_lakes] Sigma_spatial_chol = cholesky_decompose(Sigma_spatial); // Spatial correlation matrix
// 
// }

parameters {
  
  vector[K] beta_0; // intercept
  matrix[p_beta, K] beta; // abundance parameters
  matrix[p_phi, K] phi; // catchability parameters
  corr_matrix[K] Sigma_species; // species correlation
  vector<lower=0>[K] tau; // species scale
  // vector[n_lakes*K] z; // random effect sample - the matt
  matrix[n_lakes, K] z; // random effect sample - the matt
  
}

transformed parameters{
  
  matrix[n_lakes, K] omega_star; // species random effect
  
  // omega = to_matrix(kronecker_prod(cholesky_decompose(quad_form_diag(Sigma_species, tau)), Sigma_spatial_chol) * z, n_lakes, K);
  omega_star = z * quad_form_diag(Sigma_species, tau);
  
}

model {
  
  // transform variables
  matrix[N, K] B0; // intercept
  matrix[N, K] theta_star; // unbounded catchability
  matrix[N, K] theta; // bounded catchability
  matrix[N, K] gamma; // abundance
  matrix[N, K] lambda; // mean intensity function
  matrix[N, K] OMEGA; // lake-wise random effect
  // vector[n_lakes*K] zeros = rep_vector(0, n_lakes*K);
  vector[N*K] lambda_vec; // lambda vector
  matrix[n_lakes, K] omega; // species random effect

  // define intensity variables
  B0 = rep_matrix(beta_0, N)';
  
  // omega = omega_star - mean(omega_star) - variance(omega_star)/2;
  omega = omega_star - mean(omega_star);
  OMEGA = each_lake * omega;

  theta_star = Z * phi;
  // theta_star_star = Zstar * phi # zstar is all temperature doy sin cos stuff
  theta = theta_star - mean(theta_star) - variance(theta_star)/2;
  // theta = theta_star - mean(theta_star_star) - variance(theta_star_star)/2; # mean 0 for all days
  
  gamma = B0 + X * beta + OMEGA;
  lambda = log(effort) + gamma + theta;
  lambda_vec = to_vector(lambda);
  
  // priors - https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  beta_0 ~ normal(0, 10); // beto 0 prior
  to_vector(beta) ~ std_normal(); // beta prior
  to_vector(phi) ~ std_normal(); // phi prior
  tau ~ cauchy(0, 2.5); // scale prior
  Sigma_species ~ lkj_corr(1); // Sigma_species prior, parameter determines strength of correlation, https://distribution-explorer.github.io/multivariate_continuous/lkj.html
  to_vector(z) ~ std_normal(); // omega prior
  
  y_vec ~ poisson_log(lambda_vec);
  
}
