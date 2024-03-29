
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
  
  // data
  int<lower=0> Y[N, K]; // total caught
  int<lower=0> y_vec[N*K]; // vector of Y
  matrix[N, K] effort; // effort
  matrix[N, p_beta] X; // design matrix abundance
  
  int n_lakes; // number of lakes
  matrix[N, n_lakes] each_lake; // mapping for spatial random effect
  
  
}

transformed data{

  row_vector[n_lakes] mean_ones = rep_row_vector(1, n_lakes);

}

parameters {
  
  row_vector[K] beta_0; // intercept
  matrix[p_beta, K] beta; // abundance parameters
  corr_matrix[K] Sigma_species; // species correlation
  vector<lower=0>[K] tau; // species scale
  matrix[n_lakes, K] z; // random effect sample - the matt
  
}

transformed parameters{
  
  matrix[n_lakes, K] omega_star; // species random effect
  omega_star = z * quad_form_diag(Sigma_species, tau);
  
}

model {
  
  // transform variables
  matrix[N, K] B0; // intercept
  matrix[N, K] gamma; // abundance
  matrix[N, K] lambda; // mean intensity function
  matrix[N, K] OMEGA; // lake-wise random effect
  vector[N*K] lambda_vec; // lambda vector
  matrix[n_lakes, K] omega; // species random effect
  row_vector[K] omega_mean; // species random effect mean

  // define intensity variables
  B0 = rep_matrix(beta_0, N);
  
  omega_mean = (mean_ones * omega_star)./n_lakes; // calculate mean of omega
  omega = omega_star - rep_matrix(omega_mean, n_lakes);
  OMEGA = each_lake * omega;
  
  gamma = B0 + X * beta + OMEGA;
  lambda = log(effort) + gamma;
  lambda_vec = to_vector(lambda);
  
  // priors - https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  beta_0 ~ normal(0, 10); // beto 0 prior
  to_vector(beta) ~ std_normal(); // beta prior
  tau ~ cauchy(0, 2.5); // scale prior
  Sigma_species ~ lkj_corr(1); // Sigma_species prior, parameter determines strength of correlation, https://distribution-explorer.github.io/multivariate_continuous/lkj.html
  to_vector(z) ~ std_normal(); // omega prior
  
  y_vec ~ poisson_log(lambda_vec);
  
}
