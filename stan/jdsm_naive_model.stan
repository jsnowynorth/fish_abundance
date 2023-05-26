
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
  matrix extend_matrix(matrix A, int N){
    
    int n = rows(A);
    int K = cols(A);
    matrix[N,K] C;
    int count = 1;

    while(count <= N){
        C[count:(count+n-1),] = A;
    count = count + n;
    }
    
    return C;
  }
  
  
  matrix block_multiply_matrix(matrix A, matrix B){
    
    int N = rows(A);
    int K = cols(A);
    int n = rows(B);
    matrix[N,K] C;
    int count = 1;

    while(count <= N){
      for(k in 1:K){
        C[count:(count+n-1), k] = A[count:(count+n-1),k] .* B[,k];
      }
      count = count + n;
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
  
  // data
  int<lower=0> Y[N, K]; // total caught
  int<lower=0> y_vec[N*K]; // vector of Y
  matrix[N, K] effort; // effort
  matrix[N, p_beta] X; // design matrix abundance
  
  // indexing
  int n_lakes; // number of lakes
  matrix[N, n_lakes] each_lake; // mapping for spatial random effect
  
}

transformed data{

  vector[K] A =  rep_vector(1, K);

}

parameters {
  
  matrix[p_beta, K] beta; // relative abundance parameters
  
  cholesky_factor_corr[K] Sigma_species; // species correlation
  vector<lower=0>[K] tau; // species scale
  matrix[n_lakes, K] z; // random effect sample

}  

transformed parameters{

  // conditional sample mean zero
  cov_matrix[K] Sig = diag_pre_multiply(tau, Sigma_species) * diag_pre_multiply(tau, Sigma_species)';
  matrix[n_lakes, K] zcor = (cholesky_decompose(Sig) * z')';
  matrix[n_lakes, K] omega = zcor - rep_matrix((inverse(Sig) * A * (1/(sum(inverse(Sig)*n_lakes))) * sum(zcor))', n_lakes);



}

model {
  
  // define model parameters
  matrix[N, K] gamma; // relative abundance
  matrix[N, K] OMEGA; // lake-wise random effect
  vector[N * K] lambda;

  // random effect mapped
  OMEGA = each_lake * omega; // sample mean zero
  
  // relative abundance
  gamma = X * beta + OMEGA;
  

  lambda = to_vector(log(effort)  + gamma);


  // priors
  to_vector(beta) ~ std_normal(); // beta prior
  to_vector(z) ~ std_normal(); // omega prior
  tau ~ cauchy(0, 2.5); // scale prior
  Sigma_species ~ lkj_corr_cholesky(1); // Sigma_species prior, parameter determines strength of correlation, https://distribution-explorer.github.io/multivariate_continuous/lkj.html
  
  
  // likelihood
  y_vec ~ poisson_log(lambda);

  
}

// generated quantities{
// 
//   matrix[N, K] THETA; // bounded effort scaling
//   matrix[N, K] gamma; // relative abundance
//   matrix[N, K] lambda; // mean intensity function
// 
//   THETA = (ngears) * extend_matrix(theta, N);
//   gamma = X * beta;
//   lambda = log(effort .* THETA) + gamma;
// 
// }