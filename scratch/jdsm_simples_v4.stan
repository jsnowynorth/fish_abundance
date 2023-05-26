
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
  int<lower=0> p_phi; // number of catchability parameters
  
  // data
  int<lower=0> Y[N, K]; // total caught
  int<lower=0> y_vec[N*K]; // vector of Y
  matrix[N, K] effort; // effort
  matrix[N, p_beta] X; // design matrix abundance
  matrix[N, p_phi] Z; // design matrix catchability
  matrix[Nstar, p_phi] Zstar; // design matrix catchability
  
  // indexing
  int n_lakes; // number of lakes
  int ngears;
  matrix[N, n_lakes] each_lake; // mapping for spatial random effect
  
}

transformed data{

  vector[K] A =  rep_vector(1, K);
  vector[ngears] alpha =  rep_vector(1, ngears);

}

parameters {
  
  // vector[p_beta*K - 1] beta_eff; // relative abundance parameters
  matrix[p_beta, K] beta; // relative abundance parameters
  // simplex[ngears] theta_simplex[K]; // catchability parameters
  matrix[p_phi, K] phi; // catchability parameters
  
  cholesky_factor_corr[K] Sigma_species; // species correlation
  vector<lower=0>[K] tau; // species scale
  matrix[n_lakes, K] z; // random effect sample

}  

transformed parameters{

  // conditional sample mean zero
  cov_matrix[K] Sig = diag_pre_multiply(tau, Sigma_species) * diag_pre_multiply(tau, Sigma_species)';
  matrix[n_lakes, K] zcor = (cholesky_decompose(Sig) * z')';
  matrix[n_lakes, K] omega = zcor - rep_matrix((inverse(Sig) * A * (1/(sum(inverse(Sig)*n_lakes))) * sum(zcor))', n_lakes);

  // row_vector[K] beta_0 = append_col(beta_0_eff, 0); // relative abundance parameters
  // vector[p_beta*K] beta_tmp = append_row(0, beta_eff); // relative abundance parameters
  // matrix[p_beta, K] beta = to_matrix(beta_tmp, p_beta, K); // relative abundance parameters

  // matrix[ngears, K] theta; // catchability parameters
  // for (i in 1:K) theta[,i] = (ngears) * theta_simplex[i]; // catchability parameters


}

model {
  
  // define model parameters
  // matrix[N, K] THETA; // bounded effort scaling
  matrix[N, K] gamma; // relative abundance
  matrix[N, K] OMEGA; // lake-wise random effect
  // matrix[N, K] Etilde;   // catchability intercept
  vector[N * K] lambda;
  
  // catchability scaling centering
  vector[K] mux; // mean for log normal
  vector[K] sx; // sd for log normal
  matrix[N, K] phiz; // unbounded effort scaling
  matrix[N, K] phiz_star; // unbounded effort scaling
  matrix[Nstar, K] phiz_star_star; // unbounded effort scaling
  
  // random effect mapped
  OMEGA = each_lake * omega; // sample mean zero
  
  // relative abundance
  gamma = X * beta + OMEGA;
  // gamma = X * beta;
  
  // temporal catchability
  // phiz = Z * phi;
  phiz_star = Z * phi;
  phiz_star_star = Zstar * phi;
  for (i in 1:K){
    mux[i] = mean(exp(phiz_star_star[,i]))^2;
    sx[i] = variance(exp(phiz_star_star[,i]));
    phiz[,i] = phiz_star[,i] - log(mux[i]/sqrt(mux[i]+sx[i])) - log(1 + sx[i]/mux[i])/2;
  } // mean 1 by species
  
  // lambda = to_vector(gamma);
  // lambda = to_vector(phiz + gamma);
  lambda = to_vector(log(effort) + phiz + gamma);
  // lambda = to_vector(effort .* exp(phiz) .* exp(gamma));
  
  
  // catchability intercept
  // Etilde = block_multiply_matrix(effort, theta);

  // effort scaling
  // THETA = extend_matrix(theta, N);

  // priors
  to_vector(beta) ~ std_normal(); // beta prior
  // to_vector(beta_eff) ~ normal(0, 10); // beta prior
  to_vector(phi) ~ std_normal(); // beta prior
  // for (i in 1:K) theta_simplex[i] ~ dirichlet(alpha); // 1/(p_phi * K)
  to_vector(z) ~ std_normal(); // omega prior
  tau ~ cauchy(0, 2.5); // scale prior
  Sigma_species ~ lkj_corr_cholesky(1); // Sigma_species prior, parameter determines strength of correlation, https://distribution-explorer.github.io/multivariate_continuous/lkj.html
  
  
  // likelihood
  y_vec ~ poisson_log(lambda);
  // y_vec ~ poisson_log(to_vector(phiz + gamma));
  // y_vec ~ poisson_log(to_vector(log(Etilde) + phiz + gamma));

  
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


