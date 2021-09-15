matrix to_symmetric_matrix(matrix x) {
  return 0.5 * (x + x ');
}
matrix to_matrix_colwise(vector v, int m, int n) {
  matrix[m, n] res;
  int k;
  k = 1;
  for (j in 1:n) {
    for (i in 1:m) {
      res[i, j] = v[k];
      k = k + 1;
    }
  }
  return res;
}
matrix matrix_pow(matrix A, real n);
matrix matrix_pow(matrix A, real n) {
  real nn;
  nn = floor(n);
  if (nn == 0) {
    return diag_matrix(rep_vector(1., rows(A)));
  } else if (nn == 1) {
    return A;
  } else if (nn > 1) {
    if (fmod(nn, 2.) > 0) {
      return A * matrix_pow(A, nn - 1);
    } else {
      return matrix_pow(A, nn / 2) * matrix_pow(A, nn / 2);
    }
  } else {
    reject("Only non-negative values of n are allowed");
    return A;
  }
}
int symmat_size(int n) {
  int sz;
  sz = 0;
  for (i in 1:n) {
    sz = sz + i;
  }
  return sz;
}
int find_symmat_dim(int n) {
  int i;
  int remainder;
  remainder = n;
  i = 0;
  while (remainder > 0) {
    i = i + 1;
    remainder = remainder - i;
  }
  return i;
}
matrix vector_to_symmat(vector x, int n) {
  matrix[n, n] m;
  int k;
  k = 1;
  for (j in 1:n) {
    for (i in j:n) {
      m[i, j] = x[k];
      if (i != j) {
        m[j, i] = m[i, j];
      }
      k = k + 1;
    }
  }
  return m;
}
vector symmat_to_vector(matrix x) {
  vector[symmat_size(min(rows(x), cols(x)))] v;
  int m;
  int k;
  k = 1;
  m = min(rows(x), cols(x));
  for (j in 1:m) {
    for (i in j:m) {
      v[k] = x[i, j];
      k = k + 1;
    }
  }
  return v;
}
matrix rep_lower_triangular_matrix(real x, int m, int n, int diag) {
  matrix[m, n] A;
  for (i in 1:m) {
    for (j in 1:n) {
      if (i > j) {
        A[i, j] = x;
      } else if (i == j) {
        if (diag) {
          A[i, j] = x;
        } else {
          A[i, j] = 0.;
        }
      } else {
        A[i, j] = 0.;
      }
    }
  }
  return A;
}
matrix rep_upper_triangular_matrix(real x, int m, int n, int diag) {
  matrix[m, n] A;
  for (i in 1:m) {
    for (j in 1:n) {
      if (i < j) {
        A[i, j] = x;
      } else if (i == j) {
        if (diag) {
          A[i, j] = x;
        } else {
          A[i, j] = 0.;
        }
      } else {
        A[i, j] = 0.;
      }
    }
  }
  return A;
}
matrix rep_diagonal_matrix(real x, int m, int n, int k) {
  matrix[m, n] A;
  int mn;
  A = rep_matrix(0., m, n);
  mn = min(m, n);
  if (k >= 0) {
    for (i in 1:min(m, n - k)) {
      A[i, i + k] = x;
    }
  } else {
    for (i in 1:min(m + k, n)) {
      A[i - k, i] = x;
    }
  }
  return A;
}
matrix fill_matrix(matrix x, int m, int n, int[] i, int[] j, real a) {
  matrix[m, n] ret;
  ret = rep_matrix(a, m, n);
  ret[i, j] = x;
  return ret;
}
vector fill_vector(vector x, int n, int[] i, real a) {
  vector[n] ret;
  ret = rep_vector(a, n);
  ret[i] = x;
  return ret;
}
int int_sum_true(int[] x) {
  int n;
  n = 0;
  for (i in 1:num_elements(x)) {
    if (int_step(x[i])) {
      n = n + 1;
    }
  }
  return n;
}
int int_sum_false(int[] x) {
  int n;
  n = 0;
  for (i in 1:num_elements(x)) {
    if (! int_step(x[i])) {
      n = n + 1;
    }
  }
  return n;
}
int[] mask_indexes(int[] x, int n) {
  int idx[n];
  int j;
  j = 1;
  if (n > 0) {
    for (i in 1:num_elements(x)) {
      if (! int_step(x[i]) && j <= n) {
        idx[j] = i;
        j = j + 1;
      }
    }
  }
  return idx;
}
int[] select_indexes(int[] x, int n) {
  int idx[n];
  int j;
  j = 1;
  if (n > 0) {
    for (i in 1:num_elements(x)) {
      if (int_step(x[i]) && j <= n) {
        idx[j] = i;
        j = j + 1;
      }
    }
  }
  return idx;
}
real normal2_rng(real mu, real sigma) {
  real y;
  if (sigma <= 0) {
    y = mu;
  } else {
    y = normal_rng(mu, sigma);
  }
  return y;
}
matrix cholesky_decompose2(matrix A) {
  matrix[rows(A), cols(A)] L;
  int n;
  int nonzero[rows(A)];
  int num_nonzero;
  n = rows(A);
  for (i in 1:n) {
    nonzero[i] = (A[i, i] > 0);
  }
  num_nonzero = sum(nonzero);
  if (num_nonzero == n) {
    L = cholesky_decompose(A);
  } else if (num_nonzero == 0) {
    L = rep_matrix(0.0, n, n);
  } else {
    int idx[num_nonzero];
    vector[n] eps;
    idx = select_indexes(nonzero, num_nonzero);
    L = rep_matrix(0.0, n, n);
    L[idx, idx] = cholesky_decompose(A[idx, idx]);
  }
  return L;
}
vector multi_normal2_rng(vector mu, matrix Sigma) {
  vector[num_elements(mu)] y;
  int n;
  int nonzero[num_elements(mu)];
  int num_nonzero;
  n = num_elements(mu);
  for (i in 1:n) {
    nonzero[i] = (Sigma[i, i] > 0);
  }
  num_nonzero = sum(nonzero);
  if (num_nonzero == n) {
    y = multi_normal_rng(mu, Sigma);
  } else if (num_nonzero == 0) {
    y = mu;
  } else {
    int idx[num_nonzero];
    vector[n] eps;
    idx = select_indexes(nonzero, num_nonzero);
    eps = rep_vector(0.0, n);
    eps[idx] = multi_normal_rng(rep_vector(0.0, num_nonzero), Sigma[idx, idx]);
    y = mu + eps;
  }
  return y;
}
vector multi_normal_cholesky2_rng(vector mu, matrix L) {
  vector[num_elements(mu)] y;
  int n;
  int nonzero[num_elements(mu)];
  int num_nonzero;
  n = num_elements(mu);
  for (i in 1:n) {
    nonzero[i] = (L[i, i] > 0);
  }
  num_nonzero = sum(nonzero);
  if (num_nonzero == n) {
    y = multi_normal_cholesky_rng(mu, L);
  } else if (num_nonzero == 0) {
    y = mu;
  } else {
    int idx[num_nonzero];
    vector[n] eps;
    idx = select_indexes(nonzero, num_nonzero);
    eps = rep_vector(0.0, n);
    eps[idx] = multi_normal_cholesky_rng(rep_vector(0.0, num_nonzero),
                                         L[idx, idx]);
    y = mu + eps;
  }
  return y;
}
vector ssm_update_a(vector a, vector c, matrix T, vector v, matrix K) {
  vector[num_elements(a)] a_new;
  a_new = T * a + K * v + c;
  return a_new;
}
matrix ssm_update_P(matrix P, matrix Z, matrix T,
                           matrix RQR, matrix K) {
  matrix[rows(P), cols(P)] P_new;
  P_new = to_symmetric_matrix(T * P * (T - K * Z)' + RQR);
  return P_new;
}
vector ssm_update_v(vector y, vector a, vector d, matrix Z) {
  vector[num_elements(y)] v;
  v = y - Z * a - d;
  return v;
}
vector ssm_update_v_miss(vector y, vector a, vector d, matrix Z,
                                int p_t, int[] y_idx) {
  vector[num_elements(y)] v;
  int p;
  p = num_elements(y);
  if (p_t < p) {
    v = rep_vector(0., p);
    if (p_t > 0) {
      int idx[p_t];
      vector[p_t] y_star;
      vector[p_t] d_star;
      matrix[p_t, cols(Z)] Z_star;
      idx = y_idx[1:p_t];
      y_star = y[idx];
      d_star = d[idx];
      Z_star = Z[idx, :];
      v[idx] = ssm_update_v(y_star, a, d_star, Z_star);
    }
  } else {
    v = ssm_update_v(y, a, d, Z);
  }
  return v;
}
matrix ssm_update_F(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] F;
  F = to_symmetric_matrix(quad_form(P, Z') + H);
  return F;
}
matrix ssm_update_Finv(matrix P, matrix Z, matrix H) {
  matrix[rows(H), cols(H)] Finv;
  Finv = inverse_spd(to_symmetric_matrix(quad_form(P, Z') + H));
  return Finv;
}
matrix ssm_update_Finv_miss(matrix P, matrix Z, matrix H,
                                   int p_t, int[] y_idx) {
  matrix[rows(H), cols(H)] Finv;
  int p;
  int m;
  p = rows(H);
  m = cols(Z);
  if (p_t < p) {
    Finv = rep_matrix(0., p, p);
    if (p_t > 0) {
      matrix[p_t, m] Z_star;
      matrix[p_t, p_t] H_star;
      int idx[p_t];
      idx = y_idx[1:p_t];
      Z_star = Z[idx, :];
      H_star = H[idx, idx];
      Finv[idx, idx] = ssm_update_Finv(P, Z_star, H_star);
    }
  } else {
    Finv = ssm_update_Finv(P, Z, H);
  }
  return Finv;
}
matrix ssm_update_K(matrix P, matrix Z, matrix T, matrix Finv) {
  matrix[cols(Z), rows(Z)] K;
  K = T * P * Z' * Finv;
  return K;
}
matrix ssm_update_L(matrix Z, matrix T, matrix K) {
  matrix[rows(T), cols(T)] L;
  L = T - K * Z;
  return L;
}
real ssm_update_loglik(vector v, matrix Finv) {
  real ll;
  int p;
  p = num_elements(v);
  ll = (- 0.5 *
        (p * log(2 * pi())
         - log_determinant(Finv)
         + quad_form_sym(Finv, v)
       ));
  return ll;
}
real ssm_update_loglik_miss(vector v, matrix Finv, int p_t, int[] y_idx) {
  real ll;
  int p;
  p = num_elements(v);
  if (p_t == 0) {
    ll = 0.;
  } else if (p_t == p) {
    ll = ssm_update_loglik(v, Finv);
  } else {
    int idx[p_t];
    matrix[p_t, p_t] Finv_star;
    vector[p_t] v_star;
    idx = y_idx[1:p_t];
    Finv_star = Finv[idx, idx];
    v_star = v[idx];
    ll = ssm_update_loglik(v_star, Finv_star);
  }
  return ll;
}
int[,] ssm_filter_idx(int m, int p) {
  int sz[6, 3];
  sz[1, 1] = 1;
  sz[2, 1] = p;
  sz[3, 1] = symmat_size(p);
  sz[4, 1] = m * p;
  sz[5, 1] = m;
  sz[6, 1] = symmat_size(m);
  sz[1, 2] = 1;
  sz[1, 3] = sz[1, 2] + sz[1, 1] - 1;
  for (i in 2:6) {
    sz[i, 2] = sz[i - 1, 3] + 1;
    sz[i, 3] = sz[i, 2] + sz[i, 1] - 1;
  }
  return sz;
}
int ssm_filter_size(int m, int p) {
  int sz;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  sz = idx[6, 3];
  return sz;
}
real ssm_filter_get_loglik(vector x, int m, int p) {
  real y;
  y = x[1];
  return y;
}
vector ssm_filter_get_v(vector x, int m, int p) {
  vector[p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[2, 2], idx[2, 1]);
  return y;
}
matrix ssm_filter_get_Finv(vector x, int m, int p) {
  matrix[p, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = vector_to_symmat(segment(x, idx[3, 2], idx[3, 1]), p);
  return y;
}
matrix ssm_filter_get_K(vector x, int m, int p) {
  matrix[m, p] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = to_matrix_colwise(segment(x, idx[4, 2], idx[4, 1]), m, p);
  return y;
}
vector ssm_filter_get_a(vector x, int m, int p) {
  vector[m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = segment(x, idx[5, 2], idx[5, 1]);
  return y;
}
matrix ssm_filter_get_P(vector x, int m, int p) {
  matrix[m, m] y;
  int idx[6, 3];
  idx = ssm_filter_idx(m, p);
  y = vector_to_symmat(segment(x, idx[6, 2], idx[6, 1]), m);
  return y;
}
vector[] ssm_filter(vector[] y,
                    vector[] d, matrix[] Z, matrix[] H,
                    vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                    vector a1, matrix P1) {
  vector[ssm_filter_size(dims(Z)[3], dims(Z)[2])] res[size(y)];
  int q;
  int n;
  int p;
  int m;
  n = size(y);
  m = dims(Z)[3];
  p = dims(Z)[2];
  q = dims(Q)[2];
  {
    vector[p] d_t;
    matrix[p, m] Z_t;
    matrix[p, p] H_t;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;
    real ll;
    int idx[6, 3];
    idx = ssm_filter_idx(m, p);
    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    RQR = quad_form_sym(Q_t, R_t ');
    a = a1;
    P = P1;
    for (t in 1:n) {
      if (t > 1) {
        if (size(d) > 1) {
          d_t = d[t];
        }
        if (size(Z) > 1) {
          Z_t = Z[t];
        }
        if (size(H) > 1) {
          H_t = H[t];
        }
        if (size(c) > 1) {
          c_t = c[t];
        }
        if (size(T) > 1) {
          T_t = T[t];
        }
        if (size(R) > 1) {
          R_t = R[t];
        }
        if (size(Q) > 1) {
          Q_t = Q[t];
        }
        if (size(R) > 1 || size(Q) > 1) {
          RQR = quad_form_sym(Q_t, R_t ');
        }
      }
      v = ssm_update_v(y[t], a, d_t, Z_t);
      Finv = ssm_update_Finv(P, Z_t, H_t);
      K = ssm_update_K(P, Z_t, T_t, Finv);
      ll = ssm_update_loglik(v, Finv);
      res[t, 1] = ll;
      res[t, idx[2, 2]:idx[2, 3]] = v;
      res[t, idx[3, 2]:idx[3, 3]] = symmat_to_vector(Finv);
      res[t, idx[4, 2]:idx[4, 3]] = to_vector(K);
      res[t, idx[5, 2]:idx[5, 3]] = a;
      res[t, idx[6, 2]:idx[6, 3]] = symmat_to_vector(P);
      if (t < n) {
        a = ssm_update_a(a, c_t, T_t, v, K);
        P = ssm_update_P(P, Z_t, T_t, RQR, K);
      }
    }
  }
  return res;
}
vector[] ssm_filter_miss(vector[] y,
                          vector[] d, matrix[] Z, matrix[] H,
                          vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                          vector a1, matrix P1, int[] p_t, int[,] y_idx) {
  vector[ssm_filter_size(dims(Z)[3], dims(Z)[2])] res[size(y)];
  int q;
  int n;
  int p;
  int m;
  n = size(y);
  m = dims(Z)[3];
  p = dims(Z)[2];
  q = dims(Q)[2];
  {
    vector[p] d_t;
    matrix[p, m] Z_t;
    matrix[p, p] H_t;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;
    real ll;
    int idx[6, 3];
    idx = ssm_filter_idx(m, p);
    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    RQR = quad_form_sym(Q_t, R_t ');
    a = a1;
    P = P1;
    for (t in 1:n) {
      if (t > 1) {
        if (size(d) > 1) {
          d_t = d[t];
        }
        if (size(Z) > 1) {
          Z_t = Z[t];
        }
        if (size(H) > 1) {
          H_t = H[t];
        }
        if (size(c) > 1) {
          c_t = c[t];
        }
        if (size(T) > 1) {
          T_t = T[t];
        }
        if (size(R) > 1) {
          R_t = R[t];
        }
        if (size(Q) > 1) {
          Q_t = Q[t];
        }
        if (size(R) > 1 || size(Q) > 1) {
          RQR = quad_form_sym(Q_t, R_t ');
        }
      }
      v = ssm_update_v_miss(y[t], a, d_t, Z_t, p_t[t], y_idx[t]);
      Finv = ssm_update_Finv_miss(P, Z_t, H_t, p_t[t], y_idx[t]);
      K = ssm_update_K(P, Z_t, T_t, Finv);
      ll = ssm_update_loglik_miss(v, Finv, p_t[t], y_idx[t]);
      res[t, 1] = ll;
      res[t, idx[2, 2]:idx[2, 3]] = v;
      res[t, idx[3, 2]:idx[3, 3]] = symmat_to_vector(Finv);
      res[t, idx[4, 2]:idx[4, 3]] = to_vector(K);
      res[t, idx[5, 2]:idx[5, 3]] = a;
      res[t, idx[6, 2]:idx[6, 3]] = symmat_to_vector(P);
      if (t < n) {
        a = ssm_update_a(a, c_t, T_t, v, K);
        P = ssm_update_P(P, Z_t, T_t, RQR, K);
      }
    }
  }
  return res;
}
real ssm_lpdf(vector[] y,
              vector[] d, matrix[] Z, matrix[] H,
              vector[] c, matrix[] T, matrix[] R, matrix[] Q,
              vector a1, matrix P1) {
  real ll;
  int n;
  int m;
  int p;
  int q;
  n = size(y);
  m = dims(Z)[3];
  p = dims(Z)[2];
  q = dims(Q)[2];
  {
    vector[p] d_t;
    matrix[p, m] Z_t;
    matrix[p, p] H_t;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    vector[n] ll_obs;
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;
    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    RQR = quad_form_sym(Q_t, R_t ');
    a = a1;
    P = P1;
    for (t in 1:n) {
      if (t > 1) {
        if (size(d) > 1) {
          d_t = d[t];
        }
        if (size(Z) > 1) {
          Z_t = Z[t];
        }
        if (size(H) > 1) {
          H_t = H[t];
        }
        if (t < n) {
          if (size(c) > 1) {
            c_t = c[t];
          }
          if (size(T) > 1) {
            T_t = T[t];
          }
          if (size(R) > 1) {
            R_t = R[t];
          }
          if (size(Q) > 1) {
            Q_t = Q[t];
          }
          if (size(R) > 1 || size(Q) > 1) {
            RQR = quad_form_sym(Q_t, R_t ');
          }
        }
      }
      v = ssm_update_v(y[t], a, d_t, Z_t);
      Finv = ssm_update_Finv(P, Z_t, H_t);
      K = ssm_update_K(P, T_t, Z_t, Finv);
      ll_obs[t] = ssm_update_loglik(v, Finv);
      if (t < n) {
        a = ssm_update_a(a, c_t, T_t, v, K);
        P = ssm_update_P(P, Z_t, T_t, RQR, K);
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}
real ssm_miss_lpdf(vector[] y,
                   vector[] d, matrix[] Z, matrix[] H,
                   vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                   vector a1, matrix P1, int[] p_t, int[,] y_idx) {
  real ll;
  int n;
  int m;
  int p;
  int q;
  n = size(y);
  m = dims(Z)[3];
  p = dims(Z)[2];
  q = dims(Q)[2];
  {
    vector[p] d_t;
    matrix[p, m] Z_t;
    matrix[p, p] H_t;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    vector[n] ll_obs;
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;
    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    RQR = quad_form_sym(Q_t, R_t ');
    a = a1;
    P = P1;
    for (t in 1:n) {
      if (t > 1) {
        if (size(d) > 1) {
          d_t = d[t];
        }
        if (size(Z) > 1) {
          Z_t = Z[t];
        }
        if (size(H) > 1) {
          H_t = H[t];
        }
        if (t < n) {
          if (size(c) > 1) {
            c_t = c[t];
          }
          if (size(T) > 1) {
            T_t = T[t];
          }
          if (size(R) > 1) {
            R_t = R[t];
          }
          if (size(Q) > 1) {
            Q_t = Q[t];
          }
          if (size(R) > 1 || size(Q) > 1) {
            RQR = quad_form_sym(Q_t, R_t ');
          }
        }
      }
      v = ssm_update_v_miss(y[t], a, d_t, Z_t, p_t[t], y_idx[t]);
      Finv = ssm_update_Finv_miss(P, Z_t, H_t, p_t[t], y_idx[t]);
      K = ssm_update_K(P, Z_t, T_t, Finv);
      ll_obs[t] = ssm_update_loglik_miss(v, Finv, p_t[t], y_idx[t]);
      if (t < n) {
        a = ssm_update_a(a, c_t, T_t, v, K);
        P = ssm_update_P(P, Z_t, T_t, RQR, K);
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}
real matrix_diff(matrix A, matrix B) {
  real eps;
  real norm_AB;
  real norm_A;
  real a;
  real ab;
  int m;
  int n;
  m = rows(A);
  n = cols(A);
  eps = 0.0;
  norm_A = 0.0;
  norm_AB = 0.0;
  for (i in 1:m) {
    for (j in 1:n) {
      a = fabs(A[i, j]);
      ab = fabs(A[i, j] - B[i, j]);
      if (a > norm_A) {
        norm_A = a;
      }
      if (ab > norm_AB) {
        norm_AB = ab;
      }
    }
  }
  eps = norm_AB / norm_A;
  return eps;
}
real ssm_constant_lpdf(vector[] y,
                        vector d, matrix Z, matrix H,
                        vector c, matrix T, matrix R, matrix Q,
                        vector a1, matrix P1) {
  real ll;
  int n;
  int m;
  int p;
  n = size(y);
  m = cols(Z);
  p = rows(Z);
  {
    vector[n] ll_obs;
    vector[m] a;
    matrix[m, m] P;
    vector[p] v;
    matrix[p, p] Finv;
    matrix[m, p] K;
    matrix[m, m] RQR;
    int converged;
    matrix[m, m] P_old;
    real tol;
    real matdiff;
    converged = 0;
    tol = 1e-7;
    RQR = quad_form_sym(Q, R ');
    a = a1;
    P = P1;
    for (t in 1:n) {
      v = ssm_update_v(y[t], a, d, Z);
      if (converged < 1) {
        Finv = ssm_update_Finv(P, Z, H);
        K = ssm_update_K(P, Z, T, Finv);
      }
      ll_obs[t] = ssm_update_loglik(v, Finv);
      if (t < n) {
        a = ssm_update_a(a, c, T, v, K);
        if (converged < 1) {
          P_old = P;
          P = ssm_update_P(P, Z, T, RQR, K);
          matdiff = matrix_diff(P, P_old);
          if (matdiff < tol) {
            converged = 1;
          }
        }
      }
    }
    ll = sum(ll_obs);
  }
  return ll;
}
vector[] ssm_smooth_states_mean(vector[] filter,
                              matrix[] Z,
                              vector[] c, matrix[] T, matrix[] R, matrix[] Q) {
  vector[dims(Z)[3]] alpha[size(filter)];
  int n;
  int m;
  int p;
  int q;
  n = size(filter);
  m = dims(Z)[3];
  p = dims(Z)[2];
  q = dims(Q)[2];
  {
    vector[m] r[n + 1];
    vector[p] u;
    vector[m] a1;
    matrix[m, m] P1;
    vector[p] v;
    matrix[m, p] K;
    matrix[p, p] Finv;
    matrix[p, m] Z_t;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[m, m] RQR;
    if (size(c) == 1) {
      c_t = c[1];
    }
    if (size(Z) == 1) {
      Z_t = Z[1];
    }
    if (size(T) == 1) {
      T_t = T[1];
    }
    if (size(R) == 1) {
      R_t = R[1];
    }
    if (size(Q) == 1) {
      Q_t = Q[1];
    }
    if (size(Q) == 1 && size(R) == 1) {
      RQR = quad_form_sym(Q[1], R[1]');
    }
    r[n + 1] = rep_vector(0.0, m);
    for (s in 0:(n - 1)) {
      int t;
      t = n - s;
      if (size(Z) > 1) {
        Z_t = Z[t];
      }
      if (size(T) > 1) {
        T_t = T[t];
      }
      K = ssm_filter_get_K(filter[t], m, p);
      v = ssm_filter_get_v(filter[t], m, p);
      Finv = ssm_filter_get_Finv(filter[t], m, p);
      u = Finv * v - K ' * r[t + 1];
      r[t] = Z_t ' * u + T_t ' * r[t + 1];
    }
    a1 = ssm_filter_get_a(filter[1], m, p);
    P1 = ssm_filter_get_P(filter[1], m, p);
    alpha[1] = a1 + P1 * r[1];
    for (t in 1:(n - 1)) {
      if (size(c) > 1) {
        c_t = c[t];
      }
      if (size(T) > 1) {
        T_t = T[t];
      }
      if (size(Q) > 1) {
        Q_t = Q[t];
      }
      if (size(R) > 1) {
        R_t = R[t];
      }
      if (size(Q) > 1 || size(R) > 1) {
        RQR = quad_form_sym(Q_t, R_t');
      }
      alpha[t + 1] = c_t + T_t * alpha[t] + RQR * r[t + 1];
    }
  }
  return alpha;
}
int[,] ssm_sim_idx(int m, int p, int q) {
  int sz[2, 3];
  sz[1, 1] = p;
  sz[2, 1] = m;
  sz[1, 2] = 1;
  sz[1, 3] = sz[1, 2] + sz[1, 1] - 1;
  sz[2, 2] = sz[2 - 1, 3] + 1;
  sz[2, 3] = sz[2, 2] + sz[2, 1] - 1;
  return sz;
}
int ssm_sim_size(int m, int p, int q) {
  int sz;
  sz = ssm_sim_idx(m, p, q)[2, 3];
  return sz;
}
vector ssm_sim_get_y(vector x, int m, int p, int q) {
  vector[p] y;
  int idx[2, 3];
  idx = ssm_sim_idx(m, p, q);
  y = x[idx[1, 2]:idx[1, 3]];
  return y;
}
vector ssm_sim_get_a(vector x, int m, int p, int q) {
  vector[m] a;
  int idx[2, 3];
  idx = ssm_sim_idx(m, p, q);
  a = x[idx[2, 2]:idx[2, 3]];
  return a;
}
vector[] ssm_sim_rng(int n,
                    vector[] d, matrix[] Z, matrix[] H,
                    vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                    vector a1, matrix P1) {
  vector[ssm_sim_size(dims(Z)[3], dims(Z)[2], dims(Q)[2])] ret[n];
  int p;
  int m;
  int q;
  m = dims(Z)[3];
  p = dims(Z)[2];
  q = dims(Q)[2];
  {
    vector[p] d_t;
    matrix[p, m] Z_t;
    matrix[p, p] H_t;
    matrix[p, p] HL;
    vector[m] c_t;
    matrix[m, m] T_t;
    matrix[m, q] R_t;
    matrix[q, q] Q_t;
    matrix[q, q] QL;
    vector[p] y;
    vector[p] eps;
    vector[m] a;
    vector[q] eta;
    vector[p] zero_p;
    vector[q] zero_q;
    vector[m] zero_m;
    int idx[2, 3];
    d_t = d[1];
    Z_t = Z[1];
    H_t = H[1];
    HL = cholesky_decompose2(H_t);
    c_t = c[1];
    T_t = T[1];
    R_t = R[1];
    Q_t = Q[1];
    QL = cholesky_decompose2(Q_t);
    idx = ssm_sim_idx(m, p, q);
    zero_p = rep_vector(0.0, p);
    zero_q = rep_vector(0.0, q);
    zero_m = rep_vector(0.0, m);
    a = multi_normal2_rng(a1, P1);
    for (t in 1:n) {
      ret[t, idx[2, 2]:idx[2, 3]] = a;
      if (t > 1) {
        if (size(d) > 1) {
          d_t = d[t];
        }
        if (size(Z) > 1) {
          Z_t = Z[t];
        }
        if (size(H) > 1) {
          H_t = H[t];
          HL = cholesky_decompose2(H_t);
        }
      }
      eps = multi_normal_cholesky2_rng(zero_p, HL);
      y = d_t + Z_t * a + eps;
      ret[t, idx[1, 2]:idx[1, 3]] = y;
      if (t < n) {
        if (size(c) > 1) {
          c_t = c[t];
        }
        if (size(T) > 1) {
          T_t = T[t];
        }
        if (size(R) > 1) {
          R_t = R[t];
        }
        if (size(Q) > 1) {
          Q_t = Q[t];
          QL = cholesky_decompose2(Q_t);
        }
        eta = multi_normal_cholesky2_rng(zero_q, QL);
        a = c_t + T_t * a + R_t * eta;
      }
    }
  }
  return ret;
}
vector[] ssm_simsmo_states_rng(vector[] filter,
                               vector[] d, matrix[] Z, matrix[] H,
                               vector[] c, matrix[] T, matrix[] R, matrix[] Q,
                               vector a1, matrix P1) {
    vector[dims(Z)[3]] draws[size(filter)];
    int n;
    int p;
    int m;
    int q;
    n = size(filter);
    m = dims(Z)[3];
    p = dims(Z)[2];
    q = dims(Q)[2];
    {
      vector[ssm_filter_size(m, p)] filter_plus[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[p] y[n];
      vector[m] alpha_hat_plus[n];
      vector[m] alpha_hat[n];
      alpha_hat = ssm_smooth_states_mean(filter, Z, c, T, R, Q);
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      filter_plus = ssm_filter(y, d, Z, H, c, T, R, Q, a1, P1);
      alpha_hat_plus = ssm_smooth_states_mean(filter_plus, Z, c, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_a(sims[i], m, p, q)
                    - alpha_hat_plus[i]
                    + alpha_hat[i]);
      }
    }
    return draws;
}
vector[] ssm_simsmo_states_miss_rng(vector[] filter,
                                    vector[] d, matrix[] Z, matrix[] H,
                                    vector[] c, matrix[] T,
                                    matrix[] R, matrix[] Q,
                                    vector a1, matrix P1,
                                    int[] p_t, int[,] y_idx) {
    vector[dims(Z)[3]] draws[size(filter)];
    int n;
    int p;
    int m;
    int q;
    n = size(filter);
    m = dims(Z)[3];
    p = dims(Z)[2];
    q = dims(Q)[2];
    {
      vector[ssm_filter_size(m, p)] filter_plus[n];
      vector[ssm_sim_size(m, p, q)] sims[n];
      vector[p] y[n];
      vector[m] alpha_hat_plus[n];
      vector[m] alpha_hat[n];
      alpha_hat = ssm_smooth_states_mean(filter, Z, c, T, R, Q);
      sims = ssm_sim_rng(n, d, Z, H, c, T, R, Q, a1, P1);
      for (i in 1:n) {
        y[i] = ssm_sim_get_y(sims[i], m, p, q);
      }
      filter_plus = ssm_filter_miss(y, d, Z, H, c, T, R, Q,
                                    a1, P1, p_t, y_idx);
      alpha_hat_plus = ssm_smooth_states_mean(filter_plus, Z, c, T, R, Q);
      for (i in 1:n) {
        draws[i] = (ssm_sim_get_a(sims[i], m, p, q)
                    - alpha_hat_plus[i]
                    + alpha_hat[i]);
      }
    }
    return draws;
}
vector pacf_to_acf(vector x) {
  vector[num_elements(x)] x_new;
  vector[num_elements(x)] work;
  real a;
  int p;
  p = num_elements(x);
  work = x;
  x_new = x;
  if (p > 1) {
    for (j in 2:p) {
      a = x_new[j];
      for (k in 1:(j - 1)) {
        work[k] = work[k] - a * x_new[j - k];
      }
      for (k in 1:j) {
        x_new[k] = work[k];
      }
    }
  }
  return x_new;
}
vector constrain_stationary(vector x) {
  vector[num_elements(x)] r;
  int n;
  n = num_elements(x);
  for (i in 1:n) {
    r[i] = tanh(x[i]);
  }
  return pacf_to_acf(r);
}
vector acf_to_pacf(vector x) {
  vector[num_elements(x)] x_new;
  vector[num_elements(x)] work;
  real a;
  int p;
  p = num_elements(x);
  work = x;
  x_new = x;
  if (p > 1) {
    for(i in 0:(p - 2)) {
      int j;
      j = p - i;
      a = x_new[j];
      for(k in 1:(j - 1)) {
        work[k]  = (x_new[k] + a * x_new[j - k]) / (1 - pow(a, 2));
      }
      for (k in 1:j) {
        x_new[k] = work[k];
      }
    }
  }
  return x_new;
}
vector unconstrain_stationary(vector x) {
  matrix[num_elements(x), num_elements(x)] y;
  vector[num_elements(x)] r;
  vector[num_elements(x)] z;
  int n;
  n = num_elements(x);
  r = acf_to_pacf(x);
  for (i in 1:n) {
    z[i] = atanh(r[i]);
  }
  return z;
}
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
matrix stationary_cov(matrix T, matrix RQR) {
  matrix[rows(T), cols(T)] P;
  int m;
  m = rows(T);
  if (m == 1) {
    P[1, 1] = RQR[1, 1] / (1.0 - pow(T[1, 1], 2));
  } else {
    matrix[rows(T) * rows(T), rows(T) * rows(T)] TT;
    vector[rows(T) * rows(T)] RQR_vec;
    int m2;
    m2 = m * m;
    RQR_vec = to_vector(RQR);
    TT = - kronecker_prod(T, T);
    for (i in 1:m2) {
      TT[i, i] = 1.0 + TT[i, i];
    }
    P = to_matrix_colwise(inverse(TT) * RQR_vec, m, m);
  }
  return P;
}

