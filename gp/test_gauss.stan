functions {
  vector gp_pred_rng(real[] x2,
                     vector y1, real[] x1,
                     real alpha, real rho, real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   cov_exp_quad(x1, alpha, rho)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =   cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}

data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
  
  int<lower=1> N_test;
  real x_test[N_test];
  real y_test[N_test];

  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

transformed data {
  matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)
  + diag_matrix(rep_vector(1e-10, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
}

parameters {}
model {}

generated quantities {
  vector[N_test] f_predict = gp_pred_rng(x_test, y, x, alpha, rho, sigma, 1e-10);
  vector[N_test] log_lik_test;
  vector[N] log_lik;
  vector[N] f;
  {
    vector[N] zeros=to_vector(rep_array(0.0,N));
    vector[N] ones=to_vector(rep_array(1.0, N));
    vector[N] eta = to_vector( normal_rng(zeros,ones));
    f = L_cov * to_vector(eta);
  }
  for (n in 1:N)
    log_lik[n] = normal_lpdf (y[n]|f[n], sigma);
  for (n in 1:N_test)
    log_lik_test[n] = normal_lpdf(y_test[n]|f_predict[n], sigma);
}
