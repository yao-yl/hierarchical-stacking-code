functions {
  vector gp_pred_rng(real[] x2,
                     vector y1, real[] x1,
                     real alpha, real rho,  real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   cov_exp_quad(x1, alpha, rho) + diag_matrix(rep_vector(delta, N1));
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
  int<lower=0> N;
  real  x[N];
  matrix[N,2] lpd_point;
  
  int<lower=1> N_test;
  real x_test[N_test];
  real<lower=0> delta;
  real gamma_a;
  real gamma_b;
}

transformed data{
  matrix[N,2] exp_lpd_point=exp(lpd_point);
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  vector[N] eta;
}


transformed parameters{
  vector[N] w_vec;
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(x, alpha, rho);
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;
    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }
  w_vec=inv_logit(f);
}

model {
  rho ~ inv_gamma(gamma_a, gamma_b);
  alpha ~ std_normal();
  eta ~ std_normal();
  for (i in 1: N) {
    target += log(exp_lpd_point[i,1] * w_vec[i]+ exp_lpd_point[i,2]* (1-w_vec[i]));
  }
}
generated quantities{
  vector[N_test] f_predict = gp_pred_rng(x_test, f, x, alpha, rho,  delta);
  vector[N_test] w_predict=inv_logit(f_predict);
}


