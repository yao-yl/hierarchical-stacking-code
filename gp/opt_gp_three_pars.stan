data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
}


parameters {
  real log_rho;
  real log_alpha;
  real log_sigma;
}

transformed parameters{
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  
  rho=exp(log_rho);
  alpha=exp(log_alpha);
  sigma=exp(log_sigma);
}


model {
  matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)
  + diag_matrix(rep_vector(square(sigma), N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
  y ~ multi_normal_cholesky(rep_vector(0, N), L_cov);
  rho ~ cauchy(0,36);
  alpha ~ cauchy(0,36);
  sigma~ cauchy(0,36);
}

