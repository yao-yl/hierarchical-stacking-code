data {
  int<lower=0> N;             // number of data points
    int<lower=0> N_test;             // number of data points
  int<lower=0> P;             // number of predictors (including intercept)
  matrix[N,P] X;              // predictors (including 1s for intercept)
  int<lower=0,upper=1> y[N];  // binary outcome
  matrix[N_test,P] X_test;              // predictors (including 1s for intercept)
  int<lower=0,upper=1> y_test[N_test];  // binary outcome
}
parameters {
  vector[P] beta;
}
model {
  beta ~ normal(0, 3);
  y ~ bernoulli_logit(X * beta);
}
generated quantities {
  vector[N] log_lik;
  vector[N_test] log_lik_test;
  for (n in 1:N) 
    log_lik[n] = bernoulli_logit_lpmf(y[n] | X[n] * beta);
  for (n in 1:N_test) 
    log_lik_test[n] = bernoulli_logit_lpmf(y_test[n] | X_test[n] * beta);
}
