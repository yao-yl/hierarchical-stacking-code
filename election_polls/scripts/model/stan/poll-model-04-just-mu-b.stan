data{
  int N_national_polls;    // Number of polls
  int N_state_polls;    // Number of polls
  int T;    // Number of days
  int S;    // Number of states (for which at least 1 poll is available) + 1
  int P;    // Number of pollsters
  int M;    // Number of poll modes
  int Pop;    // Number of poll populations
  int<lower = 1, upper = S + 1> state[N_state_polls]; // State index
  int<lower = 1, upper = T> day_state[N_state_polls];   // Day index
  int<lower = 1, upper = T> day_national[N_national_polls];   // Day index
  int<lower = 1, upper = P> poll_state[N_state_polls];  // Pollster index
  int<lower = 1, upper = P> poll_national[N_national_polls];  // Pollster index
  int<lower = 1, upper = M> poll_mode_state[N_state_polls];  // Poll mode index
  int<lower = 1, upper = M> poll_mode_national[N_national_polls];  // Poll mode index
  int<lower = 1, upper = Pop> poll_pop_state[N_state_polls];  // Poll mode index
  int<lower = 1, upper = Pop> poll_pop_national[N_national_polls];  // Poll mode index
  int n_democrat_national[N_national_polls];
  int n_two_share_national[N_national_polls];
  int n_democrat_state[N_state_polls];
  int n_two_share_state[N_state_polls];
  vector<lower = 0, upper = 1.0>[N_national_polls] unadjusted_national;
  vector<lower = 0, upper = 1.0>[N_state_polls] unadjusted_state;
  // cov_matrix[S] ss_cov_mu_b_walk;
  // cov_matrix[S] ss_cov_mu_b_T;
  // cov_matrix[S] ss_cov_poll_bias;
  //*** prior input
  vector[S] mu_b_prior;
  vector[S] state_weights;
  real sigma_c;
  real sigma_m;
  real sigma_pop;
  real sigma_measure_noise_national;
  real sigma_measure_noise_state;
  real sigma_e_bias;
  // covariance matrix and scales
  cov_matrix[S] state_covariance_0;
  real random_walk_scale;
  real mu_b_T_scale;
  real polling_bias_scale;


  // For evaluation
  int N_national_polls_test;    // Number of polls
  int N_state_polls_test;    // Number of polls
  int<lower = 1, upper = S + 1> state_test[N_state_polls_test]; // State index
  int<lower = 1, upper = T> day_state_test[N_state_polls_test];   // Day index
  int<lower = 1, upper = T> day_national_test[N_national_polls_test];   // Day index
  int<lower = 1, upper = P> poll_state_test[N_state_polls_test];  // Pollster index
  int<lower = 1, upper = P> poll_national_test[N_national_polls_test];  // Pollster index
  int<lower = 1, upper = M> poll_mode_state_test[N_state_polls_test];  // Poll mode index
  int<lower = 1, upper = M> poll_mode_national_test[N_national_polls_test];  // Poll mode index
  int<lower = 1, upper = Pop> poll_pop_state_test[N_state_polls_test];  // Poll mode index
  int<lower = 1, upper = Pop> poll_pop_national_test[N_national_polls_test];  // Poll mode index
  int n_democrat_national_test[N_national_polls_test];
  int n_two_share_national_test[N_national_polls_test];
  int n_democrat_state_test[N_state_polls_test];
  int n_two_share_state_test[N_state_polls_test];
  vector<lower = 0, upper = 1.0>[N_national_polls_test] unadjusted_national_test;
  vector<lower = 0, upper = 1.0>[N_state_polls_test] unadjusted_state_test;
}
transformed data {
  real national_cov_matrix_error_sd = sqrt(transpose(state_weights) * state_covariance_0 * state_weights);
  cholesky_factor_cov[S] cholesky_ss_cov_poll_bias;
  cholesky_factor_cov[S] cholesky_ss_cov_mu_b_T;
  cholesky_factor_cov[S] cholesky_ss_cov_mu_b_walk;
  // scale covariance
  matrix[S, S] ss_cov_poll_bias = state_covariance_0 * square(polling_bias_scale/national_cov_matrix_error_sd);
  matrix[S, S] ss_cov_mu_b_T = state_covariance_0 * square(mu_b_T_scale/national_cov_matrix_error_sd);
  matrix[S, S] ss_cov_mu_b_walk = state_covariance_0 * square(random_walk_scale/national_cov_matrix_error_sd);
  // transformation
  cholesky_ss_cov_poll_bias = cholesky_decompose(ss_cov_poll_bias);
  cholesky_ss_cov_mu_b_T = cholesky_decompose(ss_cov_mu_b_T);
  cholesky_ss_cov_mu_b_walk = cholesky_decompose(ss_cov_mu_b_walk);
}
parameters {
  vector[S] raw_mu_b_T;
  matrix[S, T] raw_mu_b;
  
  vector[N_national_polls_test] raw_measure_noise_national_test;
  vector[N_state_polls_test] raw_measure_noise_state_test;
}
transformed parameters {
  //*** parameters
  matrix[S, T] mu_b;

  vector[T] national_mu_b_average;

  //*** containers
  vector[N_state_polls] logit_pi_democrat_state;
  vector[N_national_polls] logit_pi_democrat_national;
  //*** construct parameters
  mu_b[:,T] = cholesky_ss_cov_mu_b_T * raw_mu_b_T + mu_b_prior;  // * mu_b_T_model_estimation_error
  for (i in 1:(T-1)) mu_b[:, T - i] = cholesky_ss_cov_mu_b_walk * raw_mu_b[:, T - i] + mu_b[:, T + 1 - i];
  national_mu_b_average = transpose(mu_b) * state_weights;

  //*** fill pi_democrat
  for (i in 1:N_state_polls){
    logit_pi_democrat_state[i] =
      mu_b[state[i], day_state[i]];
  }
  logit_pi_democrat_national =
    national_mu_b_average[day_national];
}

model {
  //*** priors
  raw_measure_noise_national_test ~ std_normal();
  raw_measure_noise_state_test ~ std_normal();
  
  raw_mu_b_T ~ std_normal();
  to_vector(raw_mu_b) ~ std_normal();

  //*** likelihood
  n_democrat_state ~ binomial_logit(n_two_share_state, logit_pi_democrat_state);
  n_democrat_national ~ binomial_logit(n_two_share_national, logit_pi_democrat_national);
}

generated quantities {
  matrix[T, S] predicted_score;
  vector[N_national_polls_test] predicted_national_polls;
  vector[N_state_polls_test] predicted_state_polls;
  vector[N_national_polls_test] lpd_national_polls;
  vector[N_state_polls_test] lpd_state_polls;
  for (s in 1:S){
    predicted_score[1:T, s] = inv_logit(to_vector(mu_b[s, 1:T]));
  }
  for (n in 1:N_state_polls_test) {
    real logit_pi_tmp;
    logit_pi_tmp = mu_b[state_test[n], day_state_test[n]];
    predicted_state_polls[n] = inv_logit(logit_pi_tmp +
      raw_measure_noise_state_test[n] * sigma_measure_noise_state);
    lpd_state_polls[n]       = binomial_lpmf(n_democrat_state_test[n] |
                                                  n_two_share_state_test[n],
                                                  inv_logit(logit_pi_tmp));
  }
  for (n in 1:N_national_polls_test) {
    real logit_pi_tmp;
    logit_pi_tmp = national_mu_b_average[day_national_test[n]] +
      raw_measure_noise_national_test[n] * sigma_measure_noise_national;
    predicted_national_polls[n] = inv_logit(logit_pi_tmp);
    lpd_national_polls[n]       = binomial_lpmf(n_democrat_national_test[n] |
                                                  n_two_share_national_test[n],
                                                  inv_logit(logit_pi_tmp));
  }
}

