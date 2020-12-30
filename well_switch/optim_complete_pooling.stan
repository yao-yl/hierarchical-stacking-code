data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N,K] lpd_point;
  vector[K] lambda;
}

transformed data{
  matrix[N,K] exp_lpd_point; 
  
  exp_lpd_point=exp(lpd_point);
}

parameters {
   simplex[K] w;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
transformed parameters{
  vector[K] w_vec;
  w_vec=w;
}

model {
  for (i in 1: N) {
    target += log( exp_lpd_point[i,] * w_vec );
  }
  w~dirichlet(lambda);
}


