data {
	int<lower=1> N;  // number of data points
  int<lower=1> N_test;  // number of data points
	int<lower=1> d; //number of input variables
	int<lower=1> d_discrete; // number of discrete dummy inputs
	int<lower=3> K;  // number of models  when K=2, replace softmax by logistic for higher efficiency 
	matrix[N,d] X;   // predictors (including continous and discrete in dummy variables, no constant)
	matrix[N_test,d] X_test;             
	matrix[N,K] lpd_point;
	real tau_mu;
	real<lower=0> tau_sigma;
}

transformed data{
	matrix[N,K] exp_lpd_point=exp(lpd_point);
}

parameters {
	vector[K-1] mu;
	vector<lower=0>[K-1] sigma;
	vector[d-d_discrete] beta_con[K-1];
	vector[d_discrete] tau[K-1];
}


transformed parameters{
	vector[d] beta[K-1];
	simplex[K] w[N];
 matrix[N,K] f;
	  for(k in 1:(K-1))
	  {
	  	beta[k]= append_row(mu[k]+ sigma[k]*tau[k], beta_con[k]);
	  }
		for(k in 1:(K-1))
			f[,k]= X * beta[k];
		f[,K]=rep_vector(0, N);
    for (n in 1:N)
		  w[n]=softmax( to_vector(f[n, 1:K])  );
	}
model {
	for(k in 1:(K-1)){
  	tau[k]~std_normal();
  	beta_con[k]~normal(0,1);
	}
  	mu~normal(0,tau_mu);
  	sigma~normal(0,tau_sigma);
	for (i in 1: N) 
		target += log( exp_lpd_point[i,] * w[i] );
}
generated quantities{
	matrix[N_test,K] f_test;
		for(k in 1:(K-1))
			f_test[,k]= X_test * beta[k];
		f_test[,K]=rep_vector(0, N_test);
}
 
