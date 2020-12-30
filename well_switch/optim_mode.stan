// Stan code for h-stacking, this time we use optimizing to derive mode. To this end we remove all hyper-priors (otherwise the joint mode will be attained at sigma=0) 
data {
	int<lower=1> N;  // number of data points
	int<lower=1> N_test;  // number of data points
	int<lower=1> d; //number of input variables
	int<lower=1> d_discrete; // number of discrete dummy inputs
	int<lower=3> K;  // number of models  when K=2, replace softmax by logistic for higher efficiency 
	matrix[N,d] X;   // predictors (including continous and discrete in dummy variables, no constant)
	matrix[N_test,d] X_test;             
	matrix[N,K] lpd_point;
}

transformed data{
	matrix[N,K] exp_lpd_point=exp(lpd_point);
}

parameters {
	vector[d] beta[K-1];
}


transformed parameters{
	simplex[K] w[N];
	matrix[N,K] f;
	for(k in 1:(K-1))
		f[,k]= X * beta[k];
	f[,K]=rep_vector(0, N);
	for (n in 1:N)
		w[n]=softmax(to_vector(f[n, 1:K]));
}
model {
	for (i in 1: N) 
		target += log( exp_lpd_point[i,] * w[i] );
}
generated quantities{
	matrix[N_test,K] f_test;
	for(k in 1:(K-1))
		f_test[,k]= X_test * beta[k];
	f_test[,K]=rep_vector(0, N_test);
}
