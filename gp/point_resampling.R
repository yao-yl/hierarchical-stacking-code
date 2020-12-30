#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## stacking in GP for mode-based approximation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
setwd("~/Desktop/hstacking/code/gp")
source("stacking_utlitity.R")
## data from Neal
data.neal=read.table("odata.txt")
# odata.txt: Data from Neal's Software 
# http://llrforge.in2p3.fr/svn/cms/cms/trunk/Software/fbm.2004-11-10/doc/manual
colnames(data.neal)=c("x","y")

## The first half: training
x=data.neal$x [1:100]
y=data.neal$y [1:100]
# We use the second half of training data as  validation to avoid LOO in mode-searching
x_val=data.neal$x [101:200]
y_val=data.neal$y [101:200]
# independent test data 
set.seed(100)
n_test=200
x_test=sort(rnorm(n_test, 0, 1))
f_test=0.3 + 0.4*x_test + 0.5 * sin(2.7*x_test) + 1.1/(1 + x_test^2 )
id.out=sort(sample(1:n_test, 0.05*n_test))
y_test=f_test+rnorm(n_test, 0,0.1)
y_test[id.out]=f_test[id.out]+rnorm(length(id.out), 0,1)
plot(x_test, y_test)
x_test=c(x_val, x_test) 
# For code simplicty combine these two vector,
# but the actual test data x_test[201:400] shall not be revealed until the evaluation phase
y_test=c(y_val, y_test)
n_test=300
n=length(x)
f=0.3 + 0.4*x + 0.5 * sin(2.7*x) + 1.1/(1 + x^2 )
m_three_pars=stan_model("opt_gp_three_pars.stan")
fit_three_pars <- sampling(m_three_pars, data=list(N=length(x), x=x, y=y), chains = 1, iter=10, seed=5838298)
#  you should find two modes:
opt_1= optimizing(m_three_pars, 
									data=list(N=n, x=x, y=y),
									init=list(log_rho= (1), 
														log_alpha=0.7, 
														log_sigma=0.1),
									hessian=T,
									seed=2020 ) 
opt_1_par=opt_1$par
print(opt_1_par)
##
##rho     alpha     sigma 
##0.9884577 1.8481367 0.2576907 

opt_2= optimizing(m_three_pars, 
									data=list(N=n, x=x, y=y),
									init=list(log_rho= -1, 
														log_alpha=-5, 
														log_sigma=log(5)),
									hessian=T,
									seed=2020)
opt_2_par=opt_2$par
print(opt_2_par)
##
# rho     alpha     sigma 
# 0.4766149 1.1859285 0.2358455 

### Laplace approx ==================================
grid_normal=function( n=2000){  # generate 3-D N(0,I)
	temp=matrix(NA,  3, n)
	for(i in 1:3)
		temp[i,]=rnorm(n,0,1)
	return(temp)
}
set.seed(2020)
row_grid_nomral=grid_normal()
find_grid=function(row_grid=row_grid_nomral,  opt){
	x1=opt$par[1:3]
	S=dim(row_grid_nomral)[2]
	shift=rbind( rep(x1[1],S), rep(x1[2],S) , rep(x1[3],S) )
	h1=opt$hessian
	cov1=solve(-h1)	
	eig_d= eigen(cov1)
	grid= eig_d$vectors %*% diag( sqrt(eig_d$values))    %*%    row_grid
	grid=t(grid)+t(shift)
	rownames(grid)=c()
	colnames(grid)=names(x1)
	return(grid)
}
grid_Laplace_1=find_grid(row_grid=row_grid_nomral,  opt=opt_1)
grid_Laplace_2=find_grid(row_grid=row_grid_nomral,  opt=opt_2)
S=dim(grid_Laplace_1)[1]
grid_lp1=grid_lp2=c()
for(i in 1: S)
{
	grid_lp1[i]=log_prob(fit_three_pars, upars=grid_Laplace_1[i,]) 
	grid_lp2[i]=log_prob(fit_three_pars, upars=grid_Laplace_2[i,]) 
}
exp(max(grid_lp1)-max(grid_lp2))
m_loo=stan_model("loo.stan")
laplace_gqs1=gqs(m_loo, data=list(N=length(x), x=x, y=y, 
																	N_test=length(x_test), 
																	x_test=x_test, y_test=y_test),
								 draws= expand_draws(grid_Laplace_1) )

laplace_gqs2=gqs(m_loo, data=list(N=length(x), x=x, y=y, 
																	N_test=length(x_test), 
																	x_test=x_test, y_test=y_test),
								 draws= expand_draws(grid_Laplace_2) )


fit_draw_lik1=extract(laplace_gqs1, pars="log_lik")$log_lik
fit_draw_lik2=extract(laplace_gqs2, pars="log_lik")$log_lik
fit_draw_test1=extract(laplace_gqs1, pars="log_lik_test")$log_lik_test
fit_draw_test2=extract(laplace_gqs2, pars="log_lik_test")$log_lik_test

fit_draw_lik=array(NA, c(dim(fit_draw_lik1)[1],  2,    dim(fit_draw_lik1)[2]))
fit_draw_lik[,1,]=fit_draw_lik1
fit_draw_lik[,2,]=fit_draw_lik2



fit_draw_lik_test=array( NA, c(dim(fit_draw_test1)[1],  2,  dim(fit_draw_test1)[2]))
fit_draw_lik_test[,1,]=fit_draw_test1
fit_draw_lik_test[,2,]=fit_draw_test2

# find stacking weights use the validation set
log_lik_test_agg=test_aggregrate(fit_draw_lik_test)
validata_agg=log_lik_test_agg[1:100,]   
st_weight=stacking_weights(validata_agg)


hierstack=hier_stacking(x=x_test[1:100], x_test =x_test[101:300], plpd= validata_agg)
w_test_mat=cbind(colMeans(hierstack$w_test), 1-colMeans(hierstack$w_test))
test_lpd_hier=log_score_loo_pointwise( w_mat=w_test_mat, lpd_point= log_lik_test_agg[101:300,])


w_mode=c(exp(max(grid_lp1)-max(grid_lp2)), 1)
w_mode=w_mode/sum(w_mode)

w_is=c(exp(mean(grid_lp1)-mean(grid_lp2)), 1)
w_is=w_is/sum(w_is)



log_score_test_laplace=c(test_lpd_hier,
	log_score_loo(w=st_weight, lpd_point= log_lik_test_agg[101:300,]),
													log_score_loo(w=w_mode, lpd_point= log_lik_test_agg[101:300,]),
													log_score_loo(w=w_is, lpd_point= log_lik_test_agg[101:300,]))

# grid_approx and inportance resampling==================================
grid_uniform=function(lower=-3.2, upper=3.2, n_0=30 ){
	temp=matrix(NA,  3, n_0^3)
	xx=seq(lower,upper,length.out = n_0)
	for (i in c(1:n_0))
		for (j in c(1:n_0))
			for (l in c(1:n_0))
			   temp[,(i-1)*n_0^2+(j-1)*n_0+l]= c(xx[i], xx[j], xx[l])
	return(temp)
}
row_grid=grid_uniform()

find_grid=function(row_grid=row_grid,  opt){
	x1=opt$par[1:3]
	S=dim(row_grid)[2]
	shift=rbind( rep(x1[1],S), rep(x1[2],S) , rep(x1[3],S) )
	h1=opt$hessian
	cov1=solve(-h1)	
	eig_d= eigen(cov1)
	grid= eig_d$vectors %*% diag( sqrt(eig_d$values))    %*%    row_grid
	grid=t(grid)+t(shift)
	rownames(grid)=c()
	colnames(grid)=names(x1)
	return(grid)
}


grid_unif_1=find_grid(row_grid=row_grid,  opt=opt_1)
grid_unif_2=find_grid(row_grid=row_grid,  opt=opt_2)
S=dim(grid_unif_1)[1]
grid_lp1_unif=grid_lp2_unif=c()
for(i in 1: S)
{
	grid_lp1_unif[i]=log_prob(fit_three_pars, upars=grid_unif_1[i,]) 
	grid_lp2_unif[i]=log_prob(fit_three_pars, upars=grid_unif_2[i,]) 
}

threshold=log(max(grid_lp1_unif, grid_lp2_unif))-log(10) #at least 0.1 * mode height
id.remain1= sample(1:S, 1000,  prob =  exp(grid_lp1_unif-max(grid_lp1_unif))* c(grid_lp1_unif> threshold)  )
id.remain2= sample(1:S, 1000,  prob =  exp(grid_lp2_unif-max(grid_lp2_unif))* c(grid_lp2_unif> threshold)   )

grid_unif_resample1= grid_unif_1[id.remain1,]
grid_unif_resample2=grid_unif_2[id.remain2,]

expand_draws=function(grid)
{
	expand_draw_matrix=grid
	temp= dimnames (as.matrix(fit_three_pars))
	temp[[2]]=temp[[2]][1:3]
	dimnames(expand_draw_matrix) = temp
	return(expand_draw_matrix)
}
unif_gqs1=gqs(m_loo, data=list(N=length(x), x=x, y=y, 
																	N_test=length(x_test), 
																	x_test=x_test, y_test=y_test),
								 draws= expand_draws(grid=grid_unif_resample1) )

unif_gqs2=gqs(m_loo, data=list(N=length(x), x=x, y=y, 
																	N_test=length(x_test), 
																	x_test=x_test, y_test=y_test),
								 draws= expand_draws(grid=grid_unif_resample2) )
fit_draw_lik1=extract(unif_gqs1, pars="log_lik")$log_lik
fit_draw_lik2=extract(unif_gqs2, pars="log_lik")$log_lik
fit_draw_test1=extract(unif_gqs1, pars="log_lik_test")$log_lik_test
fit_draw_test2=extract(unif_gqs2, pars="log_lik_test")$log_lik_test
fit_draw_lik=array(NA, c(dim(fit_draw_lik1)[1],  2,    dim(fit_draw_lik1)[2]))
fit_draw_lik[,1,]=fit_draw_lik1
fit_draw_lik[,2,]=fit_draw_lik2
fit_draw_lik_test=array( NA, c(dim(fit_draw_test1)[1],  2,  dim(fit_draw_test1)[2]))
fit_draw_lik_test[,1,]=fit_draw_test1
fit_draw_lik_test[,2,]=fit_draw_test2

log_lik_test_agg=test_aggregrate(fit_draw_lik_test)
validata_agg=log_lik_test_agg[1:100,]
st_weight=stacking_weights(validata_agg)
w_mode=c(exp(max(grid_lp1)-max(grid_lp2)), 1)
w_mode=w_mode/sum(w_mode)
w_is=c(exp(mean(grid_lp1)-mean(grid_lp2)), 1)
w_is=w_is/sum(w_is)

hierstack=hier_stacking(x=x_test[1:100], x_test =x_test[101:300], plpd= validata_agg)
w_test_mat=cbind(colMeans(hierstack$w_test), 1-colMeans(hierstack$w_test))
test_lpd_hier=log_score_loo_pointwise( w_mat=w_test_mat, lpd_point= log_lik_test_agg[101:300,])

log_score_test_unif=c( test_lpd_hier, log_score_loo(w=st_weight, lpd_point= log_lik_test_agg[101:300,]),
											 log_score_loo(w=w_mode, lpd_point= log_lik_test_agg[101:300,]),
											 log_score_loo(w=w_is, lpd_point= log_lik_test_agg[101:300,]))
### MAP-II =======================================================
test_stan=stan_model('test_gauss.stan')
mode1_data <- list(alpha= opt_1_par[5], rho= opt_1_par[4], sigma= opt_1_par[6], N=length(x), x=x, y=y,
									 N_test=length(x_test), x_test=x_test, y_test=y_test)
fit1 <- sampling(test_stan, data=mode1_data, iter=3000, warmup=0,
								 chains=1, seed=5838298, refresh=1000, algorithm="Fixed_param")
fit_draw_lik_test1=extract(fit1, pars="log_lik_test")$log_lik_test
fit_draw_lik1=extract(fit1, pars="log_lik")$log_lik
fit_draw_f=extract(fit1, pars="f_predict")$f_predict
mode2_data <- list(alpha= opt_2_par[5], rho= opt_2_par[4], sigma= opt_2_par[6], N=length(x), x=x, y=y,N_test=length(x_test), x_test=x_test, y_test=y_test)
fit2 <- sampling(test_stan, data=mode2_data, iter=3000, warmup=0,
									chains=1, seed=5838298, refresh=1000, algorithm="Fixed_param")
fit_draw_lik_test2=extract(fit2, pars="log_lik_test")$log_lik_test
fit_draw_lik2=extract(fit2, pars="log_lik")$log_lik
fit_draw_f2=extract(fit2, pars="f_predict")$f_predict
fit_draw_lik_test=array( NA, c(dim(fit_draw_lik_test1)[1],  2,  dim(fit_draw_lik_test1)[2]))
fit_draw_lik_test[,1,]=fit_draw_lik_test1
fit_draw_lik_test[,2,]=fit_draw_lik_test2
fit_draw_lik=array( NA, c(dim(fit_draw_lik1)[1],  2,  dim(fit_draw_lik1)[2]))
fit_draw_lik[,1,]=fit_draw_lik1
fit_draw_lik[,2,]=fit_draw_lik2

log_lik_test_agg=test_aggregrate(fit_draw_lik_test)
validata_agg=log_lik_test_agg[1:100,]
st_weight=stacking_weights(validata_agg, lambda=1.01)
w_mode= exp (c(opt_1$value, opt_2$value))
w_mode=w_mode/sum(w_mode)

hierstack=hier_stacking(x=x_test[1:100], x_test =x_test[101:300], plpd= validata_agg)
w_test_mat=cbind(colMeans(hierstack$w_test), 1-colMeans(hierstack$w_test))
test_lpd_hier=log_score_loo_pointwise( w_mat=w_test_mat, lpd_point= log_lik_test_agg[101:300,])
log_score_test_mode=c(test_lpd_hier,  log_score_loo(w=st_weight, lpd_point= log_lik_test_agg[101:300,]),
													log_score_loo(w=w_mode, lpd_point= log_lik_test_agg[101:300,]),
													log_score_loo(w=w_is, lpd_point= log_lik_test_agg[101:300,]))
mlpd= rbind(log_score_test_mode,  log_score_test_laplace, log_score_test_unif)/200
colnames(mlpd)=c("hstacking" , "stacking", "mode", "IS" )
print(mlpd) #: the last panel of gp_point.pdf


save(x, y, x_test, y_test, log_score_test_mode,  log_score_test_laplace, log_score_test_unif, file="hier_stacking.RData")


## exact Bayesian inference: in this case NUTS can explore both modes	 
full_bayes_stan=stan_model('predict_gauss.stan')
full_bayes=sampling(full_bayes_stan,
										data=list(N=length(x), x=x, y=y,
															N_predict=length(x_test), x_predict=x_test, y_test=y_test))
full_bayes_extract=extract(full_bayes)
test_lpd=full_bayes_extract$log_lik_test
mean(apply(full_bayes_extract$log_lik_test,  2, log_mean_exp)[101:300])




########## graph 1. ####
par( oma=c(0,1.3,1,0), pty='m',mar=c(2,1,1,1) , lwd=0.5,tck=-0.01, mgp=c(1.5,0.25,0), cex.axis=0.8, cex.lab=0.8, cex.main=0.9,xpd=F)
plot(x_test[1:100],log_lik_test_agg[1:100,2]-log_lik_test_agg[1:100,1],  pch=20, cex=0.8,col=alpha(1, alpha=0.7 ), main="",axes=FALSE,xlab="",ylab="", ylim=c(-1.2,1.2),  xlim=c(-3,3),xaxs='i') 
abline(v=c(-2,-1,0,1,2), col='grey50', lty=2)
abline(h=0,col='darkred')
axis(1,  lwd=0.5, at=c(-3,0,3))
axis(2,  lwd=0.5, at=c(-1,0,1), las=2)
mtext(3, cex=0.7, text = "elpd 1 - elpd 2, pointwise", line=0.8)
box(lwd=0.5, bty='l')
rr=c(101:700)
plot(x_test[rr],1-colMeans(hierstack$w_test)[rr-100],col=1, type='l' , lwd=1,
		 main="",axes=FALSE,xlab="",ylab="", ylim=c(0,1), xlim=c(-3,3),yaxs='i', xaxs='i') 
abline(v=c(-2,-1,0,1,2), col='grey50', lty=2)
lines(x_test[rr],1-colMeans(hierstack$w_test)[rr-100],col=1, lwd=2)
f_9=apply(1-hierstack$w_test[,rr-100],  2, quantile, 0.90)
f_1=apply(1-hierstack$w_test[,rr-100],  2, quantile, 0.1)
f_7=apply(1-hierstack$w_test[,rr-100],  2, quantile, 0.75)
f_2=apply(1-hierstack$w_test[,rr-100],  2, quantile, 0.25)
f_5=apply(1-hierstack$w_test[,rr-100],  2, mean)
polygon(c(x_test[rr] ,rev(x_test[rr]) ), c(f_1 ,rev(f_9)) , col=alpha("grey", alpha = 0.3), border=NA)
polygon(c(x_test[rr] ,rev(x_test[rr]) ), c(f_2 ,rev(f_7)) , col=alpha("grey", alpha = 0.5), border=NA)
lines(x_test[rr],1-colMeans(hierstack$w_test)[rr-100],col=1, lwd=2)
axis(1,  lwd=0.5, at=c(-3,0,3))
axis(2,  lwd=0.5, at=c(0,0.5,1), las=2)
box(lwd=0.5, bty='l')
mtext(1, cex=0.7, text = "x", line=1)
mtext(2, cex=0.7, text = "w", line=1.5, las=2)
legend("bottomright", legend = c("50% CI", "90% CI"),  fill=alpha(c("grey"), alpha=c( 0.8,0.3)), cex=0.8,border = NA,box.lty=0)
mtext(3, cex=0.7, text = "weight of model 1", line=0.3)
abline(h=0.5, col='darkred')


col_vec=c("darkred", "darkgreen", "orange",    "darkblue")
pch_vec=c(8,20,18,15)
par( mar=c(3,2,1,4) )
plot(c(1:3),mlpd[,1], type='n', main="",axes=FALSE,xlab="",ylab=" ", ylim=c( -0.3380185, -0.205))
abline(h=-0.21, col='grey30', lty=2)
text(1.15, -0.216, col='grey30', labels = "exact\n Bayes", cex=0.9)
for(i in 4:1){
	points(c(1:3),mlpd[,i], col=col_vec[i], pch=pch_vec[i], cex=0.9)
	lines(c(1:3),mlpd[,i], col=col_vec[i], lwd=1 ,cex=0.5)
}
axis(1,  lwd=0.5, at=c(1,2,3), labels = c("type-II\n MAP", "Laplace\n approx.", "importance\n resample"),cex.axis=0.9, padj = 0.8 , xpd=T)
axis(2,  lwd=0.5, at=c(-0.32, -0.28, -0.24), las=2)
box(lwd=0.5, bty='l')
#mtext(2, cex=0.7, text = "mean \n test \n lpd", line=1.5, las=2)
mtext(3, cex=0.7, text = "predictive performace of ensemble \n mean  test log predictive density")
text(x=rep(3.45,4),y=c(-0.21, -0.228,-0.242, -0.264), col=col_vec,  cex=1, labels = c("hierarchical\n stacking", "stacking",  "mode height",  "importance\n weighting"), xpd=T)



####### Part 2: covariate shift #########
### we generate another independent test set, this time we use equally-spaced X to reduce variance. The training X is not changed.
### we create 10 bins and 60 data in each bin
set.seed(100)
n_bin=10 # hold-out test data
n_in_bin=60 # hold-out test data
n_test=n_bin*n_in_bin
x_test=seq(-3, 3, length.out = n_test)
f_test=0.3 + 0.4*x_test + 0.5 * sin(2.7*x_test) + 1.1/(1 + x_test^2 )
id.out=c()
count=0
for(i in 1:n_bin){
	id.out=sort(c(id.out, count+sample(1:n_in_bin, 0.05*n_in_bin)))
	count=count+n_in_bin
}
y_test=f_test+rnorm(n_test, 0,0.1)
y_test[id.out]=f_test[id.out]+rnorm(length(id.out), 0,1)
bin_index=rep(1:10, each=n_in_bin)
x_test=c(x_val, x_test) 
y_test=c(y_val, y_test)
n_test=n_test+100

##====laplace approx
grid_normal=function(n=2000){  # generate 3-D N(0,I)
	temp=matrix(NA,  3, n)
	for(i in 1:3)
		temp[i,]=rnorm(n,0,1)
	return(temp)
}
set.seed(2020)
row_grid_nomral=grid_normal()
find_grid=function(row_grid=row_grid_nomral,  opt){
	x1=opt$par[1:3]
	S=dim(row_grid_nomral)[2]
	shift=rbind( rep(x1[1],S), rep(x1[2],S) , rep(x1[3],S) )
	h1=opt$hessian
	cov1=solve(-h1)	
	eig_d= eigen(cov1)
	grid= eig_d$vectors %*% diag( sqrt(eig_d$values))    %*%    row_grid
	grid=t(grid)+t(shift)
	rownames(grid)=c()
	colnames(grid)=names(x1)
	return(grid)
}
grid_Laplace_1=find_grid(row_grid=row_grid_nomral,  opt=opt_1)
grid_Laplace_2=find_grid(row_grid=row_grid_nomral,  opt=opt_2)
S=dim(grid_Laplace_1)[1]
grid_lp1=grid_lp2=c()
for(i in 1: S)
{
	grid_lp1[i]=log_prob(fit_three_pars, upars=grid_Laplace_1[i,]) 
	grid_lp2[i]=log_prob(fit_three_pars, upars=grid_Laplace_2[i,]) 
}
exp(max(grid_lp1)-max(grid_lp2))
m_loo=stan_model("loo.stan")
laplace_gqs1=gqs(m_loo, data=list(N=length(x), x=x, y=y, 
																	N_test=length(x_test), 
																	x_test=x_test, y_test=y_test),
								 draws= expand_draws(grid_Laplace_1) )
laplace_gqs2=gqs(m_loo, data=list(N=length(x), x=x, y=y, 
																	N_test=length(x_test), 
																	x_test=x_test, y_test=y_test),
								 draws= expand_draws(grid_Laplace_2) )
fit_draw_lik1=extract(laplace_gqs1, pars="log_lik")$log_lik
fit_draw_lik2=extract(laplace_gqs2, pars="log_lik")$log_lik
fit_draw_test1=extract(laplace_gqs1, pars="log_lik_test")$log_lik_test
fit_draw_test2=extract(laplace_gqs2, pars="log_lik_test")$log_lik_test
fit_draw_lik=array(NA, c(dim(fit_draw_lik1)[1],  2,    dim(fit_draw_lik1)[2]))
fit_draw_lik[,1,]=fit_draw_lik1
fit_draw_lik[,2,]=fit_draw_lik2
fit_draw_lik_test=array( NA, c(dim(fit_draw_test1)[1],  2,  dim(fit_draw_test1)[2]))
fit_draw_lik_test[,1,]=fit_draw_test1
fit_draw_lik_test[,2,]=fit_draw_test2
# find stacking weights use the validation set
log_lik_test_agg=test_aggregrate(fit_draw_lik_test)
validata_agg=log_lik_test_agg[1:100,]   
st_weight=stacking_weights(validata_agg)

# hierarchical stacking 
hierstack=hier_stacking(x=x_test[1:100], x_test =x_test[101:n_test], plpd= validata_agg)
w_test_mat=cbind(colMeans(hierstack$w_test), 1-colMeans(hierstack$w_test))
test_lpd_hier=log_score_loo_pointwise( w_mat=w_test_mat, lpd_point= log_lik_test_agg[101:n_test,])

#full_bayes_stan=stan_model('predict_gauss.stan')
full_bayes=sampling(full_bayes_stan,
										data=list(N=length(x), x=x, y=y,
															N_predict=length(x_test), x_predict=x_test, y_test=y_test))
full_bayes_extract=extract(full_bayes)
test_lpd=full_bayes_extract$log_lik_test

#== compute point wise lpd
test_lpd_hier_vec_pw=c()
test_lpd_cp_vec_pw=c()
test_lpd_full_bay_pw=apply(full_bayes_extract$log_lik_test,  2, log_mean_exp)
for(i in 1:(n_bin*n_in_bin)){
	test_lpd_hier_vec_pw[i]= weight_log_sum_exp_w(log_lik_test_agg[100+i, ], w_test_mat[i,])
	test_lpd_cp_vec_pw[i]= log(exp(log_lik_test_agg[100+i, ]) %*% st_weight)
}


#== the difference (set h.stacking = benchmark)
win1=test_lpd_hier_vec_pw-test_lpd_full_bay_pw[101:700]
win2=test_lpd_hier_vec_pw-test_lpd_cp_vec_pw

#== the binned mean and sd
win_in_bin_mean=win_in_bin_sd=matrix(NA, n_bin, 2)
for( i in 1:n_bin){
	point.index=which(abs(bin_index-i)<2) 
	win_in_bin_mean[i, 1]=mean(win1[point.index])
	win_in_bin_sd[i,1]=sd(win1[point.index])/sqrt(length(point.index))
	win_in_bin_mean[i, 2]=mean(win2[point.index])
	win_in_bin_sd[i,2]=sd(win2[point.index])/sqrt(length(point.index))
}
bin_center=c() #== x center of each bin
for( i in 1:n_bin){
	bin_center[i]=mean(x_test[100+which(bin_index==i)])
}
#== graph 2. ====
par(mfrow=c(1,2), oma=c(1,2.2,1.1,0), pty='m',mar=c(1,1,1,1) , lwd=0.5,tck=-0.01, mgp=c(1.5,0.25,0), cex.axis=0.75, cex.lab=0.8, cex.main=0.9,xpd=F)
plot(bin_center, win_in_bin_mean[,1], col=1, pch=20 ,  main="",axes=FALSE,xlab="",ylab="", ylim=c(-0.019,0.06))
for( i in 1:n_bin)
	lines(x=rep(bin_center[i],2), y=c( win_in_bin_mean[i,1]-win_in_bin_sd[i,1], win_in_bin_mean[i,1]+win_in_bin_sd[i,1]))
abline(h=0, col='darkred', lwd=0.8)
abline(h=c(0.02,0.04), col='grey50', lty=2)
lines(x_test[rr],1-colMeans(hierstack$w_test)[rr-100],col=1, lwd=2)
axis(1,  lwd=0.5, at=bin_center, labels =NA)
axis(1,  lwd=0.5, at=bin_center[c(1,4,7,10)], labels = round(bin_center[c(1,4,7,10)], 2) ,padj = -0.8   )
axis(2,  lwd=0.5, at=c(-0.02,0,0.02,0.04), las=2)
box(lwd=0.5, bty='l')
mtext(1, cex=0.8, text = "bin center of x", line=0.7)
mtext(2, cex=0.8, text = "test\n lpd\n diff.", line=2, las=2)
mtext(3, cex=0.8, text = "Binned test-lpd difference\n hierarchical stacking minus\n   exact Bayes", line=-0.3)


plot(bin_center, win_in_bin_mean[,2], col=1, pch=20 ,  main="",axes=FALSE,xlab="",ylab="", ylim=c(-0.01,0.2))
for( i in 1:n_bin)
	lines(x=rep(bin_center[i],2), y=c( win_in_bin_mean[i,2]-win_in_bin_sd[i,2], win_in_bin_mean[i,2]+win_in_bin_sd[i,2]))
abline(h=0, col='darkred', lwd=0.8)
abline(h=c(-0.1,0.1,0.2), col='grey50', lty=2)
lines(x_test[rr],1-colMeans(hierstack$w_test)[rr-100],col=1, lwd=2)
axis(1,  lwd=0.5, at=bin_center, labels =NA)
axis(1,  lwd=0.5, at=bin_center[c(1,4,7,10)], labels = round(bin_center[c(1,4,7,10)], 2) , padj = -0.8 )
axis(2,  lwd=0.5, at=c(0,0.1,0.2), las=2)
box(lwd=0.5, bty='l')
mtext(1, cex=0.8, text = "bin center of x", line=0.7)
mtext(3, cex=0.8, text = "Binned test-lpd difference\n hierarchical stacking minus\n complete-pooling stacking", line=-0.3)





