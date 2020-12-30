# Prepare data #
# The standard model the and data, adapted from 
# Aki Vehtari and Jonah Gabry's case study 
# https://cran.r-project.org/web/packages/loo/vignettes/loo2-with-rstan.html

## Running on a cluter with arridid 1:280 to repeat training-test split 
## On your local machine set:
args=1

## begin here
libPaths("PATH_TO_LIB")
args <-  Sys.getenv("SLURM_ARRAY_TASK_ID")
print(Sys.getenv("SLURM_ARRAY_JOB_ID"))
arrayid <- as.integer(args[1])
set.seed(as.integer(arrayid))
library(rstan)
library(loo)
library("splines")
S=3000
log_mean_exp=function(v)
{
	max(v)+ log(mean(exp(v- max(v) )))
}
log_score_loo <- function(w, lpd_point) {
	N=dim(lpd_point)[1]
	weight_log_sum_exp<- function(v)
	{
		return(max(v)+ log( exp(v- max(v))%*% w    ))  
	}
	return (sum(apply(lpd_point, 1,  weight_log_sum_exp )) )
}
log_score_loo_vector <- function(w, lpd_point) {
	N=dim(lpd_point)[1]
	weight_log_sum_exp<- function(v)
	{
		return(max(v)+ log( exp(v- max(v))%*% w    ))  
	}
	return(apply(lpd_point, 1,  weight_log_sum_exp )) 
}

m=stan_model("optim_complete_pooling.stan")
stacking_weights=function(lpd_point, lambda=1.0001, stan_model_object=m)
{
	S=dim(lpd_point)[2]
	s_w=optimizing(stan_model_object,  data = list(N=dim(lpd_point)[1], K=S, lpd_point=lpd_point, lambda=rep(lambda, dim(lpd_point)[2])), iter=100000)$par[1:S] 
	return(s_w)
}  

test_aggregrate=function(lpd_point_test)
{
	test_elpd=matrix(NA, ncol=nrow(lpd_point_test[1,,]), nrow=ncol((lpd_point_test[1,,])))
	S=dim(lpd_point_test)[2]
	for(i in 1:S)
		test_elpd[,i]= apply(lpd_point_test[,i,],  2, log_mean_exp)
	return(test_elpd)
}


weight_log_sum_exp_w<- function(v, w)
{
	return(max(v)+ log( exp(v- max(v))%*% w))  
}


log_score_loo_pointwise <- function(w_mat, lpd_point) {
	sum <- 0
	N=dim(lpd_point)[1]
	for (i in 1: N  ) {
		sum <- sum + weight_log_sum_exp_w(lpd_point[i, ], w_mat[i,])
	}
	return(as.numeric(sum))
}


log_score_loo_pointwise_matrix <- function(w_mat, lpd_point) {
	sum=c()
	N=dim(lpd_point)[1]
	for (i in 1: N  ) {
		sum[i]=weight_log_sum_exp_w(lpd_point[i, ], w_mat[i,])
	}
	return(sum)
}



softmax <- function(v){
	v=v-max(v)
	exp(v)/ sum(exp(v))
}



load(file = "wells.RData") # data
n=dim(wells)[1]
#n_train=2000
n_train_v=rep(c(100,200,300,500,700,1000,1200), each=40)
n_train=n_train_v[as.integer(arrayid)] # random split based on input seed
n_test=n-n_train
set.seed(as.integer(arrayid)%%40)
random_permute=sample(1:n)
train_id=random_permute[1:n_train]
test_id=setdiff((1:n), train_id)
set.seed(as.integer(arrayid))
if(n_test>1500) #  test data max n =1500
{
	n_test=1500
	test_id=sample(test_id, n_test)
}



B1 <-as.matrix((bs(wells$logarsenic, knots=quantile(wells$logarsenic[train_id], seq(0.1,0.9,length.out = 10)), degree=3, intercept = FALSE)))  
B2 <-as.matrix((bs(wells$logarsenic, knots=quantile(wells$logarsenic[train_id], seq(0.05,0.95,length.out = 10)), degree=3, intercept = FALSE)))  
B3 <-as.matrix((bs(wells$dist100, knots=quantile(wells$dist100[train_id], seq(0.1,0.9,length.out = 10)), degree=3, intercept = FALSE)))  

X <- model.matrix(~ dist100 + arsenic+assoc+edu1+edu2+edu3,wells[train_id,])
X2 <- model.matrix(~ dist100 + logarsenic+ assoc+edu1+edu2+edu3, wells[train_id,])
X3  <- model.matrix(~ dist100 + arsenic+asthird+ asSquare+ assoc+edu1+edu2+edu3,wells[train_id,])
X4 <- cbind(model.matrix(~ dist100 + assoc+edu1+edu2+edu3,wells[train_id,]),  B1[train_id,])
X5 <- cbind(model.matrix(~ logarsenic + assoc+edu1+edu2+edu3,wells[train_id,]),  B3[train_id,])
X6 <- model.matrix(~ dist100 + logarsenic+ assoc+educ,wells[train_id,])

X_test <- model.matrix(~ dist100 + arsenic+assoc+edu1+edu2+edu3,wells[test_id,])
X2_test <- model.matrix(~ dist100 + logarsenic+ assoc+edu1+edu2+edu3,wells[test_id,])
X3_test <- model.matrix(~ dist100 + logarsenic+asthird+ asSquare+ assoc+edu1+edu2+edu3,wells[test_id,])
X4_test <- cbind(model.matrix(~ dist100 + assoc+edu1+edu2+edu3,wells[test_id,]),  B1[test_id,])
X5_test <- cbind(model.matrix(~ logarsenic + assoc+edu1+edu2+edu3,wells[test_id,]),  B3[test_id,])
X6_test <- model.matrix(~ dist100 + logarsenic+ assoc+educ,wells[test_id,])

dist100_median=median(wells$dist100[train_id] )
arsenic_median=median(wells$logarsenic[train_id] )
wells$dist100_l= pmin(wells$dist100-dist100_median, 0)
wells$dist100_r= pmax(wells$dist100-dist100_median, 0)
wells$arsenic_l=pmin(wells$logarsenic-arsenic_median, 0)
wells$arsenic_r=pmax(wells$logarsenic-arsenic_median, 0)


X_stacking_train <- model.matrix(~ edu0+ edu1+edu2+edu3 + assoc_half+ dist100_l + dist100_r + arsenic_l+ arsenic_r+0, wells[train_id,])
X_stacking_test <- model.matrix(~ edu0+ edu1+edu2+edu3 + assoc_half+ dist100_l + dist100_r + arsenic_l+ arsenic_r+0, wells[test_id,])
K=6
train_x_list=list(X, X2, X3, X4, X5, X6)
train_x_test_list=list(X_test, X2_test, X3_test, X4_test, X5_test, X6_test)
fit_list=list()
for(i in 1:K)
{
	x=train_x_list[[i]]
	x_test=train_x_test_list[[i]]
	fit_list[[i]]= stan("logistic.stan",  chains = 1,  iter = S*2,
											data = list(y = wells[train_id,]$switch, 
																	X = x, N = nrow(x), P = ncol(x), 
																	N_test=  nrow(x_test),  X_test=x_test, 
																	y_test= wells[test_id,]$switch ))
}
find_point_wise_loo_score=function(fit){
	log_lik <- extract_log_lik(fit, merge_chains = FALSE)
	r_eff <- relative_eff(exp(log_lik))
	loo_return=	loo(log_lik, r_eff = r_eff, cores = 1)
	return(loo_return$pointwise)
}
### extract loo
lpd_point=matrix(NA, nrow(X), K)
fit_draw_lik_test=array( NA, c(S,  K,  nrow(X_test) ))
for(i in 1:K)
{
	fit_draw_lik_test[,i,]=extract(fit_list[[i]], pars="log_lik_test")$log_lik_test
	lpd_point[,i]=find_point_wise_loo_score(fit_list[[i]])[,1]
}
log_lik_test=test_aggregrate(fit_draw_lik_test)


### proposed method
fit_stacking <- stan("optim.stan", chain=1, iter=4000,
										 data = list(X = X_stacking_train, N = nrow(X_stacking_train), 
										 						d = ncol(X_stacking_train), d_discrete=4,  
										 						X_test = X_stacking_test, N_test = nrow(X_stacking_test), 
										 						lpd_point=lpd_point, K=ncol(lpd_point),tau_mu=1, tau_sigma=0.5))
m_stacking_mode=stan_model("optim_mode.stan")
f_fit=extract(fit_stacking, pars='f_test')$f_test
S=dim(f_fit)[1]
n_test=dim(f_fit)[2]
w_test=matrix(NA, n_test, K)
for( i in 1:n_test)
{
	temp=matrix(NA, S,K)
	for(s in 1:S){
		temp[s,]=softmax(f_fit[s,i,])
	}
	w_test[i,]=apply(temp, 2, mean)
}

### no-pooling 
fit_stacking_mode <- optimizing(m_stacking_mode, 
																data = list(X = X_stacking_train, N = nrow(X_stacking_train), 
																						d = ncol(X_stacking_train), d_discrete=4, 																						X_test = X_stacking_test, N_test = nrow(X_stacking_test), 
																						lpd_point=lpd_point, K=ncol(lpd_point) ))
name_start=which(names(fit_stacking_mode$par)=="f_test[1,1]")
name_end=which(names(fit_stacking_mode$par)==paste("f_test[",  n_test, ",", K,"]", sep=""))

fit_stacking_mode_mat=matrix(fit_stacking_mode$par[name_start:name_end], n_test, K)
w_test_mode=t(apply(fit_stacking_mode_mat, 1, softmax))


best_model_select=which.max(apply(lpd_point, 2, sum))
cp_stacking=stacking_weights(lpd_point)
v=matrix(NA, n_test, 12)
v[,1]= log_score_loo_pointwise_matrix(w_test, log_lik_test)  # proposed method
v[,2]= log_score_loo_pointwise_matrix(w_test_mode, log_lik_test) # no-pooling 
v[,3]=log_score_loo_vector(cp_stacking,log_lik_test )  # complete-pooling 
v[,4]=log_lik_test[, best_model_select]  # model selection 
v[,5]=best_model_select # index of model selected 


save(v, file=paste("arg_", arrayid, ".RData", sep=""))  
rm(list=ls())


## post-processing and graphing:

col_v=c("darkred", "darkgreen", "darkblue", "darkorange")
method_name=c("hierarchical\n stacking", "stacking", "no-pooling\n stacking", "model\nselection")

n_train_v=rep(c(100,200,300,500,700,1000,1200), each=40)
n=nrow(wells)
repp=40
n_vec=c(100,200,300,500,700,1000,1200)
LLL=length(c(100,200,300,500,700,1000,1200))

## obtain each cluster run:
vm2=array(NA, c(1500, 5, repp, LLL))
k=0
for(i in 1:length(n_train_v)){
	n_train=n_train_v[i]
	n_test=1500
	k=i%%repp
	if(k==0)
		k=repp
	u=ceiling(i/repp)
	file_name=paste("arg_", i, ".RData", sep="")
	if(file.exists(file_name)){
		load(file_name)
		vm2[,, k,u]=v[,c(1,3,2,4,5)]
	}
}
n_test=1500
v1=matrix(NA, LLL,4)
for(u in 1:LLL){	
	v1[u,]=c(mean(vm2[,1,,u], na.rm=T), mean(vm2[,2,,u], na.rm=T), mean(vm2[,3,,u], na.rm=T), mean(vm2[,4,,u], na.rm=T) )
}
plot(v1[,1], ylim=range(v1))
for( i in 2:4)
	lines(v1[,i], col=i)


v1_cal=matrix(NA, LLL, 4) ## compute calibration error
for(u in 1:LLL){
	cal_error=matrix(NA, repp, 4)
	for(j in 1:repp){
		II=20
		prob_interval = seq(0, 1, length.out = II)
		prob_interval[II]=1.01
		cal=array(NA, c(II-1, 2, 4))	
		for(i in 1:4 )
		{
			prob=exp(vm2[,i,j, u]) # emperical error
			for(k in 1:II-1){
				index=(1:n_test)[prob>=prob_interval[k] & prob<prob_interval[k+1]]
				if(length(index)>=1){
					cal[k,1,i]=  mean(prob[index], na.rm=T)
					cal[k,2,i]= 	mean(wells$switch[test_id][index] , na.rm=T)
				}
			}
		}
		cal_error[j,]= apply( abs(cal[,1,]-cal[,2,]), 2, mean, na.rm=T)
	}
	v1_cal[u,]=apply(cal_error, 2, mean,  na.rm=T)
}



v1_worst=v2_worst=matrix(NA, LLL, 4)  ## the most shocking data

for(k in 1:8)
{
	n_worst_vec= n_vec[k]*0.1
	ww=matrix(NA,repp,4)
	for(j in 1:4)
		for(i in 1:repp){
			ww[i,j]=mean(vm2[,j,i,k][order(vm2[,j,i,k]) [1:n_worst_vec]])
		}
	v1_worst[k,]=apply(ww,2, mean, na.rm=T)
	
}

pdf("~/Desktop/well_n.pdf",width=6.9,height=2.1)
layout(matrix(c(1:3),nrow=1), width = c(1,1,1),height = c(1))
par(oma=c(1.5,1,1,0), pty='m',mar=c(1,1,1,2) , lwd=1,tck=-0.01, mgp=c(1.5,0.25,0), cex.axis=0.9, cex.lab=0.8, cex.main=0.9,xpd=F)
plot(n_vec,  v1[,1], col=col_v[1], main="",pch=19, axes=FALSE,xlab="",ylab="", type='n', ylim=c(-0.8,-0.6), xlim=c(0,1220),xaxs='i')
for(i in 1:4){
	points(n_vec,  v1[,i], col=col_v[i], pch=20)
	lines( n_vec,  v1[,i], col=col_v[i])
}
abline(h=seq(-0.75,-0.65, by=0.05), lty=2, col='grey60', lwd=0.5)
mtext(3, cex=0.7, text = "mean test log predictive density\n higher is better", line=-0.5)
axis(1,  lwd=0.5, at=c(0,400,800,1200),  padj = -0.1)
abline(v=seq(200, 1000, by=200), lty=2, col='grey60', lwd=0.4)
axis(2,  lwd=0.5, at=c(-0.8,-0.7,-0.6), las=2)
box(lwd=0.5, bty='l')
mtext(1, cex=0.7, text = "training sample size", line=1.2)


plot(n_vec, v1_cal[,1], ylim=c(0.15,0.26), type='n', col=col_v,  main="",axes=FALSE,xlab="",ylab="", xlim=c(0,1220),xaxs='i' )
for(i in 1:4){
	points(n_vec,  v1_cal[,i], col=col_v[i], pch=20)
	lines( n_vec,  v1_cal[,i], col=col_v[i])
}
abline(h=c(0.2,0.225,0.175), lty=2,col='grey60', lwd=0.5)
mtext(3, cex=0.7, text = "mean absolute calibration error\n lower is better", line=-0.5)
axis(1,  lwd=0.5, at=c(0,400,800,1200),  padj = -0.1)
abline(v=seq(200, 1000, by=200), lty=2, col='grey60', lwd=0.4)
mtext(1, cex=0.7, text = "training sample size", line=1.2)
axis(2,  lwd=0.5, at=c(0.15,0.2, 0.25), las=2)
box(lwd=0.5, bty='l')
par(mar=c(1,1,1,3.2))
plot(n_vec, v1_worst[,1], type='n', main="",pch=19, axes=FALSE,xlab="",ylab=" ", xlim=c(0,1220),xaxs='i', ylim=c(-3,-1))
for(i in 1:4){
	points(n_vec,  v1_worst[,i], col=col_v[i], pch=20)
	lines( n_vec,  v1_worst[,i], col=col_v[i])
}
axis(2,  lwd=0.5, las=2, at=c(-3,  -2, -1))
abline(h=seq(-2.5, -1.5, 0.5), lty=2, col='grey60', lwd=0.4)
axis(1,  lwd=0.5, at=c(0,400,800,1200),  padj = -0.1)
abline(v=seq(200, 1000, by=200), lty=2, col='grey60', lwd=0.4)
mtext(3, cex=0.7, text = "mean test log predictive density\n amomg the worst 10% data points", line=-0.5)
mtext(1, cex=0.7, text = "training sample size", line=1.2)
text(x=c(1330, 1217,950,1280),y=c(-1.1,  -1.45,  -1.93,  -1.8 ), labels = method_name, las=2, col=col_v, xpd=T, cex=0.9)
box(lwd=0.5, bty='l')
dev.off()



