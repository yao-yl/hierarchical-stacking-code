##
##required functions for stacking
library(rstan)

options(mc.cores = parallel::detectCores())

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

# should be equal to
# log_score_loo <- function(w, lpd_point) {
#   sum <- 0
#   N=dim(lpd_point)[1]
#   for (i in 1: N  ) {
#   
#       
#     sum <- sum + log(exp(lpd_point[i, ]) %*% w)
#   }
#   return(as.numeric(sum))
# }


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



  
  
 


stack_with_na=function(lpd_point, lambda=1.0001, lpd_point_test=NULL){
  S=dim(lpd_point)[2]
  N=dim(lpd_point)[1]
  flag=rep(0,S)
  for( i in 1:S)
    if(  sum( is.na( lpd_point[, i])) ==0)
      flag[i]=1
  lpd_point=lpd_point[,which(flag==1)]
  st_weight=stacking_weights( lpd_point=lpd_point, lambda=lambda)
  full_weight=rep(NA,S)
  full_weight[which(flag==1)]=st_weight
  loo_score=log_score_loo(st_weight, lpd_point)
  if (!is.null(lpd_point_test) ){
    lpd_point_test=lpd_point_test[,which(flag==1)]
    test_score=log_score_loo(st_weight, lpd_point_test)
  }
  else
    test_score=NULL
  return(list(loo_score=loo_score, full_weight=full_weight, test_score=test_score, flag=flag))
}


remove_dso_filename=function(stan_model_name){
  dso_filename = stan_model_name@dso@dso_filename
  loaded_dlls = getLoadedDLLs()
  if (dso_filename %in% names(loaded_dlls)) {
    message("Unloading DLL for model dso ", dso_filename)
    model.dll = loaded_dlls[[dso_filename]][['path']]
    dyn.unload(model.dll)
  } else {
    message("No loaded DLL for model dso ", dso_filename)
  }
  
  
  loaded_dlls = getLoadedDLLs()
  loaded_dlls = loaded_dlls[str_detect(names(loaded_dlls), '^file')]
  if (length(loaded_dlls) > 10) {
    for (dll in head(loaded_dlls, -10)) {
      message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
      dyn.unload(dll[['path']])
    }
  }
  message("DLL Count = ", length(getLoadedDLLs()), ": [", str_c(names(loaded_dlls), collapse = ","), "]")
}


stacked_effective_sample_size=function(n_eff_per_chain,full_weight)
{
  full_weight= full_weight/sum( full_weight)
  return(   1/sum(  1/n_eff_per_chain *  full_weight^2)   )
}


# log_marginal_likehood_per_chain=function(fit_draw_lik_per_chain){
#   marginal_likehood_per_chain_per_iter=rowSums(fit_draw_lik_per_chain)
#   return(log_mean_exp (marginal_likehood_per_chain_per_iter))
# }
# 
# log_marginal_likehood_sampled=function(fit_draw_lik){
#   return( apply(fit_draw_lik, 2, log_marginal_likehood_per_chain))
# }

BMA_weight=function(fit_draw_lp){
  log_marginal_likehood=apply(fit_draw_lp,2,log_mean_exp)   
  temp=exp(log_marginal_likehood-max(log_marginal_likehood))
  return(temp/sum(temp))
}

pesudo_BMA_weight_chain=function(fit_draw_lik){
  v=apply(loo_elpd,  2, sum )
  temp=exp(v-max(v))
  return(temp/sum(temp))
  
}


stacking_stan_model=stan_model("optim.stan")

hier_stacking=function(x,x_test, plpd){
  stacking_data=list(
    N=length(x),
    x=x,
    lpd_point= plpd,
    N_test=length(x_test),
    x_test= x_test,
    delta=1e-10,
    gamma_a=4,
    gamma_b=1)
  fit_stacking=sampling(stacking_stan_model, data=stacking_data)
  fit_extaract= extract(fit_stacking)
  w_train=fit_extaract$w_vec
  w_test=fit_extaract$w_test
  return(list(w_train, w_test))
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



