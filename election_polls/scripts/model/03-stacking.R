my_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(my_path)
setwd("..")
setwd("..")

library(tidyverse, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(rstan, quietly = TRUE)
library(stringr, quietly = TRUE)
library(lubridate, quietly = TRUE)
library(gridExtra, quietly = TRUE)
library(pbapply, quietly = TRUE)
library(parallel, quietly = TRUE)
library(boot, quietly = TRUE)
library(lqmm, quietly = TRUE)
library(gridExtra, quietly = TRUE)
library(ggrepel, quietly = TRUE)
library(doParallel)


options(mc.cores = 4)
n_chains <- 1
n_cores <- 1
n_sampling <- 500
n_warmup <- 500
# n_refresh <- n_sampling*0.25
n_refresh <- 50


softmax <- function (x) {
  x_tilde <- c(x, 0)
  return (exp(x_tilde) / sum(exp(x_tilde)))
}

## Define dates
RUN_DATE      <- ymd("2016-11-08")
start_date    <- as.Date("2016-03-01")
election_date <- ymd("2016-11-08")
upper_date    <- "2016-11-08"


## Define which models to stack
models_to_stack <- list(
  c(
    "fundamentals_model",
    "poll_model_01",
    "poll_model_02_no_fundamentals",
    "poll_model_03_stock_index",
    "poll_model_04_just_mu_b",
    "poll_model_05_no_sharing_states",
    "poll_model_06_AR2",
    "poll_model_07_no_mode_adjustment"
  )
)


## Priors (first is tau_mu, second is tau_sigma)
tau_prior <- matrix(
  c(1, 1),
  ncol = 2,
  byrow = T
)



## Forecasting splits (start later to have some forecasting data for stacking)
splits <- seq(start_date, election_date, by = "week")
splits <- splits + 3


## Read forecasts
df <- readRDS("./results/forecast-parallel/joined/forecast-2016.rds")
df$model <- unlist(df$model)
df$model <- gsub("-", "_", df$model)
res_df <- df %>%
  pivot_wider(names_from = model, values_from = lpd,
              id_cols = c("state",
                          "end",
                          "begin",
                          "pollster",
                          "polltype",
                          "method",
                          "n_respondents"))

res_df <- res_df %>%
  filter(begin <= upper_date)
state_levels <- levels(factor(res_df$state))
res_df <- res_df %>% filter(state != "--")


## Prepare expand.grid for parallel ############################################
model1 <- cmdstanr::cmdstan_model(paste0("scripts/model/stan/", "partial_pooling_weights_softmax",".stan"),compile=T)
model2 <- cmdstanr::cmdstan_model(paste0("scripts/model/stan/", "partial_pooling_weights_softmax_correlations",".stan"),compile=T)


model_names <- c("partial_pooling_softmax",
                 "partial_pooling_softmax_corr",
                 "complete_pooling_standard",
                 "no_pooling_standard",
                 "model_selection")

model_list <- list(model1, model2)
eg_par <- expand.grid(splits[-length(splits)], 1:length(model_names), 
                      1:nrow(tau_prior), 1:length(models_to_stack))
eg_par$Var5 <- splits[-1]
eg_par <- eg_par[ ,c("Var1", "Var2", "Var5", "Var3", "Var4")]
colnames(eg_par) <- c("run_date", "model", "forecast_date", "prior", 
                      "model_split")



## Run parallel ################################################################

## For parallelization
nc <- 30 # No. of cores
mc <- makeCluster(nc)
clusterEvalQ(mc, sink(paste0("./logs/log", Sys.getpid(), ".txt")))
ifelse(nrow(eg_par) == 1, registerDoSEQ(), registerDoParallel(mc))

# for (i in 1:nrow(eg_par)) {
foreach(i = 1:(nrow(eg_par))) %dopar% {
  softmax <- function (x) {
    x_tilde <- c(x, 0)
    return (exp(x_tilde) / sum(exp(x_tilde)))
  }
  library(tidyverse, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(tidyverse, quietly = TRUE)
  library(rstan, quietly = TRUE)
  library(stringr, quietly = TRUE)
  library(lubridate, quietly = TRUE)
  library(gridExtra, quietly = TRUE)
  library(pbapply, quietly = TRUE)
  library(parallel, quietly = TRUE)
  library(boot, quietly = TRUE)
  library(lqmm, quietly = TRUE)
  library(gridExtra, quietly = TRUE)
  library(ggrepel, quietly = TRUE)
  library(doParallel)
  library(cmdstanr)
  print(eg_par[i, ])
  file_name  <- paste0("./results/stacking-parallel/2016-results/", 
                       eg_par[i,1], "_", eg_par[i,3], "_", model_names[eg_par[i,2]], "_prior", eg_par[i,4], "_split", eg_par[i,5], ".rds")
  file_name_w  <- paste0("./results/stacking-parallel/2016-weights/", eg_par[i,1], "_", eg_par[i,3], "_", model_names[eg_par[i,2]], "_prior", eg_par[i,4], "_split", eg_par[i,5], "_weights.rds")
  
  if (file.exists(file_name)) {
    print("File already exists. Skipping.")
  } else {
    tmp_splits <- eg_par[i,1]
    train_df   <- dplyr::filter(res_df, end <= tmp_splits, end >= (tmp_splits - 25))
    test_df    <- dplyr::filter(res_df, begin > tmp_splits, begin <= eg_par[i,3])
    to_stack   <- as.matrix(dplyr::select(train_df, models_to_stack[[eg_par[i,5]]]))
    to_stack2  <- as.matrix(dplyr::select(test_df, models_to_stack[[eg_par[i,5]]]))
    tmp_states <- levels(factor(res_df$state))
    tmp_elpds  <- to_stack
    tmp_elpds2 <- to_stack2
    tmp_g      <- as.numeric(factor(train_df$state, levels = tmp_states))
    tmp_g2     <- as.numeric(factor(test_df$state, levels = tmp_states))
    
    if (eg_par[i,2] %in% 1:2) { # stacking in softmax space
      X          <- matrix(0, ncol = length(tmp_states), nrow = nrow(to_stack))
      tmp_ind    <- as.matrix(data.frame(seq_along(tmp_g), tmp_g))
      X[tmp_ind] <- 1
      X2           <- matrix(0, ncol = length(tmp_states), nrow = nrow(to_stack2))
      tmp_ind2     <- as.matrix(data.frame(seq_along(tmp_g2), tmp_g2))
      X2[tmp_ind2] <- 1
      
      corr_matrix <- readRDS("./results/state-covariance/state-covariance-matrix-2016.rds")
      tmp_ind     <- which(colnames(corr_matrix) == "DC")
      corr_matrix <- corr_matrix[-tmp_ind, -tmp_ind]
      corr_matrix <- cov2cor(corr_matrix)
      
      stan_data <- list(
        N = nrow(to_stack),
        N_test = nrow(to_stack2),
        d = length(tmp_states),
        d_discrete = length(tmp_states),
        K = ncol(to_stack),
        X = X,
        X_test = X2,
        lpd_point = tmp_elpds,
        tau_mu = tau_prior[eg_par[i,4],1],
        tau_sigma = tau_prior[eg_par[i,4],2],
        
        state_correlations = corr_matrix,
        elpds2 = tmp_elpds2
      )
      fit <- model_list[[eg_par[i,2]]]$sample(
        data = stan_data,
        seed = 1843,
        parallel_chains = n_cores,
        chains = n_chains,
        iter_warmup = n_warmup,
        iter_sampling = n_sampling,
        refresh = n_refresh
      )
      out <- rstan::read_stan_csv(fit$output_files())
      ext <- rstan::extract(out)
      test_df$elpd_partial_stacking <- log(apply(exp(ext$elpd2), 2, mean))
      test_df$model <- model_names[eg_par[i,2]]
      test_df$prior <- eg_par[i,4]
      test_df$split <- eg_par[i,5]
      
      if(length(dim(ext$w)) == 2) {
        weights <- apply(ext$w, c(2), mean)
      } else {
        weights <- apply(ext$w, c(2,3), mean)
      }
      if(length(dim(ext$w)) == 2) {
        weights2 <- apply(ext$f, c(2), mean)
        weights2 <- softmax(weights2[-length(weights2)])
      } else {
        weights2 <- apply(ext$f, c(2,3), mean)
        for (j in 1:nrow(weights2)) {weights2[j, ] <- softmax(weights2[j,-ncol(weights2)])}
      }
      weights <- list(weights, weights2)
    } else if (eg_par[i,2]  == 3) { # complete pooling standard
      w <- loo::stacking_weights(to_stack)
      stacked_res <- log(exp(to_stack2) %*% w)
      test_df$elpd_partial_stacking <- stacked_res
      test_df$model <- model_names[eg_par[i,2]]
      test_df$prior <- eg_par[i,4]
      test_df$split <- eg_par[i,5]
      
      weights <- w
    } else if (eg_par[i,2] == 4) { # no pooling standard
      w_mat <- matrix(1 / ncol(to_stack), ncol = ncol(to_stack), nrow = length(tmp_states))
      for (g in unique(tmp_g)) { 
        tmp_elpds1 <- to_stack[tmp_g == g, ]
        if (is.null(nrow(tmp_elpds1))) {
          tmp_elpds1 <- matrix(tmp_elpds1, ncol = length(tmp_elpds1))
        }
        tmp_w <- loo::stacking_weights(tmp_elpds1)
        w_mat[g, ] <- tmp_w
      }
      stacked_res <- NULL
      for (j in 1:length(tmp_g2)) {
        stacked_res <- c(stacked_res, log(exp(to_stack2[j, ]) %*% w_mat[tmp_g2[j], ]))
      }
      test_df$elpd_partial_stacking <- stacked_res
      test_df$model <- model_names[eg_par[i,2]]
      test_df$prior <- eg_par[i,4]
      test_df$split <- eg_par[i,5]
      
      weights <- w_mat
    } else if (eg_par[i,2] == 5) { # model selection
      tmp_apply <- apply(to_stack, 2, mean)
      tmp_ind   <- which.max(tmp_apply)
      stacked_res <- to_stack2[ ,tmp_ind]
      
      test_df$elpd_partial_stacking <- stacked_res
      test_df$model <- model_names[eg_par[i,2]]
      test_df$prior <- eg_par[i,4]
      test_df$split <- eg_par[i,5]
      
      weights <- NULL
    }
    saveRDS(test_df, file_name)
    saveRDS(weights, file_name_w)
  }
}
stopCluster(mc)

