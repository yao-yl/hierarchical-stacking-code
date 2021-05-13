my_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(my_path)
setwd("..")
setwd("..")

## Setup
#rm(list = ls())
# options(mc.cores = 4)
n_chains <- 1
n_cores <- 1
n_sampling <- 500
n_warmup <- 500
n_refresh <- 250

## Libraries
{
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
}


## Functions
source("./scripts/utilities/common_utils.R")
source("./scripts/utilities/get_stan_data2.R")


## Master variables
RUN_DATE      <- ymd("2016-11-08")
start_date    <- as.Date("2016-02-01")
election_date <- ymd("2016-11-08")


## Select models
model_names <- c(
  "poll-model-01",
  "poll-model-02-no-fundamentals",
  "poll-model-03-stock-index",
  "poll-model-04-just-mu-b",
  "poll-model-05-no-sharing-states",
  "poll-model-06-AR2",
  "poll-model-07-no-mode-adjustment",
  "fundamentals-model"
)


## Compile models
models <- list()
for (i in 1:length(model_names)) {
  models[[i]] <- cmdstanr::cmdstan_model(paste0("./scripts/model/stan/", 
                                                model_names[i],
                                                ".stan"), compile = T)
}
names(models) <- model_names


## Forecasting splits
splits1 <- seq(start_date + 14, election_date, by = "week") # dates that separate train and forecast data
splits1 <- splits1 + 4 # So that we are forecasting Friday.

splits2 <- seq(start_date + 14, election_date, by = 14) # dates that separate train and forecast data
splits2 <- splits2 + 4

splits_df1 <- data.frame(x = c(as.Date(NA), splits1), y = c(splits1, NA), type = "week")
splits_df1 <- splits_df1[!is.na(splits_df1$x) & !is.na(splits_df1$y), ]
splits_df2 <- data.frame(x = c(as.Date(NA), splits2), y = c(splits2, NA), type = "2-week")
splits_df2 <- splits_df2[!is.na(splits_df2$x) & !is.na(splits_df2$y), ]
splits_df  <- rbind(splits_df1, splits_df2)
splits_df$ind <- 1:nrow(splits_df)

eg <- expand_grid(splits_df$ind, model_names)


## For parallelization
# nc <- 30 # No. of cores
nc <- detectCores() - 2
if (nrow(eg) < nc) {
  nc <- nrow(eg)
}
mc <- makeCluster(nc)
clusterEvalQ(mc, sink(paste0("./logs/log", Sys.getpid(), ".txt")))
ifelse(nrow(eg) == 1, registerDoSEQ(), registerDoParallel(mc))


## Run models in parallel ------------------------------------------------------
# for (i in 1:nrow(eg)) {
foreach(i = 1:nrow(eg)) %dopar% {
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
  tmp_i <- unlist(eg[i,1])
  print(tmp_i)
  tmp_name       <- paste0(eg[i,2], "-", splits_df[tmp_i,3], "-", splits_df[tmp_i,1], "-", splits_df[tmp_i,2], ".rds")
  print(tmp_name)
  tmp_test_file  <- paste0("./results/forecast-parallel/test-lpd/", tmp_name)
  
  if (file.exists(tmp_test_file)) {
    print("Exists, skipping.")
  } else {
    tmp_run_date      <- splits_df[tmp_i,1]
    tmp_forecast_date <- splits_df[tmp_i,2]
    train_data <- suppressMessages(suppressWarnings(
      get_stan_data2(tmp_run_date, tmp_forecast_date, election_date, start_date)))
    tmp_model <- models[[unlist(eg[i,2])]]
    fit <- tmp_model$sample(
      data = train_data$stan_data,
      seed = 1843,
      parallel_chains = n_cores,
      chains = n_chains,
      iter_warmup = n_warmup,
      iter_sampling = n_sampling,
      refresh = n_refresh
    )
    out <- rstan::read_stan_csv(fit$output_files())
    rm(fit)
    ext <- rstan::extract(out)
    
    if (!is.null(ext$predicted_state_polls)) { # If state polls exist in test
      predicted_state_polls    <- apply(ext$predicted_state_polls, 2, mean)
      lpd_state_polls          <- log(apply(exp(ext$lpd_state_polls), 2, mean))
      state_prediction_df <- data.frame(
        index_s   = train_data$stan_data$state_test,
        poll_day  = train_data$stan_data$day_state_test,
        index_p   = train_data$stan_data$poll_state_test,
        index_m   = train_data$stan_data$poll_mode_state_test,
        index_pop = train_data$stan_data$poll_pop_state_test,
        predicted = predicted_state_polls,
        lpd       = lpd_state_polls
      )
    } else {
      state_prediction_df <- NULL
    }
    if (!is.null(ext$predicted_national_polls)) { # If national polls exist in test
      predicted_national_polls <- apply(ext$predicted_national_polls, 2, mean)
      lpd_national_polls       <- log(apply(exp(ext$lpd_national_polls), 2, mean))
      national_prediction_df <- data.frame(
        index_s   = 52,
        poll_day  = train_data$stan_data$day_national_test,
        index_p   = train_data$stan_data$poll_national_test,
        index_m   = train_data$stan_data$poll_mode_national_test,
        index_pop = train_data$stan_data$poll_pop_national_test,
        predicted = predicted_national_polls,
        lpd       = lpd_national_polls
      )
    } else {
      national_prediction_df <- NULL
    }
    prediction_df     <- rbind(state_prediction_df, national_prediction_df)
    tmp_orig          <- filter(train_data$orig_data, begin > tmp_run_date)
    tmp_orig$poll_day <- as.integer(tmp_orig$poll_day)
    merged   <- left_join(tmp_orig, prediction_df, by = c("index_s",
                                                          "poll_day",
                                                          "index_p",
                                                          "index_m",
                                                          "index_pop")) %>%
      select(-one_of("index_s",
                     "poll_day",
                     "index_p",
                     "index_m",
                     "index_pop",
                     "index_t"))
    merged$model <- eg[i,2]
    
    
    ## Save to file
    saveRDS(merged, tmp_test_file)
  }
}
stopCluster(mc)


