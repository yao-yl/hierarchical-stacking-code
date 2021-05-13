my_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(my_path)
setwd("..")
setwd("..")

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
}

## Functions
source("./scripts/utilities/common_utils.R")
source("./scripts/utilities/get_stan_data.R")

## Run on some arbitrary data to get the correlation matrices
start_date    <- as.Date("2016-02-01")
election_date <- ymd("2016-11-08")
my_data <- get_stan_data(RUN_DATE = ymd("2016-10-01"), FORECAST_DATE = ymd("2016-10-08"),
                         election_date,
                         start_date)
my_cov <- my_data$stan_data$state_covariance_0
saveRDS(my_cov, "./results/state-covariance/state-covariance-matrix-2016.rds")
