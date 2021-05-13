my_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(my_path)
setwd("..")
setwd("..")

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


lf      <- list.files("./results/forecast-parallel/test-lpd")
all_res <- NULL
for (f in lf) {
  tmp <- readRDS(paste0("./results/forecast-parallel/test-lpd/", f))
  all_res <- bind_rows(all_res, tmp)
}
all_res <- all_res %>%
  mutate(lpd_Bernoulli = (- lchoose((n_trump + n_clinton), n_clinton) + lpd) /
           (n_clinton + n_trump))

saveRDS(all_res, "./results/forecast-parallel/joined/forecast-2016.rds")