my_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(my_path)
setwd("..")
setwd("..")

library(tidyverse, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(GGally)


lf <- list.files("./results/stacking-parallel/2016-results/")
all_res <- NULL
for (f in lf) {
  tmp_res <- readRDS(paste0("./results/stacking-parallel/2016-results/", f))
  all_res <- dplyr::bind_rows(all_res, tmp_res)
}
df_longer <- all_res %>%
  pivot_wider(names_from = model, values_from= elpd_partial_stacking) %>%
  pivot_longer(-c(state:n_respondents, prior, split),
               names_to = "model", values_to = "lpd")
saveRDS(df_longer, "./results/stacking-parallel/2016-results/joined.rds")

