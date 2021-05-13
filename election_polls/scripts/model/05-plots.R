my_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(my_path)
setwd("..")
setwd("..")

library(tidyverse, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(GGally)
library(ggrepel)
library(lubridate)

folder_name <- "2016"
ndays       <- 25
df_orig <- readRDS(paste0("./results/stacking-parallel/",
                          folder_name,
                          "-results/joined.rds"))

## Weekly ----------------------------------------------------------------------
df <- df_orig
df <- df %>%
  filter(model %in% c("poll_model_01", "complete_pooling_standard", 
                      "no_pooling_standard",
                      "partial_pooling_softmax_corr", 
                      "partial_pooling_softmax", "model_selection"),
         prior == 1)
df$model <- factor(df$model,
                   levels = c(
                     "poll_model_01", "complete_pooling_standard", 
                     "no_pooling_standard",
                     "partial_pooling_softmax_corr", 
                     "partial_pooling_softmax", "model_selection"
                   ))
df$model <- as.character(df$model)

df_ref <- df %>%
  pivot_wider(names_from = model, values_from = lpd) %>%
  mutate(ref_lpd = partial_pooling_softmax) %>%
  pivot_longer(c("poll_model_01","no_pooling_standard", 
                 "complete_pooling_standard",
                 "partial_pooling_softmax_corr", 
                 "partial_pooling_softmax", "model_selection"),
               names_to = "model", values_to = "lpd") %>%
  mutate(lpd_diff = lpd - ref_lpd)


agg_df <- df %>%
  group_by(model, split, prior) %>%
  summarize(mean_lpd = mean(lpd),
            SE_lpd = sd(lpd) / sqrt(n()))


agg_df <- df_ref %>%
  group_by(model, split, prior) %>%
  summarize(mean_lpd_diff = mean(lpd_diff),
            SE_lpd_diff = sd(lpd_diff) / sqrt(n()))
agg_df$model[agg_df$model == "complete_pooling_standard"] <- "stacking"
agg_df$model[agg_df$model == "no_pooling_standard"] <- "no-pooling \n stacking"
agg_df$model[agg_df$model == "model_selection"] <- "model \n selection"
agg_df$model[agg_df$model == "partial_pooling_softmax_corr"] <- "hierarchical \nstacking corr"
agg_df$model[agg_df$model == "partial_pooling_softmax"] <- "hierarchical stacking"

agg_df2 <- agg_df %>%
  filter(!(model %in% c("hierarchical stacking", "poll_model_01")))

agg_df2$model <- factor(agg_df2$model, levels = c("hierarchical \nstacking corr", 
                                                  "stacking", "no-pooling \n stacking",
                                                  "model \n selection"))

df_sub <- df
df_sub$running_mean_lpd <- NA
for (i in 1:nrow(df_sub)) {
  tmp_mod   <- df_sub$model[i]
  tmp_split <- df_sub$split[i]
  tmp_prior <- df_sub$prior[i]
  tmp_end   <- df_sub$end[i]
  df_sub$running_mean_lpd[i] <- mean(df_sub$lpd[df_sub$model == tmp_mod &
                                                  df_sub$split == tmp_split &
                                                  df_sub$prior == tmp_prior &
                                                  df_sub$end <= tmp_end &
                                                  df_sub$end > tmp_end - 7])
}


# Get differences!
df_wider <- pivot_wider(df_sub, names_from = model, values_from = running_mean_lpd,
                        id_cols = state:split) %>%
  mutate(complete_pooling = complete_pooling_standard - partial_pooling_softmax,
         partial_pooling = partial_pooling_softmax - partial_pooling_softmax,
         partial_pooling_corr = partial_pooling_softmax_corr - partial_pooling_softmax,
         no_pooling = no_pooling_standard - partial_pooling_softmax,
         model_selection2 = model_selection - partial_pooling_softmax)

df_longer <- df_wider %>%
  pivot_longer(complete_pooling:model_selection2, names_to = "model")



df_longer$model[df_longer$model == "complete_pooling"] <- "stacking"
df_longer$model[df_longer$model == "no_pooling"] <- "no-pooling \n stacking"
df_longer$model[df_longer$model == "model_selection2"] <- "model \n selection"
df_longer$model[df_longer$model == "partial_pooling_corr"] <- "hierarchical \nstacking corr"
df_longer$model[df_longer$model == "partial_pooling"] <- "hierarchical stacking"

df_longer <- df_longer %>%
  filter(!(model %in% c("hierarchical stacking", "poll_model_01")))

df_longer$model <- factor(df_longer$model, levels = c("hierarchical \nstacking corr", "stacking", "no-pooling \n stacking",
                                                      "model \n selection"))

df_sub2 <- df_longer


my_df2 <- df_sub2 %>%
  filter(split == 1, end >= "2016-04-01") %>%
  mutate(label = if_else(end == max(end), as.character(model), NA_character_)) %>%
  select(end, label, value, model) %>%
  unique()

# OK, I'd do the full one.
p1 <- ggplot(filter(df_sub2, end >= "2016-03-01", split == 1), aes(x = end, y = value, color = model)) +
  geom_line(size = 0.4) +
  geom_line(data = filter(df_sub2, end >= "2016-03-01", split == 1, model == "stacking"), aes(x = end, y = value, color = model), size = 0.4) +
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 7),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line()) +
  geom_vline(xintercept = as.Date("2016-11-07", format = "%Y-%m-%d"),
             linetype = "dashed",
             alpha = 0.2) +
ylab("cumulative mean test lpd (higher better)") +
  scale_color_manual(values = c("purple", "darkgreen", "darkblue", "darkorange")) +
  geom_segment(aes(y = 0, yend = 0, x = as.Date("2016-03-01", format = "%Y-%m-%d"),
                   xend = as.Date("2016-11-07", format = "%Y-%m-%d")), color = "darkred",
               size = 0.4) +
  xlab("date") +
  ylab("") +
  labs(subtitle = "pointwise differences in 7-days running mean \ntest log predictive density") +
  coord_cartesian(xlim = c(as.Date("2016-03-01", format = "%Y-%m-%d"),
                           as.Date("2016-11-07", format = "%Y-%m-%d")), 
                  # ylim = c(-0.5,0.1),
                  clip = 'off') +
  scale_x_date(breaks = c(as.Date("2016-04-01", format = "%Y-%m-%d"),
                          as.Date("2016-07-01", format = "%Y-%m-%d"),
                          as.Date("2016-10-01", format = "%Y-%m-%d"),
                          as.Date("2016-11-07", format = "%Y-%m-%d")),
               labels = c("Apr", "Jul", "Oct", "election"))

## Cumulative ------------------------------------------------------------------
df <- df_orig
df <- df %>%
  filter(model %in% c("poll_model_01", "complete_pooling_standard", 
                      "no_pooling_standard",
                      "partial_pooling_softmax_corr", 
                      "partial_pooling_softmax", "model_selection"),
         prior == 1)
df$model <- factor(df$model,
                   levels = c(
                     "poll_model_01", "complete_pooling_standard", 
                     "no_pooling_standard",
                     "partial_pooling_softmax_corr", 
                     "partial_pooling_softmax", "model_selection"
                   ))
df$model <- as.character(df$model)

df_ref <- df %>%
  pivot_wider(names_from = model, values_from = lpd) %>%
  mutate(ref_lpd = partial_pooling_softmax) %>%
  pivot_longer(c("poll_model_01","no_pooling_standard", "complete_pooling_standard",
                 "partial_pooling_softmax_corr", 
                 "partial_pooling_softmax", "model_selection"),
               names_to = "model", values_to = "lpd") %>%
  mutate(lpd_diff = lpd - ref_lpd)


agg_df <- df %>%
  group_by(model, split, prior) %>%
  summarize(mean_lpd = mean(lpd),
            SE_lpd = sd(lpd) / sqrt(n()))


agg_df <- df_ref %>%
  group_by(model, split, prior) %>%
  summarize(mean_lpd_diff = mean(lpd_diff),
            SE_lpd_diff = sd(lpd_diff) / sqrt(n()))
agg_df$model[agg_df$model == "complete_pooling_standard"] <- "stacking"
agg_df$model[agg_df$model == "no_pooling_standard"] <- "no-pooling \n stacking"
agg_df$model[agg_df$model == "model_selection"] <- "model \n selection"
agg_df$model[agg_df$model == "partial_pooling_softmax_corr"] <- "hierarchical \nstacking corr"
agg_df$model[agg_df$model == "partial_pooling_softmax"] <- "hierarchical stacking"

agg_df2 <- agg_df %>%
  filter(!(model %in% c("hierarchical stacking", "poll_model_01")))

agg_df2$model <- factor(agg_df2$model, levels = c("hierarchical \nstacking corr", 
                                                  "stacking", "no-pooling \n stacking",
                                                  "model \n selection"))

df_sub <- df
df_sub$running_mean_lpd <- NA
for (i in 1:nrow(df_sub)) {
  tmp_mod   <- df_sub$model[i]
  tmp_split <- df_sub$split[i]
  tmp_prior <- df_sub$prior[i]
  tmp_end   <- df_sub$end[i]
  df_sub$running_mean_lpd[i] <- mean(df_sub$lpd[df_sub$model == tmp_mod &
                                                  df_sub$split == tmp_split &
                                                  df_sub$prior == tmp_prior &
                                                  df_sub$end <= tmp_end])
}

# Get differences!
df_wider <- pivot_wider(df_sub, names_from = model, values_from = running_mean_lpd,
                        id_cols = state:split) %>%
  mutate(complete_pooling = complete_pooling_standard - partial_pooling_softmax,
         partial_pooling = partial_pooling_softmax - partial_pooling_softmax,
         partial_pooling_corr = partial_pooling_softmax_corr - partial_pooling_softmax,
         no_pooling = no_pooling_standard - partial_pooling_softmax,
         model_selection2 = model_selection - partial_pooling_softmax)

df_longer <- df_wider %>%
  pivot_longer(complete_pooling:model_selection2, names_to = "model")



df_longer$model[df_longer$model == "complete_pooling"] <- "stacking"
df_longer$model[df_longer$model == "no_pooling"] <- "no-pooling \n stacking"
df_longer$model[df_longer$model == "model_selection2"] <- "model \n selection"
df_longer$model[df_longer$model == "partial_pooling_corr"] <- "hierarchical \nstacking corr"
df_longer$model[df_longer$model == "partial_pooling"] <- "hierarchical stacking"

df_longer <- df_longer %>%
  filter(!(model %in% c("hierarchical stacking", "poll_model_01")))

df_longer$model <- factor(df_longer$model, levels = c("hierarchical \nstacking corr", 
                                                      "stacking", "no-pooling \n stacking",
                                                      "model \n selection"))

df_sub2 <- df_longer


my_df2 <- df_sub2 %>%
  filter(split == 1, end >= "2016-04-01") %>%
  mutate(label = if_else(end == max(end), as.character(model), NA_character_)) %>%
  select(end, label, value, model) %>%
  unique()

p2 <- ggplot(filter(df_sub2, end >= "2016-04-01", split == 1), aes(x = end, y = value, color = model)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 7),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line()) +
  geom_vline(xintercept = as.Date("2016-11-07", format = "%Y-%m-%d"),
             linetype = "dashed",
             alpha = 0.2) +
  geom_label_repel(aes(label = label, color = label),
                   data = my_df2,
                   nudge_x = 13,
                   label.size=NA,
                   size = 2,
                   label.padding=.1, 
                   min.segment.length = 0.1,
                   na.rm=TRUE,
                   fill = alpha(c("white"),0.5)) +
  ylab("cumulative mean test lpd (higher better)") +
  scale_color_manual(values = c("purple", "darkorange", "darkblue", "darkgreen")) +
  geom_segment(aes(y = 0, yend = 0, x = as.Date("2016-04-01", format = "%Y-%m-%d"),
                   xend = as.Date("2016-11-07", format = "%Y-%m-%d")), color = "darkred") +
  annotate("text", x=as.Date("2016-07-01", format = "%Y-%m-%d"), y=0.03, label="hierarchical stacking = 0", color = "darkred",
           size = 2.5) +
  xlab("date") +
  ylab("") +
  labs(subtitle = "pointwise differences in mean cumulative \ntest log predictive density") +
  coord_cartesian(xlim = c(as.Date("2016-04-01", format = "%Y-%m-%d"),
                           as.Date("2016-12-30", format = "%Y-%m-%d")), 
                  clip = 'off') +
  scale_x_date(breaks = c(as.Date("2016-04-01", format = "%Y-%m-%d"),
                          as.Date("2016-07-01", format = "%Y-%m-%d"),
                          as.Date("2016-10-01", format = "%Y-%m-%d"),
                          as.Date("2016-11-07", format = "%Y-%m-%d")),
               labels = c("Apr", "Jul", "Oct", "election"))
egg::ggarrange(p1, p2, nrow = 1)



## Cumulative by state ---------------------------------------------------------
df <- df_orig
df <- df %>%
  filter(model %in% c("poll_model_01", "complete_pooling_standard", 
                      "no_pooling_standard",
                      "partial_pooling_softmax_corr",
                      "partial_pooling_softmax", "model_selection"),
         prior == 1,
         split == 1)
df$model <- factor(df$model,
                   levels = c(
                     "poll_model_01", "complete_pooling_standard", 
                     "no_pooling_standard",
                     "partial_pooling_softmax_corr",
                     "partial_pooling_softmax", "model_selection"
                   ))
df$model <- as.character(df$model)

df_ref <- df %>%
  pivot_wider(names_from = model, values_from = lpd) %>%
  mutate(ref_lpd = complete_pooling_standard) %>%
  pivot_longer(c("poll_model_01","no_pooling_standard", "complete_pooling_standard",
                 "partial_pooling_softmax_corr",
                 "partial_pooling_softmax", "model_selection"),
               names_to = "model", values_to = "lpd") %>%
  mutate(lpd_diff = lpd - ref_lpd)


agg_df <- df %>%
  group_by(model, split, prior) %>%
  summarize(mean_lpd = mean(lpd),
            SE_lpd = sd(lpd) / sqrt(n()))


agg_df <- df_ref %>%
  group_by(model, split, prior) %>%
  summarize(mean_lpd_diff = mean(lpd_diff),
            SE_lpd_diff = sd(lpd_diff) / sqrt(n()))

df_sub <- df

df_sub$running_mean_lpd <- NA
for (i in 1:nrow(df_sub)) {
  tmp_mod   <- df_sub$model[i]
  tmp_split <- df_sub$split[i]
  tmp_prior <- df_sub$prior[i]
  tmp_end   <- df_sub$end[i]
  tmp_state <- df_sub$state[i]
  df_sub$running_mean_lpd[i] <- mean(df_sub$lpd[df_sub$model == tmp_mod &
                                                  df_sub$split == tmp_split &
                                                  df_sub$prior == tmp_prior &
                                                  df_sub$end <= tmp_end &
                                                  df_sub$state == tmp_state])
}

# Get differences!
df_wider <- pivot_wider(df_sub, names_from = model, values_from = running_mean_lpd,
                        id_cols = state:split) %>%
  mutate(complete_pooling = complete_pooling_standard - partial_pooling_softmax,
         partial_pooling = partial_pooling_softmax - partial_pooling_softmax,
         partial_pooling_corr = partial_pooling_softmax_corr - partial_pooling_softmax,
         no_pooling = no_pooling_standard - partial_pooling_softmax,
         model_selection2 = model_selection - partial_pooling_softmax)

df_longer <- df_wider %>%
  pivot_longer(complete_pooling:model_selection2, names_to = "model")


df_sub2 <- df_longer

df_sub3 <- df %>%
  filter(end >= "2016-04-01", model == "poll_model_01") %>%
  group_by(state) %>%
  mutate(n = n()) %>%
  select(state, n) %>%
  unique() %>%
  arrange(n)

df_sub2$model[df_sub2$model == "complete_pooling"] <- "stacking"
df_sub2$model[df_sub2$model == "no_pooling"] <- "no-pooling \nstacking"
df_sub2$model[df_sub2$model == "model_selection2"] <- "model \nselection"
df_sub2$model[df_sub2$model == "partial_pooling_corr"] <- "hierarchical \nstacking corr"
df_sub2$model[df_sub2$model == "partial_pooling"] <- "hierarchical \nstacking"

df_sub2$model <- factor(df_sub2$model, levels = c("hierarchical \nstacking corr", 
                                                  "hierarchical \nstacking",
                                                  "stacking", "no-pooling \nstacking",
                                                  "model \nselection"))

df_sub2 <- left_join(df_sub2, df_sub3) %>%
  mutate(staten = paste0(state, " (n=", n, ")"))

tmp_df <- df_sub2 %>%
  select(state, n) %>%
  unique() %>%
  arrange(n)
df_sub$state <- factor(df_sub$state, levels = tmp_df$state)

df_paper <- df_sub2 %>%
  dplyr::filter(state %in% c("GA", "FL", "NC", "RI", "SD", "WV"), split == 1, end >= "2016-04-01") %>%
  filter(model %in% c("stacking", "no-pooling \nstacking", "model \nselection", "hierarchical \nstacking", "hierarchical \nstacking corr"))
tmp_df <- df_paper %>%
  select(staten, n) %>%
  unique() %>%
  arrange(n)
df_paper$staten <- factor(df_paper$staten, levels = tmp_df$staten)


ggplot(df_paper, aes(x = end, y = value, color = model)) +
  geom_line(size = 0.4) +
  theme_classic() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(linetype = "dashed"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 7),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank()) +
  ylab("cumulative mean test lpd (higher better)") +
  scale_color_manual(values = c("purple", "darkred", "darkgreen", "darkblue", "darkorange")) +
  xlab("date") +
  ylab("") +
  labs(subtitle = "pointwise differences in mean cumulative test log predictive density by state",
       color = "method") +
  facet_wrap(~ staten, nrow = 1, scales = "free_y") +
  scale_x_date(breaks = c(as.Date("2016-04-01", format = "%Y-%m-%d"),
                          as.Date("2016-07-01", format = "%Y-%m-%d"),
                          as.Date("2016-10-01", format = "%Y-%m-%d")),
               labels = c("Apr", "Jul", "Oct"))


df_paper2 <- df_sub2 %>%
  dplyr::filter(split == 1, end >= "2016-04-01") %>%
  filter(model %in% c("stacking", "no-pooling \nstacking", "model \nselection", 
                      "hierarchical \nstacking", "hierarchical \nstacking corr"))
tmp_df <- df_paper2 %>%
  select(staten, n) %>%
  unique() %>%
  arrange(n)
df_paper2$staten <- factor(df_paper2$staten, levels = tmp_df$staten)

ggplot(df_paper2, aes(x = end, y = value, color = model)) +
  geom_line(size = 0.4) +
  theme_classic() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(linetype = "dashed"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 7),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank()) +
  ylab("cumulative mean test lpd (higher better)") +
  scale_color_manual(values = c("purple", "darkred", "darkgreen", "darkblue", "darkorange")) +
  xlab("date") +
  ylab("") +
  labs(subtitle = "pointwise differences in mean cumulative test log predictive density by state (higher is better)",
       color = "method") +
  facet_wrap(~ staten, ncol = 7, scales = "free_y") +
  scale_x_date(breaks = c(as.Date("2016-04-01", format = "%Y-%m-%d"),
                          as.Date("2016-07-01", format = "%Y-%m-%d"),
                          as.Date("2016-10-01", format = "%Y-%m-%d")),
               labels = c("Apr", "Jul", "Oct"))


## By size ---------------------------------------------------------------------
df <- readRDS("./results/stacking-parallel/2016-results/joined.rds")
df <- df %>%
  filter(model %in% c("poll_model_01", "complete_pooling_standard", 
                      "no_pooling_standard",
                      "partial_pooling_softmax_corr", 
                      "partial_pooling_softmax", "model_selection"),
         prior == 1)
df$model <- factor(df$model,
                   levels = c(
                     "poll_model_01", "complete_pooling_standard", "no_pooling_standard",
                     "partial_pooling_softmax_corr", 
                     "partial_pooling_softmax", "model_selection"
                   ))
df$model <- as.character(df$model)

df_ref <- df %>%
  pivot_wider(names_from = model, values_from = lpd) %>%
  mutate(ref_lpd = partial_pooling_softmax_corr) %>%
  pivot_longer(c("poll_model_01","no_pooling_standard", "complete_pooling_standard",
                 "partial_pooling_softmax_corr", 
                 "partial_pooling_softmax", "model_selection"),
               names_to = "model", values_to = "lpd") %>%
  mutate(lpd_diff = lpd - ref_lpd)

df_group <- df_ref %>%
  filter(model == "poll_model_01", prior == 1, split == 1) %>%
  group_by(state) %>%
  summarize(n = n())

df_ref <- left_join(df_ref, df_group)

lu <- quantile(df_group$n, probs = c(0.33, 0.66))

df_ref$state_size <- "few~polls~(n~phantom()<~15)"
df_ref$state_size[df_ref$n > lu[1]] <- "moderate~polls~(15~phantom()<=~n~phantom()<~25)"
df_ref$state_size[df_ref$n > lu[2]] <- "many~polls~(25~phantom()<=~n)"

df_ref2 <- df_ref
df_ref2$state_size <- "all~polls"

df_ref <- rbind(df_ref, df_ref2)

df_ref$state_size <- factor(df_ref$state_size, levels = c(
  "few~polls~(n~phantom()<~15)",
  "moderate~polls~(15~phantom()<=~n~phantom()<~25)",
  "many~polls~(25~phantom()<=~n)",
  "all~polls"
))


agg_df <- df_ref %>%
  group_by(model, split, prior, state_size) %>%
  summarize(mean_lpd_diff = mean(lpd_diff),
            SE_lpd_diff = sd(lpd_diff) / sqrt(n()))
agg_df$model[agg_df$model == "complete_pooling_standard"] <- "stacking"
agg_df$model[agg_df$model == "no_pooling_standard"] <- "no-pooling \n stacking"
agg_df$model[agg_df$model == "model_selection"] <- "model \nselection"
agg_df$model[agg_df$model == "partial_pooling_softmax_corr"] <- "hierarchical \nstacking corr"
agg_df$model[agg_df$model == "partial_pooling_softmax"] <- "hierarchical \nstacking"

agg_df2 <- agg_df %>%
  filter(!(model %in% c("hierarchical \nstacking corr", "poll_model_01")))

agg_df2$model <- factor(agg_df2$model, levels = c("hierarchical \nstacking", 
                                                  "stacking", "no-pooling \n stacking",
                                                  "model \nselection"))

annot_df <- data.frame(x = 2.5, y = 0.015, state_size = factor("all~polls"),
                       levels = c(
                         "few~polls~(n~phantom()<~15)",
                         "moderate~polls~(15~phantom()<=~n~phantom()<~25)",
                         "many~polls~(25~phantom()<=~n)",
                         "all~polls"
                       ))
p1 <- ggplot(filter(agg_df2, split == 1, model != "poll_model_01"), 
             aes(x = model, y = mean_lpd_diff)) +
  geom_hline(yintercept = 0, color = "purple") +
  geom_point(aes(color = model), size = 1) +
  geom_errorbar(aes(ymin = mean_lpd_diff - 1.96 * SE_lpd_diff,
                    ymax = mean_lpd_diff + 1.96 * SE_lpd_diff,
                    color = model),
                width = 0,
                lwd = 0.5) +
  geom_errorbar(aes(ymin = mean_lpd_diff - 0.68 * SE_lpd_diff,
                    ymax = mean_lpd_diff + 0.68 * SE_lpd_diff,
                    color = model),
                width = 0,
                lwd = 1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylab("pointwise test lpd differences to reference model") +
  ylab("") +
  xlab("") +
  scale_color_manual(values=c("darkred", "darkgreen", "darkblue", "darkorange")) +
  geom_text(data = annot_df, mapping = aes(x = x, y = y), 
            label = "hierarchical stacking corr = 0",
            size = 3, check_overlap = TRUE, color = "purple") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_line(linetype = "dashed"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 11.8),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        strip.background =element_rect(fill="white"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,10,10,10)) +
  labs(subtitle = "mean pointwise difference in test log predictive density\nby number of polls in state",
       color = "method") +
  facet_wrap(~ state_size, nrow = 1, labeller=label_parsed) +
  coord_cartesian(ylim = c(-0.2, 0.025))

