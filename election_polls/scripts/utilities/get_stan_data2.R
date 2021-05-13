get_stan_data2 <- function (RUN_DATE, FORECAST_DATE, election_date = "2016-11-08", start_date = "2016-03-01", which_stock = "DJI_2016") {
  election_day <- ymd(election_date)
  start_date   <- as.Date(start_date) # Keeping all polls after March 1, 2016
  
  # wrangle polls -----------------------------------------------------------
  all_polls <- as_tibble(read.csv("data/all_polls.csv", stringsAsFactors = FALSE, header = TRUE))
  
  # select relevant columns from HufFPost polls
  all_polls <- all_polls %>%
    dplyr::select(state, pollster, number.of.observations, population, mode,
                  start.date,
                  end.date,
                  clinton, trump, undecided, other, johnson, mcmullin) %>%
    filter(ymd(start.date) <= FORECAST_DATE)
  
  # basic mutations
  df <- all_polls %>%
    tbl_df %>%
    rename(n = number.of.observations) %>%
    mutate(begin = ymd(start.date),
           end   = ymd(end.date),
           t = end - (1 + as.numeric(end-begin)) %/% 2) %>%
    filter(begin >= start_date & !is.na(begin)
           & (population == "Likely Voters" |
                population == "Registered Voters" |
                population == "Adults") # get rid of disaggregated polls
           & n > 1)
  
  # pollster mutations
  df <- df %>%
    mutate(pollster = str_extract(pollster, pattern = "[A-z0-9 ]+") %>% sub("\\s+$", "", .),
           pollster = replace(pollster, pollster == "Fox News", "FOX"), # Fixing inconsistencies in pollster names
           pollster = replace(pollster, pollster == "WashPost", "Washington Post"),
           pollster = replace(pollster, pollster == "ABC News", "ABC"),
           pollster = replace(pollster, pollster == "DHM Research", "DHM"),
           pollster = replace(pollster, pollster == "Public Opinion Strategies", "POS"),
           undecided = ifelse(is.na(undecided), 0, undecided),
           other = ifelse(is.na(other), 0, other) +
             ifelse(is.na(johnson), 0, johnson) +
             ifelse(is.na(mcmullin), 0, mcmullin))
  
  # mode mutations
  df <- df %>%
    mutate(mode = case_when(mode == 'Internet' ~ 'Online poll',
                            grepl("live phone",tolower(mode)) ~ 'Live phone component',
                            TRUE ~ 'Other'))
  
  # vote shares etc
  df <- df %>%
    mutate(two_party_sum = clinton + trump,
           polltype = population,
           n_respondents = round(n),
           # clinton
           n_clinton = round(n * clinton/100),
           pct_clinton = clinton/two_party_sum,
           n_trump = round(n * trump/100),
           pct_trump = trump/two_party_sum)
  
  
  ## --- numerical indices
  state_abb_list <- read.csv("data/potus_results_76_16.csv", stringsAsFactors = F) %>%
    pull(state) %>% unique()
  df <- df %>%
    mutate(poll_day = begin - min(begin) + 1,
           # poll_day = t - min(t) + 1,
           # Factors are alphabetically sorted: 1 = --, 2 = AL, 3 = AK, 4 = AZ...
           index_s = as.numeric(factor(as.character(state),
                                       levels = c('--',state_abb_list))),
           index_s = ifelse(index_s == 1, 52, index_s - 1),
           index_t = 1 + as.numeric(begin) - min(as.numeric(begin)),
           index_p = as.numeric(as.factor(as.character(pollster))),
           index_m = as.numeric(as.factor(as.character(mode))),
           index_pop = as.numeric(as.factor(as.character(polltype)))) %>%
    # selections
    arrange(state, t, polltype, two_party_sum) %>%
    distinct(state, t, pollster, .keep_all = TRUE) %>%
    select(
      # poll information
      state, t, begin, end, pollster, polltype, method = mode, n_respondents,
      # vote shares
      pct_clinton, n_clinton,
      pct_trump, n_trump,
      poll_day, index_s, index_p, index_m, index_pop, index_t)
  
  # useful vectors
  all_polled_states <- df$state %>% unique %>% sort
  
  # browser()
  # day indices
  # first_day <- min(df$begin)
  first_day <- start_date
  # ndays <- max(df$begin) - min(df$begin)
  # all_t <- min(df$begin) + days(0:(ndays))
  # all_t_until_election <- min(all_t) + days(0:(election_day - min(all_t)))
  # pollster indices
  all_pollsters <- levels(as.factor(as.character(df$pollster)))
  
  # getting state contextual information from 2012 -------------------------
  states2012 <- read.csv("data/2012.csv",
                         header = TRUE, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    mutate(score = obama_count / (obama_count + romney_count),
           national_score = sum(obama_count)/sum(obama_count + romney_count),
           delta = score - national_score,
           share_national_vote = (total_count*(1+adult_pop_growth_2011_15))
           /sum(total_count*(1+adult_pop_growth_2011_15))) %>%
    arrange(state)
  
  state_abb <- states2012$state
  rownames(states2012) <- state_abb
  
  # get state indices
  all_states <- states2012$state
  state_name <- states2012$state_name
  names(state_name) <- state_abb
  
  # set prior differences
  prior_diff_score <- states2012$delta
  names(prior_diff_score) <- state_abb
  
  # set state weights
  state_weights <- c(states2012$share_national_vote / sum(states2012$share_national_vote))
  names(state_weights) <- state_abb
  
  # electoral votes, by state:
  ev_state <- states2012$ev
  names(ev_state) <- state_abb
  
  
  # create covariance matrices ----------------------------------------------
  # start by reading in data
  state_data <- read.csv("data/potus_results_76_16.csv")
  state_data <- state_data %>%
    select(year, state, dem) %>%
    group_by(state) %>%
    mutate(dem = dem ) %>% #mutate(dem = dem - lag(dem)) %>%
    select(state,variable=year,value=dem)  %>%
    ungroup() %>%
    na.omit() %>%
    filter(variable == 2016)
  
  census <- read.csv('data/acs_2013_variables.csv')
  census <- census %>%
    filter(!is.na(state)) %>%
    select(-c(state_fips,pop_total,pop_density)) %>%
    group_by(state) %>%
    gather(variable,value,
           1:(ncol(.)-1))
  
  state_data <- state_data %>%
    mutate(variable = as.character(variable)) %>%
    bind_rows(census)
  
  # add urbanicity
  urbanicity <- read.csv('data/urbanicity_index.csv') %>%
    dplyr::select(state,pop_density = average_log_pop_within_5_miles) %>%
    gather(variable,value,
           2:(ncol(.)))
  
  state_data <- state_data %>%
    bind_rows(urbanicity)
  
  # add pct white evangelical
  white_evangel_pct <- read_csv('data/white_evangel_pct.csv') %>%
    gather(variable,value,
           2:(ncol(.)))
  
  state_data <- state_data %>%
    bind_rows(white_evangel_pct)
  
  # add region, as a dummy for each region
  regions <- read_csv('data/state_region_crosswalk.csv') %>%
    select(state = state_abb, variable=region) %>%
    mutate(value = 1) %>%
    spread(variable,value)
  
  regions[is.na(regions)] <- 0
  
  regions <- regions %>%
    gather(variable,value,2:ncol(.))
  
  #state_data <- state_data %>%
  #  bind_rows(regions)
  
  # scale and spread
  state_data_long <- state_data %>%
    group_by(variable) %>%
    # scale all varaibles
    mutate(value = (value - min(value, na.rm=T)) /
             (max(value, na.rm=T) - min(value, na.rm=T))) %>%
    #mutate(value = (value - mean(value)) / sd(value)) %>%
    # now spread
    spread(state, value) %>%
    na.omit() %>%
    ungroup() %>%
    select(-variable)
  
  # compute the correlation matrix
  # formula is
  # a*(lambda*C + (1-lambda)*C_1)
  # where C is our correlation matrix with min 0
  # and C_1 is a sq matrix with all 1's
  # lambda=0 is 100% correlation, lambda=1 is our corr matrix
  
  # save correlation
  C <- cor(state_data_long)
  
  # increase the baseline correlation of the matrix to correspond to national-level error
  C[C < 0] <- 0 # baseline cor for national poll error
  
  
  tmp_C <- C
  diag(tmp_C) <- NA
  
  # mixing with matrix of 0.5s
  lambda <- 0.75
  C_1 <- matrix(data=1,nrow = 51,ncol=51)
  a <- 1
  new_C <- (lambda*C + (1-lambda)*C_1) %>% make.positive.definite()
  
  tmp <- new_C
  diag(tmp) <- NA
  
  state_correlation_polling <- new_C
  
  # make pos definite
  state_correlation_polling <- make.positive.definite(state_correlation_polling)
  
  # covariance matrix for polling error
  state_covariance_polling_bias <- cov_matrix(51, 0.078^2, 0.9) # 3.4% on elec day
  state_covariance_polling_bias <- state_covariance_polling_bias * state_correlation_polling
  
  # covariance for prior e-day prediction
  state_covariance_mu_b_T <- cov_matrix(n = 51, sigma2 = 0.18^2, rho = 0.9) # 6% on elec day
  state_covariance_mu_b_T <- state_covariance_mu_b_T * state_correlation_polling
  
  # covariance matrix for random walks
  state_covariance_mu_b_walk <- cov_matrix(51, (0.017)^2, 0.9)
  state_covariance_mu_b_walk <- state_covariance_mu_b_walk * state_correlation_polling # we want the demo correlations for filling in gaps in the polls
  
  ## MAKE DEFAULT COV MATRICES
  # we're going to make TWO covariance matrix here and pass it
  # and 3 scaling values to stan, where the 3 values are
  # (1) the national sd on the polls, (2) the national sd
  # on the prior and (3) the national sd of the random walk
  # make initial covariance matrix (using specified correlation)
  state_covariance_0 <- cov_matrix(51, 0.07^2, 0.9)
  state_covariance_0 <- state_covariance_0 * state_correlation_polling # we want the demo correlations for filling in gaps in the polls
  
  
  # save the inital scaling factor
  national_cov_matrix_error_sd <- sqrt(t(state_weights) %*% state_covariance_0 %*% state_weights) %>% as.numeric()
  
  # save the other scales for later
  fit_rmse_day_x <- function(x){0.03 +  (10^-6.6)*(x)^2} # fit to error from external script
  days_til_election <- as.numeric(difftime(election_day,RUN_DATE))
  expected_national_mu_b_T_error <- fit_rmse_day_x(days_til_election)
  
  polling_bias_scale <- 0.013 # on the probability scale -- we convert later down
  mu_b_T_scale <- expected_national_mu_b_T_error # on the probability scale -- we convert later down
  random_walk_scale <- 0.05/sqrt(300) # on the probability scale -- we convert later down
  
  # gen fake matrices, check the math (this is recreated in stan)
  national_cov_matrix_error_sd <- sqrt(t(state_weights) %*% state_covariance_0 %*% state_weights) %>% as.numeric()
  
  ss_cov_poll_bias = state_covariance_0 * (polling_bias_scale/national_cov_matrix_error_sd*4)^2
  ss_cov_mu_b_T = state_covariance_0 * (mu_b_T_scale/national_cov_matrix_error_sd*4)^2
  ss_cov_mu_b_walk = state_covariance_0 * (random_walk_scale/national_cov_matrix_error_sd*4)^2
  
  # create priors -----------------------------------------------------------
  # read in abramowitz data
  abramowitz <- read.csv('data/abramowitz_data.csv') %>% filter(year < 2016)
  prior_model <- lm(incvote ~  juneapp + q2gdp, data = abramowitz)
  
  # make predictions
  national_mu_prior <- predict(prior_model,newdata = tibble(q2gdp = 1.1, juneapp = 4))
  # on correct scale
  national_mu_prior <- national_mu_prior / 100
  # Mean of the mu_b_prior
  mu_b_prior <- logit(national_mu_prior + prior_diff_score)
  # or read in priors if generated already
  prior_in <- read_csv("data/state_priors_08_12_16.csv") %>%
    filter(date <= RUN_DATE) %>%
    group_by(state) %>%
    arrange(date) %>%
    filter(date == max(date)) %>%
    select(state,pred) %>%
    ungroup() %>%
    arrange(state)
  
  mu_b_prior <- logit(prior_in$pred + 0.0)
  names(mu_b_prior) <- prior_in$state
  national_mu_prior <- weighted.mean(inv.logit(mu_b_prior), state_weights)
  # cat(sprintf('Prior Clinton two-party vote is %s\nWith a national sd of %s\n',
  #             round(national_mu_prior,3),round(mu_b_T_scale,3)))
  
  # alpha for disconnect in state v national polls
  score_among_polled <- sum(states2012[all_polled_states[-1],]$obama_count)/
    sum(states2012[all_polled_states[-1],]$obama_count +
          states2012[all_polled_states[-1],]$romney_count)
  alpha_prior <- log(states2012$national_score[1]/score_among_polled)
  
  # firms that adjust for party (after 2016)
  adjusters <- c(
    "ABC",
    "Washington Post",
    "Ipsos",
    "Pew",
    "YouGov",
    "NBC"
  )
  
  # Stock indexes
  stock_name <- paste0("./data/", which_stock, ".csv")
  stock_data <- as_tibble(read.csv(stock_name, stringsAsFactors = F)) %>%
    select(Date, Close) %>%
    mutate(Close_scaled = scale(Close)[ ,1],
           date = as.Date(Date)) %>%
    select(date, Close_scaled)
  
  FORECAST_DATE <- as.Date(FORECAST_DATE)
  RUN_DATE      <- as.Date(RUN_DATE)
  tmp_date_sequence  <- seq(start_date - 7, FORECAST_DATE, by = "day")
  dates_df <- data.frame(date = tmp_date_sequence)
  last_stock_df <- left_join(dates_df, stock_data, by = "date") %>%
    mutate(day_name = weekdays(date))
  last_stock_df$last_week_close <- 0
  tmp_last <- 0
  for (i in 1:nrow(last_stock_df)) {
    if (last_stock_df$day_name[i] == "Friday") {
      if (!is.na(last_stock_df$Close_scaled[i])) {
        last_stock_df$last_week_close[i] <- tmp_last
        tmp_last <- last_stock_df$Close_scaled[i]
      } else {
        tmp_l   <- T
        tmp_ind <- 1
        while (tmp_l) {
          if (is.na(last_stock_df$Close_scaled[i - tmp_ind])) {
            tmp_ind <- tmp_ind + 1
            next()
          }
          tmp_last <- last_stock_df$Close_scaled[i - tmp_ind]
          tmp_l    <- F
        }
        last_stock_df$last_week_close[i] <- tmp_last
      }
    }
    else {
      last_stock_df$last_week_close[i] <- tmp_last
    }
  }
  time_diff <- as.numeric(FORECAST_DATE - RUN_DATE)
  last_stock_df$previous_running_mean <- 0
  for (i in 1:nrow(last_stock_df)) {
    min_i <- min(abs((i - 4:0)))
    # print(min_i)
    last_stock_df$previous_running_mean[i] <- mean(last_stock_df$Close_scaled[min_i:i], na.rm = T)
    # if (i == 1) browser()
  }
  new_prm <- c(rep(NA, time_diff), last_stock_df$previous_running_mean)
  new_prm <- new_prm[1:nrow(last_stock_df)]
  last_stock_df$prm2 <- new_prm
  prm_at_run <- last_stock_df$prm2[last_stock_df$date == RUN_DATE]
  last_stock_df$prm2[last_stock_df$date >= RUN_DATE] <- prm_at_run
  
  last_stock_df <- filter(last_stock_df, date >= start_date,
                          date <= min(FORECAST_DATE, election_date))
  # MAYBE TOTO: If we do not always forecast on Friday and one-week-ahead,
  # we need to fix the forecasting dates, so that we do not break the
  # time-respecting setting.
  
  # browser()
  
  df_full <- df
  df_test <- filter(df, ymd(begin) > RUN_DATE)
  df      <- filter(df, ymd(end) <= RUN_DATE)
  
  # browser()
  
  # passing data to Stan ----------------------------------------------------
  N_state_polls <- nrow(df %>% filter(index_s != 52))
  N_national_polls <- nrow(df %>% filter(index_s == 52))
  T <- as.integer(round(difftime(election_day, first_day))) + 1
  S <- 51
  P <- length(unique(df_full$pollster))
  M <- length(unique(df_full$method))
  Pop <- length(unique(df_full$polltype))
  
  current_T   <- as.integer(round(difftime(min(FORECAST_DATE, election_date), first_day))) + 1
  stock_index <- last_stock_df$prm2
  
  state <- df %>% filter(index_s != 52) %>% pull(index_s)
  day_national <- df %>% filter(index_s == 52) %>% pull(poll_day)
  day_state <- df %>% filter(index_s != 52) %>% pull(poll_day)
  poll_national <- df %>% filter(index_s == 52) %>% pull(index_p)
  poll_state <- df %>% filter(index_s != 52) %>% pull(index_p)
  poll_mode_national <- df %>% filter(index_s == 52) %>% pull(index_m)
  poll_mode_state <- df %>% filter(index_s != 52) %>% pull(index_m)
  poll_pop_national <- df %>% filter(index_s == 52) %>% pull(index_pop)
  poll_pop_state <- df %>% filter(index_s != 52) %>% pull(index_pop)
  
  n_democrat_national <- df %>% filter(index_s == 52) %>% pull(n_clinton)
  n_democrat_state <- df %>% filter(index_s != 52) %>% pull(n_clinton)
  n_two_share_national <- df %>% filter(index_s == 52) %>% transmute(n_two_share = n_trump + n_clinton) %>% pull(n_two_share)
  n_two_share_state <- df %>% filter(index_s != 52) %>% transmute(n_two_share = n_trump + n_clinton) %>% pull(n_two_share)
  unadjusted_national <- df %>% mutate(unadjusted = ifelse(!(pollster %in% adjusters), 1, 0)) %>% filter(index_s == 52) %>% pull(unadjusted)
  unadjusted_state <- df %>% mutate(unadjusted = ifelse(!(pollster %in% adjusters), 1, 0)) %>% filter(index_s != 52) %>% pull(unadjusted)
  
  
  # Forecast data
  state_test <- df_test %>% filter(index_s != 52) %>% pull(index_s)
  day_national_test <- df_test %>% filter(index_s == 52) %>% pull(poll_day)
  day_state_test <- df_test %>% filter(index_s != 52) %>% pull(poll_day)
  poll_national_test <- df_test %>% filter(index_s == 52) %>% pull(index_p)
  poll_state_test <- df_test %>% filter(index_s != 52) %>% pull(index_p)
  poll_mode_national_test <- df_test %>% filter(index_s == 52) %>% pull(index_m)
  poll_mode_state_test <- df_test %>% filter(index_s != 52) %>% pull(index_m)
  poll_pop_national_test <- df_test %>% filter(index_s == 52) %>% pull(index_pop)
  poll_pop_state_test <- df_test %>% filter(index_s != 52) %>% pull(index_pop)
  
  n_democrat_national_test <- df_test %>% filter(index_s == 52) %>% pull(n_clinton)
  n_democrat_state_test <- df_test %>% filter(index_s != 52) %>% pull(n_clinton)
  n_two_share_national_test <- df_test %>% filter(index_s == 52) %>% transmute(n_two_share = n_trump + n_clinton) %>% pull(n_two_share)
  n_two_share_state_test <- df_test %>% filter(index_s != 52) %>% transmute(n_two_share = n_trump + n_clinton) %>% pull(n_two_share)
  unadjusted_national_test <- df_test %>% mutate(unadjusted = ifelse(!(pollster %in% adjusters), 1, 0)) %>% filter(index_s == 52) %>% pull(unadjusted)
  unadjusted_state_test <- df_test %>% mutate(unadjusted = ifelse(!(pollster %in% adjusters), 1, 0)) %>% filter(index_s != 52) %>% pull(unadjusted)
  
  
  
  # priors (on the logit scale)
  sigma_measure_noise_national <- 0.04
  sigma_measure_noise_state <- 0.04
  sigma_c <- 0.06
  sigma_m <- 0.04
  sigma_pop <- 0.04
  sigma_e_bias <- 0.02
  
  polling_bias_scale <- as.numeric(polling_bias_scale) * 4
  mu_b_T_scale <- as.numeric(mu_b_T_scale) * 4
  random_walk_scale <- as.numeric(random_walk_scale) * 4
  
  # put the data in a list to export to Stan
  data <- list(
    # Number of polls, days, states, etc.
    N_national_polls = N_national_polls,
    N_state_polls = N_state_polls,
    T = T,
    S = S,
    P = P,
    M = M,
    Pop = Pop,
    state = state,
    state_weights = state_weights, # how many votes compared to national sum
    
    
    # Stock indexes
    current_T   = current_T,
    stock_index = stock_index,
    
    
    # Indexes
    day_state = as.integer(day_state),
    day_national = as.integer(day_national),
    poll_state = poll_state,
    poll_national = poll_national,
    poll_mode_national = poll_mode_national,
    poll_mode_state = poll_mode_state,
    poll_pop_national = poll_pop_national,
    poll_pop_state = poll_pop_state,
    unadjusted_national = unadjusted_national, # if the pollster adjusts
    unadjusted_state = unadjusted_state,
    
    
    # Outcomes
    n_democrat_national = n_democrat_national,
    n_democrat_state = n_democrat_state,
    n_two_share_national = n_two_share_national,
    n_two_share_state = n_two_share_state,
    
    
    # Variance scales (set manually) and priors
    sigma_measure_noise_national = sigma_measure_noise_national,
    sigma_measure_noise_state = sigma_measure_noise_state,
    mu_b_prior = mu_b_prior, # based on Abramowitz
    sigma_c = sigma_c,
    sigma_m = sigma_m,
    sigma_pop = sigma_pop,
    sigma_e_bias = sigma_e_bias,
    state_covariance_0 = state_covariance_0, # covariances between states, based on population, region, etc.
    polling_bias_scale = polling_bias_scale,
    mu_b_T_scale = mu_b_T_scale,
    random_walk_scale = random_walk_scale,
    
    
    # Forecast data
    N_national_polls_test = length(day_national_test),
    N_state_polls_test = length(state_test),
    state_test = as.array(state_test),
    day_state_test = as.array(as.integer(day_state_test)),
    day_national_test = as.array(as.integer(day_national_test)),
    poll_state_test = as.array(poll_state_test),
    poll_national_test = as.array(poll_national_test),
    poll_mode_national_test = as.array(poll_mode_national_test),
    poll_mode_state_test = as.array(poll_mode_state_test),
    poll_pop_national_test = as.array(poll_pop_national_test),
    poll_pop_state_test = as.array(poll_pop_state_test),
    n_democrat_national_test = as.array(n_democrat_national_test),
    n_democrat_state_test = as.array(n_democrat_state_test),
    n_two_share_national_test = as.array(n_two_share_national_test),
    n_two_share_state_test = as.array(n_two_share_state_test),
    unadjusted_national_test = unadjusted_national_test, # if the pollster adjusts
    unadjusted_state_test = unadjusted_state_test
  )
  orig_data <- df_full
  return (list(stan_data = data,
               orig_data = orig_data))
}
