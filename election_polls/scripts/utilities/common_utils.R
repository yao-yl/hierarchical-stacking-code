cov_matrix <- function(n, sigma2, rho){
  m <- matrix(nrow = n, ncol = n)
  m[upper.tri(m)] <- rho
  m[lower.tri(m)] <- rho
  diag(m) <- 1
  (sigma2^.5 * diag(n))  %*% m %*% (sigma2^.5 * diag(n))
}

mean_low_high <- function(draws, states, id){
  tmp <- draws
  draws_df <- data.frame(mean = inv.logit(apply(tmp, MARGIN = 2, mean)),
                         high = inv.logit(apply(tmp, MARGIN = 2, mean) + 1.96 * apply(tmp, MARGIN = 2, sd)),
                         low  = inv.logit(apply(tmp, MARGIN = 2, mean) - 1.96 * apply(tmp, MARGIN = 2, sd)),
                         state = states,
                         type = id)
  return(draws_df)
}

check_cov_matrix <- function(mat,wt=state_weights){
  # get diagnoals
  s_diag <- sqrt(diag(mat))
  # output correlation matrix
  cor_equiv <- cov2cor(mat)
  diag(cor_equiv) <- NA
  # output equivalent national standard deviation
  nat_product <- sqrt(t(wt) %*% mat %*% wt) / 4
  
  # print & output
  hist(as_vector(cor_equiv),breaks = 10)
  
  hist(s_diag,breaks=10)
  
  
  print(sprintf('national sd of %s',round(nat_product,4)))
}