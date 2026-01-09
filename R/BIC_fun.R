BIC_fun <- function(y, yfit, K0, Ky, Kx) {
  n <- nrow(y)
  
  # Compute squared residual sum for each observation
  arg_bic <- numeric(n)
  for (m in 1:n) {
    arg_bic[m] <- t(y[m, ] - yfit[m, ]) %*% (y[m, ] - yfit[m, ])
  }
  
  # Bayesian Information Criterion
  penalty_term <- (K0 + Ky + Kx + 2) * log(n)
  BIC_val <- n * log(sum(arg_bic) / n) + penalty_term
  
  return(BIC_val)
}
