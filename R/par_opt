par_opt <- function(y, x, W, gpy, gpx, K0, Ky, Kx, lam_cands) {
  
  # Candidate lambda values
  lb_vec   <- lam_cands$lb
  lrho_vec <- lam_cands$lrho
  
  # Create grid of (lb, lrho) pairs
  BICmat_l <- as.matrix(expand.grid(lb_vec, lrho_vec))
  colnames(BICmat_l) <- c("lb", "lrho")
  
  # Store BIC values for each combination
  BIC_rho <- numeric(nrow(BICmat_l))
  
  # Evaluate BIC for each combination
  for (j in seq_len(nrow(BICmat_l))) {
    
    model_j <- pen2SLS_est(
      y, x, W, gpy, gpx, K0, Ky, Kx, 
      lb    = BICmat_l[j, 1], 
      lrho  = BICmat_l[j, 2]
    )
    
    yhat_j <- model_j$fitted.values
    BIC_rho[j] <- BIC_fun(y, yhat_j, K0, Ky, Kx)
  }
  
  # Select optimal lambda values
  optBIC_idx <- which.min(BIC_rho)
  lb_opt     <- BICmat_l[optBIC_idx, 1]
  lrho_opt   <- BICmat_l[optBIC_idx, 2]
  
  return(list(
    lb   = lb_opt, 
    lrho = lrho_opt
  ))
}
