sffr_pen2SLS <- function(y, x, W, gpy, gpx, K0, Ky, Kx, lam_cands,
                         boot = FALSE, nboot = NULL, percentile = NULL) {
  
  # Step 1: Parameter selection
  par_model <- par_opt(y, x, W, gpy, gpx, K0, Ky, Kx, lam_cands)
  lb    <- par_model$lb
  lrho  <- par_model$lrho
  
  # Step 2: Final model estimation
  final_model <- pen2SLS_est(y, x, W, gpy, gpx, K0, Ky, Kx, lb, lrho)
  
  # Extract estimates
  b0hat   <- final_model$b0hat
  bhat    <- final_model$bhat
  rhohat  <- final_model$rhohat
  yhat    <- final_model$fitted.values
  resids  <- final_model$residuals
  b0_mat  <- final_model$b0_mat
  b_mat   <- final_model$b_mat
  r_mat   <- final_model$r_mat
  
  # Initialize confidence intervals
  Lbhat <- Ubhat <- Lrhohat <- Urhohat <- Lyhat <- Uyhat <- NULL
  
  # Step 3: Bootstrap confidence intervals
  if (boot) {
    cen_resids <- scale(resids, scale = FALSE)
    indx <- seq_len(nrow(y))
    
    boot_bhat   <- vector("list", nboot)
    boot_rhohat <- vector("list", nboot)
    yhatboot    <- vector("list", nboot)
    
    for (i in seq_len(nboot)) {
      boot_indx   <- sample(indx, nrow(y), replace = TRUE)
      boot_resids <- cen_resids[boot_indx, ]
      booty       <- y + boot_resids
      
      boot_model <- pen2SLS_est(booty, x, W, gpy, gpx, K0, Ky, Kx, lb, lrho)
      
      boot_bhat[[i]]   <- boot_model$bhat
      boot_rhohat[[i]] <- boot_model$rhohat
      yhatboot[[i]]    <- boot_model$fitted.values
    }
    
    # Compute confidence intervals
    alpha <- (100 - percentile) / 2
    Lbhat    <- calculate_percentile(boot_bhat, alpha)
    Ubhat    <- calculate_percentile(boot_bhat, 100 - alpha)
    Lrhohat  <- calculate_percentile(boot_rhohat, alpha)
    Urhohat  <- calculate_percentile(boot_rhohat, 100 - alpha)
    Lyhat    <- calculate_percentile(yhatboot, alpha)
    Uyhat    <- calculate_percentile(yhatboot, 100 - alpha)
  }
  
  # Return all outputs
  return(list(
    b0hat         = b0hat, 
    bhat          = bhat, 
    rhohat        = rhohat, 
    b0_mat        = b0_mat,
    b_mat         = b_mat,
    r_mat         = r_mat,
    fitted.values = yhat, 
    residuals     = resids,
    CI_bhat       = list(Lbhat, Ubhat), 
    CI_rhohat     = list(Lrhohat, Urhohat),
    CIy           = list(Lyhat, Uyhat), 
    gpy           = gpy, 
    gpx           = gpx,
    K0            = K0, 
    Ky            = Ky, 
    Kx            = Kx
  ))
}
