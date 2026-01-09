select_smoothing_params <- function(
    y, x, W,
    gpx, gpy,
    K0, Kx, Ky,
    tau,
    alpha,
    ridge_eps,
    osqp_opts,
    lam_cands,
    verbose = TRUE
) {
  if (is.null(lam_cands$lb) || is.null(lam_cands$lrho)) {
    stop("Need candidate values for lam_cands$lb and lam_cands$lrho")
  }
  
  grid <- expand.grid(lb = lam_cands$lb, lrho = lam_cands$lrho)
  ncomb <- nrow(grid)
  bic_vals <- numeric(ncomb)
  
  for (j in seq_len(ncomb)) {
    lb_j   <- grid$lb[j]
    lrho_j <- grid$lrho[j]
    if (verbose) cat(sprintf("Trying lb = %g, lrho = %g ???\n", lb_j, lrho_j))
    
    fit_j <- sff_qr_est(
      y = y, x = x, W = W,
      gpx = gpx, gpy = gpy,
      K0 = K0, Kx = Kx, Ky = Ky,
      tau = tau,
      lb = lb_j, lrho = lrho_j,
      alpha = alpha,
      ridge_eps = ridge_eps,
      osqp_opts = osqp_opts,
      verbose = verbose
    )
    
    yhat_j <- fit_j$fitted
    resid_mat <- y - yhat_j
    u <- as.numeric(resid_mat)
    loss_j <- sum( (u >= 0) * tau * u + (u < 0) * (tau - 1) * u )
    n <- nrow(y)
    M <- ncol(y)
    df_j <- NA
    df_j <- length(fit_j$coef$b0) + length(fit_j$coef$b) + length(fit_j$coef$rho)
    bic_vals[j] <- log(loss_j/(n*M)) + (log(n*M) * df_j)/(n*M)
  }
  
  best_idx <- which.min(bic_vals)
  best_lb   <- grid$lb[best_idx]
  best_lrho <- grid$lrho[best_idx]
  
  if (verbose) {
    cat(sprintf("Selected lb = %g, lrho = %g (min BIC = %g)\n",
                best_lb, best_lrho, bic_vals[best_idx]))
  }
  
  return(list(lb_opt = best_lb, lrho_opt = best_lrho,
              bic_grid = data.frame(grid, bic = bic_vals)))
}
