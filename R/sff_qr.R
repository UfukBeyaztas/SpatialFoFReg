sff_qr <- function(
    y, x, W,
    gpx, gpy,
    K0 = 10, Kx = 10, Ky = 10,
    tau = 0.5,
    lam_cands = NULL,
    lb = 0.01,
    lrho = 0.01,
    alpha = 1e-2,
    ridge_eps = 1e-8,
    osqp_opts = list(
      trace  = 0,
      maxit = 8000,
      factr  = 1e7
    ),
    verbose = TRUE,
    BIC = FALSE){

  stopifnot(is.matrix(y), is.matrix(x), is.matrix(W))
  n  <- nrow(y); py <- ncol(y); px <- ncol(x)
  stopifnot(nrow(x) == n, nrow(W) == n, ncol(W) == n,
            ncol(y) == py, ncol(x) == px)

  if(!is.null(lam_cands))
    BIC = TRUE

  if(BIC == FALSE){
    final_model <- sff_qr_est(
      y = y, x = x, W = W,
      gpx = gpx, gpy = gpy,
      K0 = K0, Kx = Kx, Ky = Ky,
      tau = tau,
      lb = lb,
      lrho = lrho,
      alpha = alpha,
      ridge_eps = ridge_eps,
      osqp_opts = osqp_opts,
      verbose = verbose
    )
  }

  if(BIC == TRUE & !is.null(lam_cands)){
    sel <- select_smoothing_params(y = y, x = x, W = W,
                                   gpx = gpx, gpy = gpy,
                                   K0 = K0, Kx = Kx, Ky = Ky,
                                   tau = tau,
                                   alpha = alpha,
                                   ridge_eps = ridge_eps,
                                   osqp_opts = osqp_opts,
                                   lam_cands = lam_cands,
                                   verbose = verbose)

    final_model <- sff_qr_est(
      y = y, x = x, W = W,
      gpx = gpx, gpy = gpy,
      K0 = K0, Kx = Kx, Ky = Ky,
      tau = tau,
      lb   = sel$lb_opt,
      lrho = sel$lrho_opt,
      alpha = alpha,
      ridge_eps = ridge_eps,
      osqp_opts = osqp_opts,
      verbose = verbose)


  }


  list(
    surfaces = list(
      grid_x = gpx, grid_y = gpy,
      b0hat = final_model$surfaces$b0hat,
      bhat  = final_model$surfaces$bhat,
      rhohat = final_model$surfaces$rhohat
    ),
    fitted = final_model$fitted, resid = final_model$resid
  )
}
