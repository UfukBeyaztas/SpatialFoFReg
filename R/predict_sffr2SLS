predict_sffr2SLS <- function(object, xnew, Wnew) {
  
  # Extract model components
  gpx     <- object$gpx
  gpy     <- object$gpy
  K0      <- object$K0
  Ky      <- object$Ky
  Kx      <- object$Kx
  b0_mat  <- object$b0_mat
  b_mat   <- object$b_mat
  r_mat   <- object$r_mat
  
  py <- length(gpy)
  px <- length(gpx)
  n  <- nrow(xnew)
  
  # Create B-spline basis functions
  bs_basis0 <- create.bspline.basis(rangeval = c(gpy[1], gpy[py]), nbasis = K0)
  bs_basisy <- create.bspline.basis(rangeval = c(gpy[1], gpy[py]), nbasis = Ky)
  bs_basisx <- create.bspline.basis(rangeval = c(gpx[1], gpx[px]), nbasis = Kx)
  
  # Evaluate basis functions on grids
  evalbase0 <- eval.basis(gpy, bs_basis0)
  evalbasey <- eval.basis(gpy, bs_basisy)
  evalbasex <- eval.basis(gpx, bs_basisx)
  
  diff.x <- gpx[2] - gpx[1]
  diff.y <- gpy[2] - gpy[1]
  
  # Construct design matrices for new X
  x.p  <- (xnew %*% evalbasex) * diff.x
  x0.k <- kronecker(rep(1, n), evalbase0)
  x1.k <- kronecker(x.p, evalbasey)
  x.k  <- cbind(x0.k, x1.k)
  
  # Initial prediction
  y0     <- x.k %*% c(b0_mat, b_mat)
  y0     <- matrix(y0, nrow = n, byrow = TRUE)
  wy     <- Wnew %*% y0
  y.p    <- (wy %*% evalbasey) * diff.y
  y.k    <- kronecker(y.p, evalbasey)
  Pi     <- cbind(x.k, y.k)
  
  yhat_pred <- Pi %*% c(b0_mat, b_mat, r_mat)
  yhat_pred <- matrix(yhat_pred, nrow = n, ncol = py, byrow = TRUE)
  
  # Fixed-point iteration for final prediction
  yit        <- yhat_pred
  yhat.diff  <- 1e6
  niter      <- 1
  
  while (yhat.diff > 0.001 && niter < 1000) {
    wy_it     <- Wnew %*% yit
    y.p_it    <- (wy_it %*% evalbasey) * diff.y
    y.k_it    <- kronecker(y.p_it, evalbasey)
    Pi_it     <- cbind(x.k, y.k_it)
    yit1      <- Pi_it %*% c(b0_mat, b_mat, r_mat)
    yit2      <- matrix(yit1, nrow = n, ncol = py, byrow = TRUE)
    yhat.diff <- max(abs(yit - yit2))
    
    yit <- yit2
    niter <- niter + 1
  }
  
  yhat_pred <- yit2
  
  return(yhat_pred)
}
