pen2SLS_est <- function(y, x, W, gpy, gpx, K0, Ky, Kx, lb, lrho) {
  
  # Dimensions
  n  <- nrow(y)
  py <- ncol(y)
  px <- ncol(x)
  
  y0 <- y
  y  <- c(t(y))  # Vectorize response
  
  # Create B-spline bases
  bs_basis0 <- create.bspline.basis(rangeval = c(gpy[1], gpy[py]), nbasis = K0)
  bs_basisy <- create.bspline.basis(rangeval = c(gpy[1], gpy[py]), nbasis = Ky)
  bs_basisx <- create.bspline.basis(rangeval = c(gpx[1], gpx[px]), nbasis = Kx)
  
  # Evaluate basis functions
  evalbase0 <- eval.basis(gpy, bs_basis0)
  evalbasey <- eval.basis(gpy, bs_basisy)
  evalbasex <- eval.basis(gpx, bs_basisx)
  
  # Grid spacing
  diff.x <- gpx[2] - gpx[1]
  diff.y <- gpy[2] - gpy[1]
  
  # Penalty matrices
  pm.b0  <- matrix(0, K0, K0)
  pm.b.t0 <- bsplinepen(bs_basisy, Lfdobj = 0)
  pm.b.s0 <- bsplinepen(bs_basisx, Lfdobj = 0)
  pm.b.t2 <- bsplinepen(bs_basisy, Lfdobj = 2)
  pm.b.s2 <- bsplinepen(bs_basisx, Lfdobj = 2)
  
  # Construct block diagonal penalty matrix
  pm.bt  <- kronecker(pm.b.s0, pm.b.t2)
  pm.bs  <- kronecker(pm.b.s2, pm.b.t0)
  pm.btu <- kronecker(pm.b.t0, pm.b.t2)
  pm.but <- kronecker(pm.b.t2, pm.b.t0)
  
  pen.mat <- bdiag(pm.b0, lb * (pm.bt + pm.bs), lrho * (pm.btu + pm.but))
  
  # Construct y-related basis terms
  wy  <- W %*% y0
  y.p <- (wy %*% evalbasey) * diff.y
  y.k <- kronecker(y.p, evalbasey)
  
  # Construct x-related basis terms
  x.p   <- (x %*% evalbasex) * diff.x
  x0.k  <- kronecker(rep(1, n), evalbase0)
  x1.k  <- kronecker(x.p, evalbasey)
  x.k   <- cbind(x0.k, x1.k)
  
  # Instruments
  z1 <- W %*% x
  z2 <- W %*% z1
  z1.p <- (z1 %*% evalbasex) * diff.x
  z2.p <- (z2 %*% evalbasex) * diff.x
  
  z1_expanded <- kronecker(z1.p, evalbasey)
  z2_expanded <- kronecker(z2.p, evalbasey)
  zstar       <- cbind(x.k, z1_expanded, z2_expanded)
  Pi          <- cbind(x.k, y.k)
  
  # Project Pi using instrument zstar
  zstar_sparse <- as(zstar, "sparseMatrix")
  Pi_sparse    <- as(Pi, "sparseMatrix")
  
  ztz     <- ginv(as.matrix(t(zstar_sparse) %*% zstar_sparse))
  zpi     <- t(zstar_sparse) %*% Pi_sparse
  Pi.hat  <- zstar_sparse %*% (ztz %*% zpi)
  
  # Penalized 2SLS estimation
  temp1     <- t(Pi.hat) %*% Pi_sparse + pen.mat
  temp2     <- t(Pi.hat) %*% y
  bhat_mat  <- ginv(as.matrix(temp1)) %*% as.matrix(temp2)
  
  # Extract coefficient components
  b0_mat <- bhat_mat[1:ncol(x0.k)]
  b_mat  <- bhat_mat[(ncol(x0.k) + 1):ncol(x.k)]
  r_mat  <- bhat_mat[(ncol(x.k) + 1):length(bhat_mat)]
  
  # Evaluate fitted components
  b0hat   <- evalbase0 %*% b0_mat
  bhat    <- evalbasey %*% matrix(b_mat, nrow = Ky, ncol = Kx) %*% t(evalbasex)
  rhohat  <- evalbasey %*% matrix(r_mat, nrow = Ky, ncol = Ky) %*% t(evalbasey)
  
  # Fitted values and residuals
  yhat    <- Pi %*% bhat_mat
  yhat    <- matrix(yhat, ncol = py, nrow = n, byrow = TRUE)
  resids  <- y0 - yhat
  
  return(list(
    b0hat         = b0hat,
    bhat          = bhat,
    rhohat        = rhohat,
    b0_mat        = b0_mat,
    b_mat         = b_mat,
    r_mat         = r_mat,
    fitted.values = yhat,
    residuals     = resids
  ))
}
