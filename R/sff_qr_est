sff_qr_est <- function(
    y, x, W,
    gpx, gpy,
    K0, Kx, Ky,
    tau,
    lb,
    lrho,
    alpha,
    ridge_eps,
    osqp_opts,
    verbose
) {

  n  <- nrow(y); py <- ncol(y); px <- ncol(x)

  diff.x <- gpx[2] - gpx[1]
  diff.y <- gpy[2] - gpy[1]

  bs0 <- create.bspline.basis(rangeval = c(gpy[1], gpy[py]), nbasis = K0)
  bsy <- create.bspline.basis(rangeval = c(gpy[1], gpy[py]), nbasis = Ky)
  bsx <- create.bspline.basis(rangeval = c(gpx[1], gpx[px]), nbasis = Kx)

  E0 <- Matrix(eval.basis(gpy, bs0), sparse = TRUE)
  EY <- Matrix(eval.basis(gpy, bsy), sparse = TRUE)
  EX <- Matrix(eval.basis(gpx, bsx), sparse = TRUE)

  # Stage 1: instruments in basis
  y_vec <- c(t(y))
  x_p  <- Matrix(x %*% EX, sparse = TRUE) * diff.x
  ones_n <- Matrix(rep(1, n), ncol = 1, sparse = TRUE)

  x0_k <- kronecker(ones_n, E0)
  x1_k <- kronecker(x_p, EY)
  Xk   <- cbind(x0_k, x1_k)

  z1   <- Matrix(W %*% x, sparse = TRUE)
  z2   <- Matrix(W %*% z1, sparse = TRUE)
  z1_p <- Matrix(z1 %*% EX, sparse = TRUE) * diff.x
  z2_p <- Matrix(z2 %*% EX, sparse = TRUE) * diff.x

  z1_exp <- kronecker(z1_p, EY)
  z2_exp <- kronecker(z2_p, EY)

  wy     <- W %*% y
  wy_vec <- c(t(wy))

  Zstar <- cbind(Xk, z1_exp, z2_exp)
  if (verbose) cat("[Stage-1] QR(WY ~ Xk + WX + W2X) ...\n")
  fit1 <- rq(as.numeric(wy_vec) ~ as.matrix(Zstar) - 1, tau = 0.5, method = "fn")
  # tau = 0.5 is for numerical stability at lower and upper tails....
  wyhat_vec <- as.numeric(fitted(fit1))
  wyhat     <- matrix(wyhat_vec, nrow = n, ncol = py, byrow = TRUE)

  yphat  <- Matrix(wyhat %*% EY, sparse = TRUE) * diff.y
  yk_hat <- kronecker(yphat, EY)

  # Second-stage design
  Pi_hat <- cbind(Xk, yk_hat)
  Pi_hat <- as(Pi_hat, "generalMatrix")

  # Penalty matrices
  pm.b0  <- Diagonal(K0, 0)
  pm.b.t0 <- as(Matrix(bsplinepen(bsy, Lfdobj = 0), sparse = TRUE), "generalMatrix")
  pm.b.s0 <- as(Matrix(bsplinepen(bsx, Lfdobj = 0), sparse = TRUE), "generalMatrix")
  pm.b.t2 <- as(Matrix(bsplinepen(bsy, Lfdobj = 2), sparse = TRUE), "generalMatrix")
  pm.b.s2 <- as(Matrix(bsplinepen(bsx, Lfdobj = 2), sparse = TRUE), "generalMatrix")

  P_beta_t <- kronecker(pm.b.s0, pm.b.t2)
  P_beta_s <- kronecker(pm.b.s2, pm.b.t0)
  P_beta   <- P_beta_t + P_beta_s

  P_rho_tu <- kronecker(pm.b.t0, pm.b.t2) + kronecker(pm.b.t2, pm.b.t0)

  p_x0 <- K0; p_x1 <- Kx * Ky; p_yk <- Ky * Ky
  stopifnot(ncol(Pi_hat) == (p_x0 + p_x1 + p_yk))

  P <- bdiag(pm.b0, lb * P_beta, lrho * P_rho_tu)
  P <- as((P + t(P))/2, "generalMatrix") + Diagonal(ncol(P), ridge_eps)

  # Stage 2:
  if (verbose) cat("[Stage-2] penalized smoothed-QR at tau =", tau, "...\n")
  obj_fn <- function(beta, y, Pi, P, tau, alpha) {
    u <- as.numeric(y - Pi %*% beta)
    smooth_part <- -alpha * plogis(u / alpha, log.p = TRUE)

    loss <- sum(tau * u + smooth_part)
    penalty <- 0.5 * as.numeric(crossprod(beta, P %*% beta))

    return(loss + penalty)
  }

  grad_fn <- function(beta, y, Pi, P, tau, alpha) {
    u <- as.numeric(y - Pi %*% beta)
    weight <- tau - plogis(-u / alpha)
    grad <- (P %*% beta) - crossprod(Pi, weight)

    return(as.vector(grad))
  }

  beta_init <- rep(0, ncol(Pi_hat))
  res_optim <- optim(
    par = beta_init,
    fn = obj_fn,
    gr = grad_fn,
    y = as.matrix(y_vec),
    Pi = Pi_hat,
    P = P,
    tau = tau,
    alpha = alpha,
    method = "L-BFGS-B",
    control = osqp_opts
  )

  bhat_all <- res_optim$par

  b0_idx <- seq_len(p_x0)
  b_idx  <- (max(b0_idx) + 1):(p_x0 + p_x1)
  r_idx  <- (max(b_idx) + 1):(p_x0 + p_x1 + p_yk)

  b0_mat <- bhat_all[b0_idx]
  b_mat  <- bhat_all[b_idx]
  r_mat  <- bhat_all[r_idx]

  b0hat  <- as.matrix(E0 %*% b0_mat)
  bhat   <- as.matrix(EY %*% matrix(b_mat, nrow = Ky, ncol = Kx) %*% t(EX))
  rhohat <- as.matrix(EY %*% matrix(r_mat, nrow = Ky, ncol = Ky) %*% t(EY))

  yhat_vec <- as.numeric(Pi_hat %*% bhat_all)
  yhat     <- matrix(yhat_vec, nrow = n, ncol = py, byrow = TRUE)
  resid    <- y - yhat

  list(
    tau = tau, alpha = alpha, lb = lb, lrho = lrho,
    dims = list(K0 = K0, Kx = Kx, Ky = Ky, n = n, py = py, px = px),
    coef = list(
      b0  = b0_mat,
      b   = matrix(b_mat, nrow = Ky, ncol = Kx),
      rho = matrix(r_mat, nrow = Ky, ncol = Ky)
    ),
    surfaces = list(
      grid_x = gpx, grid_y = gpy,
      b0hat = b0hat,
      bhat  = bhat,
      rhohat = rhohat
    ),
    basis = list(E0 = E0, EY = EY, EX = EX, bs0 = bs0, bsy = bsy, bsx = bsx),
    design = list(Pi_hat = Pi_hat),
    fitted = yhat, resid = resid
  )
}
