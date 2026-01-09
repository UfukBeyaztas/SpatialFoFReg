sff_dgp <- function(n, nphi = 10, gpy = NULL, gpx = NULL,
                    rf = 0.9, sd.error = 0.01, tol = 0.001,
                    max_iter = 1000, case = c("1","2","3")) {

  # Default grids if not provided
  if (is.null(gpy)) gpy <- seq(0, 1, length.out = 101)
  if (is.null(gpx)) gpx <- seq(0, 1, length.out = 101)

  # Define cosine basis
  X.phi1 <- function(nphi, gpx) {
    phi <- matrix(0, nphi, length(gpx))
    for (j in 1:nphi) {
      phi[j, ] <- (j^-1.5) * sqrt(2) * cos(j * pi * gpx)
    }
    return(phi)
  }

  # Define sine basis
  X.phi2 <- function(nphi, gpx) {
    phi <- matrix(0, nphi, length(gpx))
    for (j in 1:nphi) {
      phi[j, ] <- (j^-1.5) * sqrt(2) * sin(j * pi * gpx)
    }
    return(phi)
  }

  # Generate one functional covariate realization
  rX.s <- function(nphi, gpx) {
    xsi <- rnorm(2 * nphi)
    phi_all <- rbind(X.phi1(nphi, gpx), X.phi2(nphi, gpx))
    X <- xsi * phi_all
    Xs <- colSums(X)
    return(Xs)
  }

  # Generate spatial weights matrix
  generate_W <- function(n) {
    wei <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          wei[i, j] <- 1 / (1 + abs(i - j))
        }
      }
    }
    W <- matrix(0, n, n)
    for (i in 1:n) {
      W[i, ] <- wei[i, ] / sum(wei[i, ])
    }
    return(W)
  }

  # Regression coefficient surface
  beta <- function(s, t) {
    2 + s + t + 0.5 * sin(2 * pi * s * t)
  }

  # Spatial autocorrelation surface
  rho <- function(u, t) {
    rf * (1 + u * t) / (1 + abs(u - t))
  }

  # Generate functional predictor X (n x length(gpx))
  X <- t(replicate(n, rX.s(nphi, gpx)))

  # Compute spatial weight matrix
  W <- generate_W(n)

  # Evaluate beta(s, t) over grids
  beta.ts <- outer(gpx, gpy, beta)

  # Compute X \beta(t)
  delta_gpx <- gpx[2] - gpx[1]
  Xbeta.t <- X %*% beta.ts * delta_gpx

  if(case == "1"){
    # Add Gaussian noise
    eps <- matrix(rnorm(n * length(gpy), mean = 0, sd = sd.error), n, length(gpy))
  }else if(case == "2"){
    n_t   <- length(gpy)
    n_obs <- n * n_t

    # Signal magnitude used for heteroscedasticity and contamination
    mu_abs <- abs(Xbeta.t)
    mu_min <- min(mu_abs)
    mu_max <- max(mu_abs)
    mu_std <- (mu_abs - mu_min) / (mu_max - mu_min + 1e-8)

    # Heteroscedastic base scale (bounded)
    a_scale <- 3.0
    sigma_mat <- sd.error * (1 + a_scale * mu_std)

    # Base symmetric heteroscedastic noise
    z0 <- matrix(rnorm(n_obs), nrow = n, ncol = n_t)
    base_noise <- sigma_mat * z0

    # X-dependent contamination probability for upper outliers
    p0 <- 0.05
    p1 <- 0.20
    p_mat <- p0 + p1 * mu_std
    p_mat[p_mat > 0.3] <- 0.3

    B <- matrix(
      rbinom(n_obs, size = 1, prob = as.vector(p_mat)),
      nrow = n, ncol = n_t
    )

    # Outlier shift: moderate, one-sided, proportional to local sigma
    c_out <- 5
    shift_mat <- c_out * sigma_mat

    # Final errors: heteroscedastic + upper-tail contamination
    eps <- base_noise + B * shift_mat
  }else if(case == "3"){
    n_obs <- n * length(gpy)

    # Generate errors from Asymmetric Laplace Distribution

    eps_matrix <- ald::rALD(n_obs,
                            mu = 0,
                            sigma = sd.error,
                            p = 0.5)

    eps <- matrix(eps_matrix, n, length(gpy))
  }


  G_t <- Xbeta.t + eps

  # Evaluate ??(u, t) over grid
  rho.ut <- outer(gpy, gpy, rho)

  # Spatial autoregressive operator function
  fn <- function(f) {
    as.matrix(W %*% f %*% rho.ut * (gpy[2] - gpy[1]))
  }

  # Fixed-point iteration to solve spatial functional model
  Y_0 <- G_t
  Y_1 <- G_t + fn(Y_0)
  diff <- max(abs(Y_0 - Y_1))
  Y_0 <- Y_1
  L <- 1

  while (diff > tol && L < max_iter) {
    Y_1 <- G_t + fn(Y_0)
    diff <- max(abs(Y_0 - Y_1))
    Y_0 <- Y_1
    L <- L + 1
  }

  # Final outputs
  Y <- Y_1
  Y_true <- Y - eps

  return(list(
    Y       = Y,
    Y_true  = Y_true,
    X       = X,
    W       = W,
    rho     = rho.ut,
    beta    = beta.ts
  ))
}
