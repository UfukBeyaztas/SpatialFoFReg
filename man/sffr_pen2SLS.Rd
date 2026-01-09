\name{sffr_pen2SLS}
\alias{sffr_pen2SLS}
\title{Penalised Spatial Two–Stage Least Squares for SFoFR}
\description{
Fits the penalised spatial function-on-function regression (SFoFR) model
via the two–stage least–squares (Pen2SLS) estimator introduced by
Beyaztas, Shang and Sezer (2025). It selects optimal smoothing parameters, estimates regression and spatial
autocorrelation surfaces, and (optionally) builds percentile bootstrap confidence bands.
}
\usage{
sffr_pen2SLS(
  y, x, W, gpy, gpx,
  K0, Ky, Kx,
  lam_cands,
  boot       = FALSE,
  nboot      = NULL,
  percentile = NULL
)
}
\arguments{
  \item{y}{\code{n x length(gpy)} matrix of functional responses
           evaluated on grid \code{gpy}.}

  \item{x}{\code{n x length(gpx)} matrix of functional predictors
           evaluated on grid \code{gpx}.}

  \item{W}{\code{n x n} row–normalised spatial weight matrix,
           typically inverse–distance.}

  \item{gpy}{Numeric vector of response evaluation points
             \eqn{t \in [0,1]}.}

  \item{gpx}{Numeric vector of predictor evaluation points
             \eqn{s \in [0,1]}.}

  \item{K0}{Integer; number of basis functions for the intercept
            \eqn{\beta_0(t)}.}

  \item{Ky}{Integer; number of basis functions in the \emph{response}
            direction of the bivariate surfaces \eqn{\rho(t,u)} and
            \eqn{\beta(t,s)}.}

  \item{Kx}{Integer; number of basis functions in the \emph{predictor}
            direction of the regression surface \eqn{\beta(t,s)}.}

  \item{lam_cands}{Two–column matrix or data frame whose rows contain
            candidate smoothing pairs
            \eqn{(\lambda_\rho,\;\lambda_\beta)} to be ranked by BIC.}

  \item{boot}{Logical; if \code{TRUE} percentile bootstrap confidence
              intervals are produced.}

  \item{nboot}{Number of bootstrap resamples.  Required when
               \code{boot = TRUE}.}

  \item{percentile}{Desired CI nominal width in percent
                    (e.g., \code{95}).  Required when
                    \code{boot = TRUE}.}
}
\details{
The estimator minimises the penalised objective

\deqn{
\Vert Z^* \bigl\{\operatorname{vec}(Y) - \Pi\theta\bigr\}\Vert^2
\;+\;
\frac{1}{2}\lambda_\rho P(\rho)
\;+\;
\frac{1}{2}\lambda_\beta P(\beta),
}{%
||Z^*{vec}(Y)-\Pi\theta||^2 + (1/2)\lambda_\rho P(\rho) +
(1/2)\lambda_\beta P(\beta),
}

where \eqn{\theta = (\operatorname{vec}\rho,\operatorname{vec}\beta)} are
tensor-product B-spline coefficients,
\eqn{Z^*} is the projection onto instrumental variables, and
\eqn{P(\cdot)} are Kronecker-sum quadratic roughness penalties in both
surface directions.  Candidate smoothing pairs are scored by the Bayesian
Information Criterion

\deqn{
\mathrm{BIC} = -2\log\mathcal{L} + \omega \log n,
}

with log-likelihood based on squared residuals and
\eqn{\omega} equal to the effective degrees of freedom.

If \code{boot = TRUE}, residuals are centred, resampled, and the entire
estimation procedure is repeated \code{nboot} times.  Lower and upper
percentile bounds are then extracted for
\eqn{\beta(t,s)}, \eqn{\rho(t,u)}, and
\eqn{\widehat{Y}_i(t)}.
}
\value{
A named \code{list}:

\describe{
  \item{b0hat}{Estimated intercept curve \eqn{\widehat{\beta}_0(t)}.}
  \item{bhat}{Matrix of \eqn{\widehat{\beta}(t,s)} values.}
  \item{rhohat}{Matrix of \eqn{\widehat{\rho}(t,u)} values.}
  \item{b0_mat, b_mat, r_mat}{Raw coefficient matrices of B-spline
        basis weights for \eqn{\beta_0}, \eqn{\beta}, and \eqn{\rho}.}
  \item{fitted.values}{\eqn{\widehat{Y}_i(t)} matrix.}
  \item{residuals}{\eqn{\widehat{\varepsilon}_i(t)} matrix.}
  \item{CI_bhat}{Two-element list with lower/upper percentile surfaces
        (NULL unless \code{boot = TRUE}).}
  \item{CI_rhohat}{Analogous list for \eqn{\rho}.}
  \item{CIy}{Percentile bands for the fitted responses.}
  \item{gpy, gpx, K0, Ky, Kx}{Returned for convenience.}
}
}
\author{
Eylul Fidan, Ufuk Beyaztas, and Soutir Bandyopadhyay
}
\examples{
## Not run: ---------------------------------------------------------------
## 1. simulate data
# sim <- sff_dgp(n = 100, rf = 0.7)
## 2. candidate smoothing grid (four pairs)
# lam <- list(lb = c(10^{-3}, 10^{-2}, 10^{-1}),
#                       lrho = c(10^{-3}, 10^{-2}, 10^{-1}))
## 3. fit model without bootstrap
# fit <- sffr_pen2SLS(
#   y          = sim$Y,
#   x          = sim$X,
#   W          = sim$W,
#   gpy        = seq(0, 1, length.out = 101),
#   gpx        = seq(0, 1, length.out = 101),
#   K0         = 10,
#   Ky         = 10,
#   Kx         = 10,
#  lam_cands  = lam,
#  boot       = FALSE
# )
## End(Not run)
}
