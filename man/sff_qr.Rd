\name{sff_qr}
\alias{sff_qr}
\title{
Spatial Function-on-Function Quantile Instrumental-Variable Estimation
}
\description{
Performs penalised two-stage spatial function-on-function quantile regression with tensor-product B-splines and smoothing-parameter selection.
}
\usage{
sff_qr(
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
  factr  = 1e7),
  verbose = TRUE,
  BIC = FALSE
)
}
\arguments{
  \item{y}{A numeric \eqn{n \times p_y} matrix of functional responses \eqn{Y_i(t_m)}, with \eqn{n} observations and \eqn{p_y} grid-points at times given by \code{gpy}.}
  \item{x}{A numeric \eqn{n \times p_x} matrix of functional predictors \eqn{X_i(s_\ell)}, with \eqn{n} observations and \eqn{p_x} grid-points at times given by \code{gpx}.}
  \item{W}{A numeric \eqn{n \times n} spatial weight matrix (row-normalised) capturing spatial autoregression of the functional response.}
  \item{gpx}{A numeric vector of length \eqn{p_x} giving the grid points \eqn{s_\ell} at which \code{x} is observed.}
  \item{gpy}{A numeric vector of length \eqn{p_y} giving the grid points \eqn{t_m} at which \code{y} is observed.}
  \item{K0}{Integer: number of basis functions for the intercept surface \eqn{\beta_0(t)}. Default is 10.}
  \item{Kx}{Integer: number of basis functions in the \eqn{s}-direction for the tensor-product basis for \eqn{\beta(t,s)}. Default is 10.}
  \item{Ky}{Integer: number of basis functions in the \eqn{t}-direction for the tensor-product basis for \eqn{\beta(t,s)} and for \eqn{\rho(t,u)}. Default is 10.}
  \item{tau}{Quantile level \eqn{\tau\in(0,1)} for the quantile regression. Default is 0.5 (median).}
  \item{lam_cands}{A list with components \code{lb} and \code{lrho} giving candidate values of smoothing parameters \eqn{\lambda_{\beta}} and \eqn{\lambda_{\rho}}. If non-NULL then smoothing-parameter selection via BIC is activated.}
  \item{lb}{Smoothing parameter \eqn{\lambda_{\beta}} for the roughness penalty on \eqn{\beta(t,s)}. Ignored if \code{lam_cands} is provided.}
  \item{lrho}{Smoothing parameter \eqn{\lambda_{\rho}} for the roughness penalty on \eqn{\rho(t,u)}. Ignored if \code{lam_cands} is provided.}
  \item{alpha}{Smoothing parameter for the smoothed pinball (‘‘Moreau/Huberised’’) loss approximation. Smaller \eqn{\alpha} ≈ closer to exact check-loss. Default is 1e-2.}
  \item{ridge_eps}{Small ridge parameter added to the penalty matrix to ensure positive-definiteness / numeric stability. Default is 1e-8.}
  \item{osqp_opts}{List of control parameters passed to the \code{optim} solver (using the \code{"L-BFGS-B"} method) in the second stage. Key components include: \code{maxit} (maximum number of iterations), \code{factr} (convergence tolerance control), and \code{trace} (verbosity level)}
  \item{verbose}{Logical. If \code{TRUE}, progress messages are printed from Stage 1 and Stage 2 of the algorithm.}
  \item{BIC}{Logical. If \code{TRUE} and \code{lam_cands} is non-NULL, the smoothing parameters are selected via the BIC criterion; if \code{FALSE}, the user-supplied \code{lb} and \code{lrho} are used directly.}
}
\details{
This function implements the two-stage instrumental-variable quantile regression for spatial function-on-function models.
In Stage 1, a quantile regression of the spatially-lagged functional response \eqn{W\,Y_i(t_m)} on basis-expanded instruments (including \eqn{(X_i(s),\,W\,X_i(s),\,W^2\,X_i(s))}) is performed for each grid-point combination via a call to \code{rq}. This produces fitted values \eqn{\widehat{W\,Y}_i(t_m)}.
In Stage 2, the design matrix is constructed by combining the basis expansions for the intercept surface, the tensor-product basis for \eqn{\beta(t,s)}, and a basis representation of the spatial feedback surface \eqn{\rho(t,u)}. A penalised smoothed-quantile objective is then solved
\deqn{
\min_{\beta,\rho} \sum_{i,m} \eta_{\tau}\!\bigl(Y_i(t_m) - \widehat Y_i(t_m)\bigr)
\;+\;\tfrac12\,\beta^\top P\,\beta,
}
where \eqn{\eta_\tau(u)} is the check-loss (approximated by a smoothed version) and \eqn{P} is the block-diagonal penalty matrix with roughness penalties controlled by \eqn{\lambda_{\beta}} and \eqn{\lambda_{\rho}}.
If \code{lam_cands} is provided, a grid-search over \eqn{(\lambda_{\beta},\lambda_{\rho})} is performed and the optimal pair minimises the BIC criterion:
\deqn{
\mathrm{BIC}(\lambda_{\rho},\lambda_{\beta})
= \log\!\Bigl(\tfrac1{nM}\sum_{i,m}\eta_\tau(Y_i(t_m) - \widehat Y_i(t_m))\Bigr)
\;+\;\frac{\log(nM)\;\mathrm{df}(\lambda_{\rho},\lambda_{\beta})}{\,nM\,},
}
where \eqn{n} is the number of individuals, \eqn{M} the number of grid-points, and \eqn{\mathrm{df}(\cdot)} denotes the effective degrees of freedom (approximated here by the number of free basis coefficients). The selected \eqn{\widehat\lambda_{\beta}} and \eqn{\widehat\lambda_{\rho}} are then used for the final fit reported in the surfaces component of the result.
}
\value{
A list with the following components:
\item{surfaces}{A list containing:
  \itemize{
    \item \code{grid_x} — the vector of grid locations \eqn{s_\ell} (same as input \code{gpx}),
    \item \code{grid_y} — the vector of grid locations \eqn{t_m} (same as input \code{gpy}),
    \item \code{b0hat} — the estimated intercept surface \eqn{\widehat\beta_0(t)}, as a numeric vector of length \eqn{p_y},
    \item \code{bhat} — the estimated coefficient surface \eqn{\widehat\beta(t,s)}, as a numeric matrix of dimension \eqn{p_y \times p_x},
    \item \code{rhohat} — the estimated spatial-feedback surface \eqn{\widehat\rho(t,u)}, as a numeric matrix of dimension \eqn{p_y \times p_y}.
  }
}
\item{fitted}{A numeric \eqn{n \times p_y} matrix of fitted values \eqn{\widehat Y_i(t_m)}.}
\item{resid}{A numeric \eqn{n \times p_y} matrix of residuals \eqn{Y_i(t_m) - \widehat Y_i(t_m)}.}
}
\author{
Eylul Fidan, Ufuk Beyaztas & Soutir Bandyopadhyay
}
\examples{
## Not run: ---------------------------------------------------------------
## simulate data
# sim <- sff_dgp(n = 250, rf = 0.7)
## fit model
# fit <- sff_qr(
#   y          = sim$Y,
#   x          = sim$X,
#   W          = sim$W,
#   gpy        = seq(0, 1, length.out = 101),
#   gpx        = seq(0, 1, length.out = 101)
# )
## End(Not run)
}
