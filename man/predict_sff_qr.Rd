\name{predict_sff_qr}
\alias{predict_sff_qr}
\title{
Prediction for Spatial Function-on-Function Quantile IV Model
}
\description{
Generates fitted values from a previously estimated spatial function-on-function quantile instrumental-variable model (via \code{\link{sff_qr}}) for new functional predictor data.
}
\usage{
predict_sff_qr(object, xnew, Wnew)
}
\arguments{
  \item{object}{An object returned by \code{\link{sff_qr}} containing the estimated surfaces \code{b0hat}, \code{bhat}, \code{rhohat} and the grids \code{grid_x}, \code{grid_y}.}
  \item{xnew}{A numeric \eqn{n_{\rm new} \times p_x} matrix of functional predictor observations at grid points \code{grid_x}.}
  \item{Wnew}{A numeric \eqn{n_{\rm new} \times n_{\rm new}} spatial weight matrix (row-normalised) for the new observations, used for the spatial feedback component.}
}
\details{
This function computes predicted values \eqn{\widehat Y_i(t_m)} for new data as follows:

1. It uses the estimated intercept surface \eqn{\widehat \beta_0(t)}, the estimated coefficient surface \eqn{\widehat \beta(t,s)} and the new predictor matrix \eqn{X_{\rm new}} to compute the initial fitted surface
\deqn{
\widehat Y^{(0)}_i(t_m)
= \widehat \beta_0(t_m)
+ \int_{s_\ell} \widehat \beta(t_m,s_\ell)\;X_{i,\rm new}(s_\ell)\;\Delta s.
}
2. It then incorporates the spatial feedback component using the estimated feedback surface \eqn{\widehat \rho(u,t)} and the matrix \eqn{W_{\rm new}}:
\deqn{
\widehat Y^{(k+1)}_i(t_m)
= \widehat \beta_0(t_m)
+ \int_{s_\ell} \widehat \beta(t_m,s_\ell)\;X_{i,\rm new}(s_\ell)\;\Delta s
+ \sum_j W_{{\rm new},i j} \int_{t_{m'}} \widehat Y^{(k)}_j(t_{m'})\;\widehat \rho(t_{m'},t_m)\;\Delta t.
}
The algorithm iterates until convergence (maximum absolute change \(<0.001\)) or until a maximum number of iterations (1000) is reached.
}
\value{
A numeric matrix of dimension \eqn{n_{\rm new} \times p_y} giving predicted functional responses \eqn{\widehat Y_i(t_m)} on the grid \code{grid_y}.
}
\author{
Eylul Fidan, Ufuk Beyaztas & Soutir Bandyopadhyay
}
\examples{
## Not run: ---------------------------------------------------------------
## 1. Fit a model on small simulated data
# train <- sff_dgp(n = 500, rf = 0.5)
# fit <- sff_qr(
#   y   = train$Y,
#   x   = train$X,
#   W   = train$W,
#   gpy = seq(0, 1, length.out = 101),
#   gpx = seq(0, 1, length.out = 101)
# )
#
## 2. Simulate NEW covariates and a compatible weight matrix
# test <- sff_dgp(n = 1000, rf = 0.5)  ## we keep only X and W
# pred <- predict_sff_qr(fit, xnew = test$X, Wnew = test$W)
## End(Not run)
}
