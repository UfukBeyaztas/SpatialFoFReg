\name{predict_sffr2SLS}
\alias{predict_sffr2SLS}
\title{Out–of–sample prediction for Penalised Spatial FoFR models}
\description{
Given a fitted model returned by \code{\link{sffr_pen2SLS}}, this function
produces predicted functional responses at new spatial units whose functional
covariates and spatial weight matrix are supplied by the user.
A fixed–point solver enforces the spatial autoregressive feedback implicit in
the SFoFR model.
}
\usage{
predict_sffr2SLS(object, xnew, Wnew)
}
\arguments{
  \item{object}{An object of class \code{"sffr2SLS"}:
                the list produced by \code{\link{sffr_pen2SLS}}.
                At minimum it must contain
                \code{gpx}, \code{gpy},
                \code{K0}, \code{Ky}, \code{Kx},
                and the B–spline coefficient matrices
                \code{b0_mat}, \code{b_mat}, \code{r_mat}.}

  \item{xnew}{Numeric matrix of dimension
              \eqn{n_\mathrm{new} \times |\mathrm{gpx}|},
              holding the functional covariate for the \emph{new}
              spatial units, evaluated on the same predictor grid
              used during model fitting.}

  \item{Wnew}{Row–normalised \eqn{n_\mathrm{new} \times n_\mathrm{new}}
              spatial weight matrix that captures proximity among the
              \emph{new} units.  Its definition should mirror that of
              the training matrix (e.g.\ inverse distance, k–nearest
              neighbours, etc.).}
}
\details{
Let \eqn{\widehat{\beta}_0(t)}, \eqn{\widehat{\beta}(t,s)}, and
\eqn{\widehat{\rho}(t,u)} be the estimated surfaces stored in
\code{object}.  For each new unit i the algorithm first forms the
non–spatial regression prediction

\deqn{
  \widehat{G}_i(t) \;=\;
  \widehat{\beta}_0(t) +
  \int_{0}^{1} X_i(s)\,\widehat{\beta}(t,s)\,ds,
}

computed efficiently by pre-evaluated B–spline bases.
Spatial feedback is then introduced by iterating

\deqn{
  Y_i^{(\ell+1)}(t) \;=\;
  \widehat{G}_i(t) +
  \sum_{j=1}^{n_\mathrm{new}} w_{ij}
        \int_{0}^{1} Y_j^{(\ell)}(u)\,\widehat{\rho}(t,u)\,du,
}

until the sup–norm difference between successive curves falls below
\code{1e-3} or \code{1,000} iterations are reached.
Convergence is guaranteed when
\eqn{\|\widehat{\rho}\|_\infty < 1/\|Wnew\|_\infty}, a condition typically
satisfied by the fitted model if the training weight matrix met it
during estimation.
}
\value{
A numeric matrix of dimension
\eqn{n_\mathrm{new} \times |\mathrm{gpy}|}
containing the predicted functional responses evaluated on \code{gpy}.
Row \eqn{i} corresponds to the i-th row of \code{xnew}.
}
\author{
Eylul Fidan, Ufuk Beyaztas, and Soutir Bandyopadhyay
}
\note{
If the new weight matrix induces very strong dependence,
the fixed–point iterations may converge slowly.
Consider scaling \code{Wnew} to have
\eqn{\|Wnew\|_\infty \le 1} or relaxing the tolerance.
}
\seealso{
\code{\link{sff_dgp}} for simulated data generation;
\code{\link{sffr_pen2SLS}} for model fitting.
}
\examples{
## Not run: ---------------------------------------------------------------
## 1. Fit a model on small simulated data
# train <- sff_dgp(n = 500, rf = 0.5)
# lam   <- list(lb = c(10^{-3}, 10^{-2}, 10^{-1}),
#                       lrho = c(10^{-3}, 10^{-2}, 10^{-1}))

# fit <- sffr_pen2SLS(train$Y, train$X, train$W,
#                     gpy = seq(0, 1, length = 101),
#                     gpx = seq(0, 1, length = 101),
#                     K0 = 10, Ky = 10, Kx = 10,
#                     lam_cands = lam)

## 2. Simulate NEW covariates and a compatible weight matrix
# test <- sff_dgp(n = 1000, rf = 0.5)  ## we keep only X and W
# pred <- predict_sffr2SLS(fit, xnew = test$X, Wnew = test$W)
## End(Not run)
}
