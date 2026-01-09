\name{sff_dgp}
\alias{sff_dgp}
\title{Simulate data from a Spatial Function-on-Function Regression model under multiple error structures}
\description{
Generates synthetic functional predictors and responses from the spatial
function-on-function regression (SFoFR) data-generating process.
The generator supports three distinct error structures, homoscedastic Gaussian,
signal-dependent heteroscedastic Gaussian with upper-tail contamination, and
asymmetric Laplace-to enable flexible simulation scenarios.
}
\usage{
sff_dgp(
  n,
  nphi     = 10,
  gpy      = NULL,
  gpx      = NULL,
  rf       = 0.9,
  sd.error = 0.01,
  tol      = 0.001,
  max_iter = 1000,
  case     = c("1","2","3")
)
}
\arguments{
  \item{n}{Number of spatial units (curves) to generate.}

  \item{nphi}{Number of sine \emph{and} cosine basis functions used to build
              each functional predictor.  Total latent scores generated are
              therefore \code{2 * nphi}.}

  \item{gpy}{Numeric vector of evaluation points for the response
              domain \eqn{t \in [0,1]}.
              Defaults to an equally–spaced grid of 101 points.}

  \item{gpx}{Numeric vector of evaluation points for the predictor
              domain \eqn{s \in [0,1]}.
              Defaults to an equally–spaced grid of 101 points.}

  \item{rf}{Scalar in \eqn{(0,1)} controlling the strength of spatial
            autocorrelation through the surface \eqn{\rho(t,u)}.
            Values closer to 1 yield stronger dependence.}

  \item{sd.error}{Scale parameter controlling the magnitude of the noise
                  component. Its interpretation depends on the chosen case.}

  \item{tol}{Absolute tolerance used in the fixed–point iteration that solves
             the spatial autoregressive operator equation (stopping rule on the
             sup–norm of successive iterates).}

  \item{max_iter}{Maximum number of fixed–point iterations.  Prevents
                   infinite looping when strong spatial feedback and small
                   \code{tol} interact.}

  \item{case}{
    Specifies the noise generation mechanism:
    \describe{
      \item{"1"}{Homoscedastic Gaussian errors:
        \eqn{\varepsilon_i(t) \sim \mathcal{N}(0, \sigma^2)} with constant variance.}

      \item{"2"}{Signal-dependent heteroscedastic Gaussian errors with upper-tail
        contamination. The errors are generated as
        \eqn{\varepsilon_i(t) = \sigma_i(t)\eta_i(t) + B_i(t)c_{\mathrm{out}}\sigma_i(t)},
        where \eqn{\eta_i(t) \stackrel{i.i.d.}{\sim} \mathcal{N}(0,1)} and
        \eqn{\sigma_i(t)} is a time-location specific scale that increases smoothly
        with the local signal magnitude. Letting \eqn{m_i(t)} denote the absolute
        signal component (e.g., \eqn{m_i(t) = |\mu_i(t)|} with \eqn{\mu_i(t)} the
        non-error part of \eqn{Y_i(t)}), the code first rescales
        \eqn{m_i(t)} to \eqn{[0,1]} via
        \deqn{
          m_i^{\star}(t) =
          \frac{m_i(t) - \min_{j,u} m_j(u)}{\max_{j,u} m_j(u) - \min_{j,u} m_j(u) + 10^{-8}},
        }
        and sets \eqn{\sigma_i(t) = \sigma_{\varepsilon}\{1 + a_{\sigma} m_i^{\star}(t)\}},
        so that \eqn{\sigma_i(t)} ranges roughly from \eqn{\sigma_{\varepsilon}} to
        \eqn{(1 + a_{\sigma})\sigma_{\varepsilon}}. Upper-tail contamination is
        introduced through a Bernoulli indicator
        \eqn{B_i(t) \sim \mathrm{Bernoulli}\{p_i(t)\}}, with
        \eqn{p_i(t) = p_0 + p_1 m_i^{\star}(t)} (capped at a maximum probability,
        e.g. \eqn{p_i(t) \le 0.3}), and a positive shift of size
        \eqn{c_{\mathrm{out}}\sigma_i(t)} whenever \eqn{B_i(t)=1}. This yields a
        conditional error distribution that is both heteroscedastic and strongly
        right–skewed in high–signal regions.}

      \item{"3"}{Asymmetric Laplace errors:
        \eqn{\varepsilon_i(t) \sim \mathrm{ALD}(0, \sigma, p=0.5)}, introducing
        heavy tails and asymmetry.}
    }
  }
}
\details{
The generator mimics the penalised SFoFR set-up:

\deqn{
  Y_i(t) = \sum_{j=1}^{n} w_{ij}\int_0^1 Y_j(u) \rho(t,u) du +
           \int_0^1 X_i(s) \beta(t,s) ds + \varepsilon_i(t),
}{%
Y_i(t)=\sum_j w_{ij}\int Y_j(u)\rho(t,u)du + \int X_i(s)\beta(t,s)ds + \varepsilon_i(t)
}

where

\itemize{
  \item \eqn{w_{ij}} are row-normalised inverse–distance weights,
  \item \eqn{X_i(s)} is built from Fourier scores
        \eqn{\xi_{ijk} \sim \mathcal{N}(0,1)} and
        damped basis functions
        \eqn{\phi_k^{\cos}(s)=(k^{-3/2})\sqrt{2}\cos(k\pi s)} and
        \eqn{\phi_k^{\sin}(s)=(k^{-3/2})\sqrt{2}\sin(k\pi s)},
  \item the regression surface is
        \eqn{\beta(t,s)=2+s+t+0.5\sin(2\pi s t)},
  \item the spatial autocorrelation surface is
        \eqn{\rho(t,u)=rf\,(1+ut)/(1+|u-t|)},
  \item the error structure is determined by \code{case}, as described above.
}

Given the contraction condition
\eqn{\|\rho\|_\infty < 1/\|W\|_\infty}, the Neumann series defining
\eqn{(\mathbb{I}-\mathcal{T})^{-1}} converges and the solution is obtained by
simple fixed-point iterations until the change is below \code{tol}.
}
\value{
A named \code{list} with components
\describe{
  \item{Y}{\code{n x length(gpy)} matrix of observed functional
           responses on the grid \code{gpy}.}

  \item{Y_true}{Same dimension as \code{Y}; noise–free latent responses
           before adding \eqn{\varepsilon_i(t)}.}

  \item{X}{\code{n x length(gpx)} matrix of functional predictors.}

  \item{W}{\code{n x n} row-normalised spatial weight matrix based on
           inverse distances.}

  \item{rho}{\code{length(gpy) x length(gpy)} matrix containing
             \eqn{\rho(t,u)} evaluated on the response grid.}

  \item{beta}{\code{length(gpx) x length(gpy)} matrix containing
              \eqn{\beta(t,s)} evaluated on the Cartesian product of the
              predictor and response grids.}
}
}
\author{
Eylul Fidan, Ufuk Beyaztas, and Soutir Bandyopadhyay
}
\examples{
## Not run: ---------------------------------------------------------------
## generate datasets under three scenarios
# dat1 <- sff_dgp(n = 100, rf = 0.6, case = "1")  # Homoscedastic Gaussian
# dat2 <- sff_dgp(n = 100, rf = 0.6, case = "2")  # Heteroscedastic + upper-tail contamination
# dat3 <- sff_dgp(n = 100, rf = 0.6, case = "3")  # Asymmetric Laplace
## End(Not run)
}
