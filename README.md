SpatialFoFReg <img src="https://img.shields.io/badge/R-%3E=3.5.0-1f425f.svg" alt="R (>= 3.5.0)" align="right" height="20"/>

Spatial Function-on-Function Regression Models (SpatialFoFReg) provides a comprehensive suite of functions for implementing advanced regression methods where both the predictor and response are functional, and observations are spatially dependent.

The package implements the Penalized Spatial Two-Stage Least Squares (Pen2SLS) estimator and Spatial Function-on-Function Quantile Regression (SFoF-QR), allowing researchers to model complex spatial interactions in functional data under various error structures.

🚀 Key Features

Dual Estimation Frameworks: Implements both mean-based Pen2SLS regression and robust Quantile Instrumental-Variable (IV) estimation.

Spatial Autoregression: Explicitly models spatial feedback through an autoregressive operator equation.

Flexible Data Generation: A robust data-generating process (sff_dgp) supporting homoscedastic, heteroscedastic, and asymmetric error distributions.

Automated Smoothing: Built-in smoothing parameter selection using the Bayesian Information Criterion (BIC).

Inference Tools: Support for percentile bootstrap confidence bands for surfaces and fitted values.

🛠 Installation

You can install the development version of SpatialFoFReg from GitHub:

# install.packages("devtools")
library(devtools)
install_github("UfukBeyaztas/SpatialFoFReg")

📈 Core Functions
Function	Description
sffr_pen2SLS	

Fits the penalized SFoFR model via the two-stage least-squares estimator.

sff_qr	

Performs penalized two-stage spatial function-on-function quantile regression.

sff_dgp	

Generates synthetic functional data from an SFOFR process under multiple error structures.

predict_sffr2SLS	

Produces out-of-sample predictions for fitted Pen2SLS models using a fixed-point solver.

predict_sff_qr	

Generates fitted values for new data from a previously estimated quantile IV model.

🧪 Simulation Scenarios

The sff_dgp function allows for testing models under three distinct error mechanisms:

Case "1": Homoscedastic Gaussian errors with constant variance.

Case "2": Signal-dependent heteroscedastic Gaussian errors with upper-tail contamination.

Case "3": Asymmetric Laplace Distribution (ALD) errors, introducing heavy tails and asymmetry.

💻 Quick Start Example
1. Generate Spatially Correlated Functional Data

library(SpatialFoFReg)

# Simulate 250 spatial units under heteroscedastic errors (Case 2)
sim_data <- sff_dgp(n = 250, rf = 0.7, case = "2")

2. Fit a Spatial Quantile Regression Model

# Define evaluation grids
grid_pts <- seq(0, 1, length.out = 101) # [cite: 364]

# Fit the median (tau = 0.5) spatial quantile model
fit_qr <- sff_qr(
  y = sim_data$Y, 
  x = sim_data$X, 
  W = sim_data$W,
  gpy = grid_pts, 
  gpx = grid_pts,
  tau = 0.5,
  BIC = TRUE
)

3. Prediction for New Spatial Units

# Generate new observations
test_data <- sff_dgp(n = 100, rf = 0.7) # [cite: 68]

# Predict functional responses incorporating spatial feedback
preds <- predict_sff_qr(
  object = fit_qr, 
  xnew = test_data$X, 
  Wnew = test_data$W
)
