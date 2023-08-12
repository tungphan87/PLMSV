Particle Learning for Multivariate Stochastic Volatility (MSV) Models

This is a C++ implementation of the particle learning algorithm for MSV models via Rcpp. Ones need the following R packages to run this code :Rcpp, RcppArmadillo, mvtnorm (install.packages(c('Rcpp', 'RcppArmadillo', 'mvtnorm')) to install).

Model: 

$$y_{jt} = exp(x_{jt}/2) \epsilon^y_{jt}$$

$$x_{jt} = \alpha_{j} + \beta_{j} x_{j,t-1} + \epsilon^x_{jt}$$
  
where, 

$$\epsilon^y \sim MVN(0,I)$$

$$\epsilon^x \sim MVN(0, \Sigma)$$

Inputs: 

      y: J X T matrix of first differences of EEG signals at J channels over T periods (T is typically the length of the encoding period),

      y_star: log(y^2 + 1.0e-5) mu0, C0: prior mean and covariance matrix for vector of log-vols at time 0, x_0,

      n0, V0: Inverse-Wishart prior for Sigma,

      theta0, Sigma_theta: prior for the vector theta = (\alpha,\beta),

      N: number of particles (choosen to balance trade-off between speed and accuracy, N = 5000 typically),

      m, v, p: mixtrure means, variances, and probability vectors when linearizing y (page 3 eqn 11 and 12 in the manuscript).

Outputs:

      x: J X T matrix of latent log-volatilities,

      \beta: posterior estimates of persistences,

      \rho: average off-diagonal entries of Sigma,

      \alpha: posterior estimates of intercept of volatility equation (in theta variable).


To Run code:
1. Install R (https://www.r-project.org/about.html) and C++ compiler (https://clang.llvm.org/get_started.html).
2. Install required R packages: install.packages(c('Rcpp', 'RcppArmadillo', 'mvtnorm')).
3. Source PL_MSV C++ code: sourceCpp('PL_MSV.cpp') and use PL_MSV just like a regular R function. 
4. Type >> Rscript simulation.R (on mac) or >> R --no-save simulation.R (on linux) to run an example file with simulated data. 
