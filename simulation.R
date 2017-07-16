# required packages
library(Rcpp)
library(RcppArmadillo)
library(bayesm)
library(mvtnorm)

# load cpp file 
sourceCpp('PL_MSV.cpp')
set.seed(20)

# Simulate data
J = 3
Time = 2000
alpha = rnorm(J, -0.0,0.2)
b = diag(runif(J,0.95,0.98))
x0 = matrix(0,1, J)

x = matrix(0, J, Time)
y = matrix(0, J, Time)

Sigma = rWishart(1, df = J + 3, diag(J)*0.05)[,,1]

x0 = rnorm(J,alpha/(1-diag(b)), sqrt(diag(Sigma)/(1-diag(b)^2)))
x0~ rnorm(J,0,1)


for(t in 1:Time)
{
  if(t==1)
  {
    x[,t]  = alpha + b%*%x0 + mvrnorm(1, rep(0, J), Sigma)
  }else
  {
    
    x[,t] = alpha + b%*%x[,t-1] +  mvrnorm(1, rep(0, J), Sigma)
  }
  y[,t] = rnorm(J, 0, exp(x[,t]/2))
}

# Priors
mu0 = matrix(0,J,1)
C0 = diag(J)*20
theta0 = matrix(0, 2*J,1)
theta0[(J+1):(2*J)] = 0.95
Sigma_theta = diag(2*J)*0.1
k = 1
n0 = 2*J+5
V0 = diag(J)
N = 10000

y_star = log(y^2+0.00001) 

p= c(0.00609, 0.04775, 0.13057, 0.20674, 0.22715, 0.18842, 0.12047, 0.05591, 0.01575,  0.00115)
m = c( 1.92677, 1.34744, 0.73504,0.02266, -0.85173, -1.97278,  -3.46788, -5.55246, -8.68384, -14.65000 )
v= c(0.11265, 0.17788, 0.26768, 0.40611, 0.62699,  0.98583, 1.57469, 2.54498 , 4.16591,7.33342 )

# Run code
result = PL_MSV(y,y_star, mu0, C0, theta0,Sigma_theta, k, n0,V0,N,p,m,v)
run.time = proc.time()-pt


# extract posterior draws
x_track = result$x_track
rho_track = result$rho_track
theta_track = result$theta_track

# plot true vs posterior volaitlity
channel = 3
plot(x[channel,], type ='l', col ='red')
lines(x_track[1,channel,], col = 'grey') # 2.5%
lines(x_track[2,channel,], col = 'blue') # 50%
lines(x_track[3,channel,], col = 'grey') # 97.5%

