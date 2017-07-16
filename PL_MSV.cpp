/***
 Multivariate Stochastic Volatility Model using Particle Learning Methods
Model:
y_{jt} exp(x_{jt}/2) \epsilon^y_{jt}
x_{jt} = \alpha_{j} + \beta_{j} x_{j,t-1} + \epsilon^x_{jt}
where,
\epsilon^y ~ MVN(0,I)
\epsilon^x ~ MVN(0, \Sigma)


Inputs:
y: J X T matrix of first differences of EEG signals at J channels over T periods (T is typically the length of the encoding period)
mu0, Sigma0: prior mean and covariance matrix for vector of log-vols at time 0.
N: number of particles (choosen to balance trade-off between speed and accuracy, N = 5000 typically)
m, v, p are mixtrure mean, variance, and probability vectors when linearizing y (page 3 eqn 11 and 12 in the manuscript)
prior: prior for alpha and beta


Outputs:
x: J X T matrix of latent log-volatilities
\beta: posterior estimates of persistences of volatilities (in theta variable)
\alpha: posterior estimates of intercept of volatility equation (in theta variable)
\rho: average off-diagonal entries of Sigma

Required R packages:
# required packages
Rcpp, RcppArmadillo, bayesm

*/

#include <RcppArmadillo.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp ;
using namespace arma;
Rcpp::Environment bayesm("package:bayesm");
Rcpp::Function rwishart = bayesm["rwishart"];


// Helper Functions
NumericVector quantile(NumericVector x, NumericVector q) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(q - 0.000000001)];
}

vec rmvnorm( mat mu , mat Sigma){
  mat C = chol(Sigma) ;
  int n = mu.n_elem ;
  mat x_old = rnorm(n,0,1) ;
  mat x = trans(C)*x_old  + mu;
  return(x);
}

double dmvnorm(mat x, mat mu, mat Sigma)
{
  double n = mu.n_elem/1.0;
  const double pi = std::atan(1.0)*4;
  double prob = -(1.0/2.0)*as_scalar(trans(x-mu)*inv(Sigma)*(x-mu)) - 1.0/2.0*  log(det(Sigma)) - n/2*log(2.0*pi)  ;
  return(prob);
}


vec Igamma(int n, double a, double b)
{
  return(1.0/rgamma(n,a, 1.0/b));
}


mat sample_x(mat y, mat x_pred, mat sigsq_x, mat m, mat v, IntegerVector index)
{
  int J = y.size();
  mat mu = zeros(J,1);
  mat Sigma = zeros(J,J);
  for(int j=0; j < J; j++)
  {
    mu(j,0) = m(index(j));
    Sigma(j,j) = v(index(j));
  }
  mu = inv(Sigma)*(y-mu) + inv(sigsq_x)*x_pred;
  Sigma = inv(inv(Sigma) + inv(sigsq_x));
  mu = Sigma*(mu);
  mat x = rmvnorm(mu,Sigma);
  return x;
}



// [[Rcpp::export]]
double sample_mu(double y, double mu_pred, double sigsq_mu, double sigsq_x)
{
  double mu =0;
  double sigma =0;
  mu = 1.0/sigsq_x*(y) + 1.0/sigsq_mu*(mu_pred);
  sigma = 1.0/(1.0/sigsq_x + 1.0/sigsq_mu);
  mu = sigma*mu;
  return(rnorm(1,mu,sqrt(sigma))(0));
}


// Compute predictive log density of y_t | x_t-1

// [[Rcpp::export]]
double predictive(mat y, mat x, mat sigsq_x, mat v, mat m, IntegerVector index)
{
  int J = y.size();
  mat mu = zeros(J,1);
  mat Sigma = zeros(J,J);
  for(int j=0; j < J; j++)
  {
    mu(j,0) = m(index(j));
    Sigma(j,j) = v(index(j));
  }
  double  hood = dmvnorm(y, x + mu, sigsq_x  + Sigma);
  return(hood);
}

double likelihood(double y, double mu, double x)
{
  double hood = Rf_dnorm4(y, mu, exp(x/2.0),1);
  return(hood);
}

// main function
// [[Rcpp::export]]
List PL_MSV(mat y, mat y_star, mat mu0, mat C0, mat theta0, mat Sigma_theta, int k, int n0, mat V0, int N,
            NumericVector p, NumericVector m, NumericVector v)
{
  int T = y.n_cols; // number of time points
  int J = y.n_rows; // number of electrodes

  // storing draws
  cube x_track = zeros(3,J,T);  // store 5%, 50% and 97.5% quantiles of x
  mat theta_x_track = zeros(3,2*J);  // ------------------- quantiles of theta_x = (alpha,beta)
  cube rho_track = zeros(J,J,3); // -------------------- average off-diagonal entries of Sigma

  cube x = zeros(J,1,N);   // current log-vol matrix
  cube x_temp = zeros(J,1,N);  // previous log-vol matrix
  mat theta_x = zeros(2*J,N); // current theta
  mat theta_x_temp = zeros(2*J,N); // previous theta
  cube Sigma_x =zeros(J,J,N); // curent Sigma_x

  cube s_theta_x_var = zeros(2*J,2*J,N);  // current sufficient statistics for theta_x (see Supporting Information on Page 9 for details)
  cube s_theta_x_var_temp = zeros(2*J,2*J,N); // previous sufficient statistics for theta_x
  mat s_theta_x_mean = mat(2*J,N); // current posterior mean of theta_x to be updated for each time point
  cube s_sigma_x = zeros(J,J,N); // current posterior covariance of theta_x
  mat s_theta_x_mean_temp = mat(2*J,N); // pervious posterior mean of theta_x
  cube s_sigma_x_temp = zeros(J,J,N); // previous posterior covariance of theta_x
  mat temp = zeros(2*J,1);
  NumericVector q = NumericVector::create(0.025,0.5,0.975); // quantiles to be stored
  NumericVector w(N);  // weights of samples
  IntegerVector choice = seq_len(N);  // mixture index (step 2 algorithm 1 page 3)
  IntegerVector indices(N);
  int i_new =0; // new mixture index
  mat x_pred = zeros(J,N); // predicted currrent x
  mat x_t_prev = eye(J,2*J);  // previous x at time t-1
  mat x_new = zeros(J,N); // new posterior x
  x_t_prev.submat(0,0,J-1,J-1)= eye(J,J);
  mat theta_x_mat = zeros(J,J);


  // store quantile values
  NumericMatrix x_part(J,N);
  NumericMatrix theta_x_part(2*J,N);
  NumericVector rho_part(N); // average correlation
  NumericVector Quantiles_x(3);
  NumericVector Quantiles_theta(3);
  NumericVector Quantiles_theta2(3);
  NumericVector Quantiles_rho(3);
  IntegerVector index(J);

  mat Sigma_x_inv = eye(J,J);
  mat Sigma_x_inv_kron = eye(2*J,2*J);
  mat epsilon_t = zeros(J,1);
  mat y_star_t = zeros(J,1); // log(y_t^2)

  // Initialize
  for(int i=0 ; i < N; i++)
  {
    x.slice(i) = rmvnorm(mu0,C0);
    s_theta_x_mean.col(i) = rmvnorm(theta0, Sigma_theta);
    Sigma_x.slice(i) = eye(J,J);
    s_theta_x_mean.col(i) = theta0;
    s_sigma_x.slice(i) = V0;
    s_theta_x_var.slice(i) = eye(2*J,2*J);
    theta_x.col(i) = theta0;
  }

  int n_p = p.size();
  IntegerVector ind = seq_len(n_p);
  IntegerMatrix lambda_t(J,N);


  for(int t = 1 ; t < T; t++)
  {
    y_star_t = y_star.col(t);
    x_temp = x;


    // store current sufficient statistics in temp. variables
    s_theta_x_mean_temp = s_theta_x_mean;
    s_sigma_x_temp = s_sigma_x;
    s_theta_x_var_temp = s_theta_x_var;


    // draw mixture indices (step 1 in algorithm 1 page 3 manuscript)
    for(int j=0 ; j < J ; j++)
    {
      lambda_t.row(j) =  Rcpp::RcppArmadillo::sample(ind,N,true,p) - 1;
    }

    // Resample (step 2 in algorithm 1 page 3 manuscript)
    for(int i=0; i < N; i++)
    {
      for(int j=0;j< J;j ++)
      {
        theta_x_mat(j,j)= theta_x(J+j,i);
      }
      index = lambda_t( _, i);
      x_pred.col(i) = theta_x.submat(0,i,J-1,i) + theta_x_mat*x_temp.slice(i);
      w(i) = predictive(y_star_t, x_pred.col(i), Sigma_x.slice(i), v,m, index);
    }
    w = exp(w-max(w));
    w = w/sum(w);
    indices =  Rcpp::RcppArmadillo::sample(choice, N, true, w);
    indices = indices-1;

    // propagate (step 3 in algorithm 1 page 3 manuscript)
    for(int i=0; i < N;i++)
    {
      i_new = indices(i);
      index = lambda_t( _, i_new);
      x.slice(i) = sample_x(y_star_t, x_pred.col(i_new), Sigma_x.slice(i_new),m,v,index);
      for(int j=0; j < J; j++)
      {
        x_part(j,i) = x(j,0,i);
      }
    }

    // update sufficient statistics and sample parameters (step 4 and 5 in algorithm 1 page 3 manuscript)
    for(int i =0; i < N; i++)
    {
      i_new = indices(i);
      Sigma_x_inv = inv(Sigma_x.slice(i_new));
      Sigma_x_inv_kron = kron(inv(Sigma_x_inv),eye(J,J));

      index = indices(i);
      x_t_prev.submat(0,J,J-1,2*J-1)= diagmat(vectorise(x_temp.slice(i_new)));
      s_theta_x_mean.col(i) = s_theta_x_mean_temp.col(i_new) + trans(x_t_prev)*Sigma_x_inv*x.slice(i); // update mean and co-variance of theta_x given new informaiton
      s_theta_x_var.slice(i) = s_theta_x_var_temp.slice(i_new) + trans(x_t_prev)*Sigma_x_inv*x_t_prev;
      temp =  rmvnorm(inv(s_theta_x_var.slice(i))*s_theta_x_mean.col(i), inv(s_theta_x_var.slice(i))); // sample theta_x


      for(int j=0;j < J;j++)
        {
          theta_x(j,i) = temp(j,0);
          if( (0 <=temp(j+J)) && (temp(j+J) <= 1))
            {
              theta_x(j+J,i) = temp(j+J);
            }
        }

      for(int j=0; j < J;j++)
      {
        theta_x_part(j,i) = theta_x(j,i);
        theta_x_part(j+J,i) = theta_x(j+J,i);
      }

      // sample Sigma
      epsilon_t = (x.slice(i) - x_t_prev*theta_x.col(i));
      s_sigma_x.slice(i) = s_sigma_x_temp.slice(i_new) + epsilon_t*trans(epsilon_t);
      List l = rwishart(n0 + t, inv(s_sigma_x.slice(i) + V0));
      NumericMatrix Sigma_temp = l["IW"];
      mat Sigma(Sigma_temp.begin(), J,J, false);
      Sigma_x.slice(i) = Sigma;
    }

    for(int j =0; j < J; j++)
    {
      Quantiles_x = quantile(x_part.row(j),q);
      for(int ii = 0; ii < 3; ii++)
      {
        x_track(ii,j,t) = Quantiles_x(ii);
      }
    }

    // storing posterior 2.5%,50% and 97.5%
    if(t == (T-1))
    {
      for(int j = 0; j < J; j++)
      {
        Quantiles_theta = quantile(theta_x_part.row(j),q);
        Quantiles_theta2 = quantile(theta_x_part.row(j+J),q);
        for(int jj = 0; jj < J; jj++)
        {
          for(int i =0; i < N; i++)
          {
            rho_part(i) = Sigma_x(j,jj,i);
          }
          Quantiles_rho = quantile(rho_part,q);
          for(int ii =0; ii < 3; ii++)
          {
            rho_track(j,jj,ii) = Quantiles_rho(ii);
            theta_x_track(ii,j) = Quantiles_theta(ii);
            theta_x_track(ii,j+J) = Quantiles_theta2(ii);
          }
        }
      }
    }

    // print every 100 time points
    if(t%100 ==0)
    {
      cout << t << endl;
    }
  }

  return List::create( Named("x_track") = x_track, Named("theta_track") = theta_x_track, Named("rho_track") = rho_track);
}
