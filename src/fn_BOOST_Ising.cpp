#include <chrono>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

//
// TO-DO: I need to clean up the code and maybe add .h files
//        I need to figure out what (hyper)parameters should be user-specified.
//        Implement progress bar or option to monitor convergence/process
//        Pass `burn_prop` as argument to R function -> ceiling(n_iter * burn_prop)
//        I need to time the code in the main R function, instead of C++
//        I need to clean up the code and maybe add .h files
//        I need to figure out what (hyper)parameters should be user-specified.
//        Implement progress bar or option to monitor convergence/process
//
//
// prepare list of data information
// List info = List::create(
//   _["n"] = n,  // sample size
//   _["N"] = N,  // neighborhood
//   _["Q"] = Q   // number of states
// );
//
// prepare list of algorithm settings
// List settings = List::create(
//   _["n_iter"]      = n_iter,
//   _["burn"]        = burn,
//   _["tau"]         = tau,
//   _["theta_init"]  = theta_init,
//   _["omega0_init"] = omega0_init,
//   _["M"]           = M
// );
//
// list of results
// return List::create(
//   _["P"]         = P,                 // observed data
//   _["logpost_s"] = logpost_s,         // samples for log-posterior
//   _["theta_s"]   = theta_s,           // samples for interaction
//   _["omega0_s"]  = omega0_s,          // samples for first-order intensity
//   _["time"]      = duration.count(),  // execution time
//   _["accept"]    = accept ,           // acceptance probabilities
//   _["info"]      = info,              // data information
//   _["params"]    = params,            // hyperparameters
//   _["settings"]  = settings           // algorithm and model specifications
// );
//

arma::vec copy(arma::vec x);
arma::mat copy(arma::mat X);
double H(vec p, IntegerMatrix A, NumericVector omega, double theta);
arma::vec potts_c_omega(arma::vec p, IntegerMatrix A, NumericVector omega, double theta);
void update_theta(arma::vec p, IntegerMatrix A, double omega0, double &theta, double sigma_theta, int burn, int it, double tau, int M, int &try_theta, double &accept_theta);
void update_omega0(arma::vec p, IntegerMatrix A, double &omega0, double theta, double sigma_omega0, int burn, int it, double tau, int M, int &try_omega0, double &accept_omega0);
double logprior_theta(double theta, double sigma_theta);
double logprior_omega0(double omega0, double sigma_omega0);

// [[Rcpp::export(.BOOST_Ising_MCMC_cpp)]]
List BOOST_Ising_MCMC_cpp(
    arma::vec p, IntegerMatrix A,
    double mean_omega0, double sigma_omega0,
    double mean_theta, double sigma_theta,
    int n_iter, int burn,
    int M, double tau,
    double omega0_init, double theta_init
)
{
  // read data information
  int n = p.n_elem;    // sample size
  int N = A.ncol();    // neighborhood
  int Q = p.max() + 1; // number of states

  // algorithm settings
  int try_omega0 = 0, try_theta = 0;
  double accept_omega0 = 0, accept_theta = 0;

  // temporary variables
  double omega0 = omega0_init;
  double theta  = theta_init;

  // allocate storage for MCMC samples
  NumericVector omega0_s(n_iter + 1), theta_s(n_iter + 1);

  // set initial values
  omega0_s(0) = omega0;
  theta_s(0)  = theta;

  // start timer
  auto start = std::chrono::high_resolution_clock::now();

  // MCMC algorithm
  for (int it = 1; it < n_iter + 1; it++)
  {
    // update interaction parameter
    update_theta(p, A, omega0, theta, sigma_theta, burn, it, tau, M, try_theta, accept_theta);

    // update first-order intensity parameter
    update_omega0(p, A, omega0, theta, sigma_omega0, burn, it, tau, M, try_omega0, accept_omega0);

    // store updated values
    omega0_s(it) = omega0;
    theta_s(it)  = theta;
  }

  // end timer and compute execution time
  auto stop = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

  // prepare list of acceptance probabilities
  List accept = List::create(
    _["omega0"]  = accept_omega0 / try_omega0,
    _["theta"]   = accept_theta / try_theta
  );

  // list of results
  return List::create(
    _["omega0_s"] = omega0_s,     // samples for first-order intensity
    _["theta_s"]  = theta_s,      // samples for interaction
    _["accept"]   = accept,       // acceptance probabilities
    _["time"]     = time.count()  // execution time
  );

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////

/*
 * Create a copy of the given vector
 */
arma::vec copy(arma::vec x)
{
  int n = x.n_elem;

  arma::vec x_t(n);

  for (uword i = 0; i < n; i++)
  {
    x_t(i) = x(i);
  }

  return x_t;
}

/*
 * Create a copy of the given matrix
 */
arma::mat copy(arma::mat X)
{
  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat X_t(n, p);

  for (uword i = 0; i < n; i++)
  {
    for (uword j = 0; j < p; j++)
    {
      X_t(i, j) = X(i, j);
    }
  }

  return X_t;
}

/*
 * Generate a random integer
 */
int sample_int(int n, sugar::probs_t probs = R_NilValue)
{
  return as<int>(
    sample(n, 1, false, probs)
  );
}

////////////////////////////////////////////////////////////////////////////////
// Ising model functions
////////////////////////////////////////////////////////////////////////////////

/*
 * Compute the first term of the Hamiltonian, viewed as the weighted average
 *   of the numbers of spots with different spins
 */
double weighted_average(arma::vec p, NumericVector omega)
{
  int n = p.n_elem;
  int Q = omega.length();

  double total = 0;

  for (int i = 0; i < n; i++)
  {
    for (int q = 0; q < Q; q++)
    {
      total += omega(q) * (p(i) == q);
    }
  }

  return total;
}

/*
 * Find the number of spins that interact with the current vertex
 */
double interacting_spins(arma::vec p, IntegerMatrix A, int i)
{
  IntegerVector a    = A(i, _);
  double interacting = 0;

  for (int j = 0; j < a.length(); j++)
  {
    int ind = a(j) - 1;

    if (ind < 0)
    {
      break;
    }

    interacting += p(i) != p(ind);
  }

  return interacting;
}

/*
 * Find the number of neighboring spins set to a certain state
 */
double neighboring_spins(arma::vec p, IntegerMatrix A, int k, int i)
{
  IntegerVector a  = A(i, _);
  double neighbors = 0;

  for (int j = 0; j < a.length(); j++)
  {
    int ind = a(j) - 1;

    if (ind < 0)
    {
      break;
    }

    neighbors += p(ind) != k;
  }

  return neighbors;
}

/*
 * Compute the Hamiltonian of the 2-dimensional lattice
 */
double H(arma::vec p, IntegerMatrix A, NumericVector omega, double theta)
{
  int n = p.n_elem;
  int N = A.ncol();

  double energy = weighted_average(p, omega);

  for (int i = 0; i < n; i++)
  {
    energy += theta * interacting_spins(p, A, i);
  }

  return -energy;
}

////////////////////////////////////////////////////////////////////////////////
// MCMC algorithm functions
////////////////////////////////////////////////////////////////////////////////

/*
 * Generate an auxiliary variable with the current lattice state
 */
arma::vec potts_c_omega(arma::vec p, IntegerMatrix A, NumericVector omega, double theta)
{
  int n = p.n_elem;
  int N = A.ncol();
  int Q = omega.length();

  NumericVector probs_t(Q);

  for (int i = 0; i < n; i++)
  {
    for (int q = 0; q < Q; q++)
    {
      probs_t(q) = exp(omega(q) + theta * neighboring_spins(p, A, q, i));
    }

    probs_t = probs_t / sum(probs_t);
    p(i) = sample_int(Q, probs_t) - 1;
  }

  return p;
}

/*
 * Auxiliary function to update the interaction parameter
 */
void update_theta(arma::vec p, IntegerMatrix A, double omega0, double &theta, double sigma_theta, int burn, int it, double tau, int M, int &try_theta, double &accept_theta)
{
  double theta_t      = R::rnorm(theta, tau);
  arma::mat p_t       = copy(p);
  NumericVector omega = {omega0, 1};

  double u = R::runif(0, 1);

  for (int m = 0; m < M; m++)
  {
    p_t = potts_c_omega(p, A, omega, theta_t);
  }

  double R_0 = H(p, A, omega, theta) - H(p_t, A, omega, theta) +
    H(p_t, A, omega, theta_t) - H(p, A, omega, theta_t) +
    logprior_theta(theta_t, sigma_theta) - logprior_theta(theta, sigma_theta);

  if (it > burn)
  {
    try_theta++;
  }

  // update values
  if (log(u) < R_0)
  {
    theta = theta_t;

    if (it > burn)
    {
      accept_theta++;
    }
  }

  return;
}

/*
 * Auxiliary function to update the first-order intensity parameter
 */
void update_omega0(arma::vec p, IntegerMatrix A, double &omega0, double theta, double sigma_omega0, int burn, int it, double tau, int M, int &try_omega0, double &accept_omega0)
{
  double omega0_t       = R::rnorm(omega0, tau);
  arma::vec p_t         = copy(p);
  NumericVector omega   = {omega0, 1};
  NumericVector omega_t = {omega0_t, 1};

  double u = R::runif(0, 1);

  for (int m = 0; m < M; m++)
  {
    p_t = potts_c_omega(p, A, omega_t, theta);
  }

  double R_w = H(p, A, omega, theta) - H(p_t, A, omega, theta) +
    H(p_t, A, omega_t, theta) - H(p, A, omega_t, theta) +
    logprior_omega0(omega0_t, sigma_omega0) - logprior_omega0(omega0, sigma_omega0);

  if (it > burn)
  {
    try_omega0++;
  }

  // update values
  if (log(u) < R_w)
  {
    omega0 = omega0_t;

    if (it > burn)
    {
      accept_omega0++;
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Model formulation functions
////////////////////////////////////////////////////////////////////////////////

/*
 * Compute the log-prior of the interaction parameter
 */
double logprior_theta(double theta, double sigma_theta)
{
  return -(theta * theta) / (2 * sigma_theta * sigma_theta);
}

/*
 * Compute the log-prior of the first-order intensity parameter
 */
double logprior_omega0(double omega0, double sigma_omega0)
{
  return -(omega0 * omega0 - 2 * omega0 + 1) / (2 * sigma_omega0 * sigma_omega0);
}
