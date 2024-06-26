#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // extract parameters
  double intercept = params["intercept"];
  double slope = params["slope"];
  double quadratic = params["quadratic"];
  double sigma = params["sigma"];
  
  // unpack data
  std::vector<double> Girth = Rcpp::as< std::vector<double> >(data["Girth"]);
  std::vector<double> Volume = Rcpp::as< std::vector<double> >(data["Volume"]);
  
  // calculate log-likelihood
  double loglikelihood = 0.0;
  for (unsigned int i = 0; i < Girth.size(); ++i) {
    double vol_predict = intercept + slope*Girth[i] + quadratic*Girth[i]*Girth[i];
    loglikelihood += R::dnorm(Volume[i], vol_predict, sigma, true);
  }
  
  // return in SEXP format
  return Rcpp::wrap(loglikelihood);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // extract parameters
  double intercept = params["intercept"];
  double slope = params["slope"];
  double quadratic = params["quadratic"];
  double sigma = params["sigma"];
  
  // calculate log-prior
  double logprior = R::dunif(intercept, 0, 100, true) +
    R::dunif(slope, -10, 10, true) +
    R::dunif(quadratic, 0, 1, true) +
    R::dlnorm(sigma, 0, 1, true);
  
  // return in SEXP format
  return Rcpp::wrap(logprior);
}


// NOTE: Do not edit this function name
// [[Rcpp::export]]  
SEXP create_xptr(std::string function_name) {  
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);  
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);  
  
  // NOTE: If your loglikelihood function is not called "loglike" please edit:
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  } 
  // NOTE: If your logprior function is not called "logprior" please edit:
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  } 
  
  stop("cpp function %i not found", function_name);
}
