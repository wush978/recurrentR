#include <Rcpp.h>

using namespace Rcpp;

//'@title Generate Homogeneous Poisson Process
//'
//'@param lambda Numeric value, the intensity.
//'@param T_0 Positive numeric value, the upper bound.
//'
//'@description Generate a homogeneous poisson process with constant intensity \code{lambda} 
//'in time interval \code{[0,T_0]}.
//'
//'@export
//[[Rcpp::export]]
SEXP gen_homo_poisson(double lambda, double T_0) {
	BEGIN_RCPP
	if (lambda < 0) throw std::invalid_argument("lambda > 0 does not hold!");
	std::vector<double> retval;
	double u(unif_rand());
	double t(-log1p(u - 1) / lambda);
	while(t < T_0) {
		retval.push_back(t);
		u = unif_rand();
		t -= log1p(u - 1) / lambda;
	}
	return wrap(retval);
	END_RCPP
}

/*
gen_homo_poisson <- function(lambda, T_0) {
	result <- c()
	t <- 0
	while(t <= T_0) {
		u <- runif(1)
		t <- t - log(u) / lambda
		if (t <= T_0) result <- append(result, t)
	}
	result
}
*/