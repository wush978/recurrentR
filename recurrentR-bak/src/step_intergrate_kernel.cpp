#include <Rcpp.h>

using namespace Rcpp;

//[[Rcpp::export]]
SEXP step_integrate_kernel(Function g, NumericVector x, NumericVector p, double a, double b) {
	BEGIN_RCPP
	double retval = 0;
	int i;
	while(x[i] < a & i < x.size()) {i++;}
	for(;i < x.size() & x[i] <= b;i++) {
		retval += as<double>(g(wrap(x[i])));
	}
	END_RCPP
}