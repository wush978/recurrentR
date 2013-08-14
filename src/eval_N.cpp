#include <Rcpp.h>

using namespace Rcpp;


//'@title Evaluate N
//'
//'@param s Numeric vector, sorted
//'@param y Numeric vector, sorted
//'@param m Integer vector, sorted by y
//[[Rcpp::export]]
SEXP eval_N(NumericVector s, NumericVector y, IntegerVector m) {
	BEGIN_RCPP
	IntegerVector retval(s.size(), 0);
	int offset = 0, j = 0;
	for(int i = 0;i < s.size();i++) {
		while(j < y.size() & y[j] < s[i]) {
			offset -= m[j];
			j++;
		}
		retval[i] = offset;
	}
	return retval;
	END_RCPP
}