#include <Rcpp.h>

using namespace Rcpp;

inline int bupper(NumericVector x, double a) {
	return std::upper_bound(x.begin(), x.end(), a) - x.begin() + 1;
}

inline int blower(NumericVector x, double b) {
	return std::lower_bound(x.begin(), x.end(), b) - x.begin();
}

//[[Rcpp::export]]
SEXP substring_index(NumericVector x, double a, double b) {
	BEGIN_RCPP
	int begin_index = bupper(x, a);
	if (begin_index > x.size()) return R_NilValue;
	int end_index = blower(x, b);
	if (end_index == 0) return R_NilValue;
	if (begin_index > end_index) return R_NilValue;
	IntegerVector retval(end_index - begin_index + 1, 0);
	for(int i = begin_index;i < end_index;i++) {
		retval[i - begin_index] = i;
	}
	return retval;
	END_RCPP
}
