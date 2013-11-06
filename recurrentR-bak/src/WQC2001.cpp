#include <Rcpp.h>

using namespace Rcpp;

template<typename T>
bool leq(T a, T b) {
	return a <= b;
}

inline int bupper(NumericVector x, double a) {
	return std::upper_bound(x.begin(), x.end(), a, leq<double>) - x.begin() + 1;
}

inline int blower(NumericVector x, double b) {
	return std::lower_bound(x.begin(), x.end(), b) - x.begin();
}


RcppExport SEXP substring_index(SEXP Rx, SEXP Ra, SEXP Rb) {
	BEGIN_RCPP
	NumericVector x(Rx);
	double a(as<double>(Ra));
	double b(as<double>(Rb));
	int begin_index = bupper(x, a);
	if (begin_index > x.size()) return R_NilValue;
	int end_index = blower(x, b);
	if (end_index == 0) return R_NilValue;
	if (begin_index > end_index) return R_NilValue;
	IntegerVector retval(end_index - begin_index + 1, 0);
	for(int i = begin_index;i <= end_index;i++) {
		retval[i - begin_index] = i;
	}
	return retval;
	END_RCPP
}
