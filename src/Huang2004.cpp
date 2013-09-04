#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP DSSj(SEXP RSk, SEXP RX) {
	BEGIN_RCPP
	NumericMatrix X(RX);
	NumericVector Sk(RSk);
	if (Sk.size() != X.nrow()) throw std::invalid_argument("length(Sk) == nrow(X) does not hold");
	if (X.nrow() == 0) throw std::invalid_argument("Empty X");
	NumericMatrix retval(X.nrow(), X.ncol());
	retval.fill(0);
	for(int j = 0;j < X.ncol();j++) {
		retval(0, j) = Sk[0] * X(0, j);
	}
	for(int i = 1;i < X.nrow();i++) {
		for(int j = 0;j < X.ncol();j++) {
			X(i,j) = X(i-1, j) + Sk[i] * X(i, j);
		}
	}
	return retval;
	END_RCPP
}