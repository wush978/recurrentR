#include <Rcpp.h>

using namespace Rcpp;

//'@title DSSi
//'@description \eqn{\frac{d}{d \beta} SS_i} where \eqn{SS_i = \sum_{j=i}^n{S_j}}
//'and \eqn{S_i = Z_i e^{X_i \beta}}.
//'@return numeric matrix. \code{\retcal[,i]} is \eqn{DSS_i}.
RcppExport SEXP DSSi(SEXP RSk, SEXP RX) {
	BEGIN_RCPP
	NumericMatrix X(RX);
	NumericVector Sk(RSk);
	if (Sk.size() != X.nrow()) throw std::invalid_argument("length(Sk) == nrow(X) does not hold");
	if (X.nrow() == 0) throw std::invalid_argument("Empty X");
	NumericMatrix retval(X.ncol(), X.nrow());
	retval.fill(0);
	{
		int i = 0;
		int k = Sk.size() - i - 1;
		for(int j = 0;j < X.ncol();j++) {
			retval(j, k) = Sk[k] * X(k, j);
		}
	}
	
	for(int i = 1;i < Sk.size();i++) {
		int k = Sk.size() - i - 1;
		for(int j = 0;j < X.ncol();j++) {
			retval(j, k) = retval(j, k + 1) + Sk[k] * X(k, j);
		}
	}
	return retval;
	END_RCPP
}