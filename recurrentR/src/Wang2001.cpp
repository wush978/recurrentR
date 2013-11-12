#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

typedef std::vector<double> NumVec;

struct StepFunction {
  NumVec x;
	NumVec y;
	
	StepFunction(NumericVector Rx, NumericVector Ry) 
	: x(as<NumVec>(Rx)), y(as<NumVec>(Ry)) { 
		if (x.size() + 1 != y.size()) throw std::invalid_argument("length(Rx) + 1 == length(Ry) does not hold");
	}
	~StepFunction() { }
	
	double call(double input) const {
		NumVec::const_iterator i = std::upper_bound(x.begin(), x.end(), input);
		return y[i - x.begin()];
	}
	
	SEXP sort_call(NumericVector input) {
		NumVec::const_iterator i = x.begin();
		NumericVector retval(input.size());
		for(int j = 0;j < input.size();j++) {
			while(*i <= input[j] & i != x.end()) { i++; }
			retval[j] = y[i - x.begin()];
		}
		return retval;
	}
	
private:
	StepFunction(const StepFunction&);
	void operator=(const StepFunction&);
};

//double step_integrate

RCPP_MODULE(StepFunction) {
	
	class_<StepFunction>("StepFunction")
	.constructor<NumericVector, NumericVector>()
	.field_readonly("x", &StepFunction::x)
	.field_readonly("y", &StepFunction::y)
	.method("call", &StepFunction::call)
	.method("sort_call", &StepFunction::sort_call)
	;
	
}

template<class T>
T* extract_ptr(SEXP s) {
	Rcpp::S4 s4(s);
	Rcpp::Environment env(s4);
	Rcpp::XPtr<T> xptr(env.get(".pointer"));
	return static_cast<T*>(R_ExternalPtrAddr(xptr));
}

RcppExport SEXP StepFunction_sort_call(SEXP Robj, SEXP Rx) {
	BEGIN_RCPP
	StepFunction& obj(*extract_ptr<StepFunction>(Robj));
	NumericVector x(Rx);
	return obj.sort_call(x);
	END_RCPP
}

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