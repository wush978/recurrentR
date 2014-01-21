#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

typedef std::vector<double> NumVec;

struct StepFunction {
  NumVec x;
	NumVec y;
  bool lower_bound;
	
	StepFunction(NumericVector Rx, NumericVector Ry) 
	: x(as<NumVec>(Rx)), y(as<NumVec>(Ry)), lower_bound(false) { 
		if (x.size() + 1 != y.size()) throw std::invalid_argument("length(Rx) + 1 == length(Ry) does not hold");
	}
	~StepFunction() { }
	
	double call(double input) const {
    NumVec::const_iterator i;
		if (lower_bound) {
      i = std::lower_bound(x.begin(), x.end(), input);
		} else {
      i = std::upper_bound(x.begin(), x.end(), input);
		}
		return y[i - x.begin()];
	}
	
	SEXP sort_call(NumericVector input) {
		NumVec::const_iterator i = x.begin();
		NumericVector retval(input.size());
    if (lower_bound) {
    	for(int j = 0;j < input.size();j++) {
  			while(i != x.end() & *i < input[j]) { i++; }
  			retval[j] = y[i - x.begin()];
  		}
    } else {
  		for(int j = 0;j < input.size();j++) {
  			while(i != x.end() & *i <= input[j]) { i++; }
  			retval[j] = y[i - x.begin()];
  		}
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
  .field("lower_bound", &StepFunction::lower_bound)
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
