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
	
private:
	StepFunction(const StepFunction&);
	void operator=(const StepFunction&);
};

//'@export StepFunction
RCPP_MODULE(recurrentR) {
	
	class_<StepFunction>("StepFunction")
	.constructor<NumericVector, NumericVector>()
	.method("call", &StepFunction::call)
	;
	
}