#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
SEXP get_s_d(List t) {
  BEGIN_RCPP
  std::map<double, int> retval;
  for(int i = 0;i < t.size();i++) {
    NumericVector temp(wrap(t[i]));
    for(int j = 0;j < temp.size();j++) {
      retval[temp[j]] += 1;
    }
  }
  NumericVector s(retval.size());
  IntegerVector d(retval.size());
  int i = 0;
  for(std::map<double, int>::const_iterator element = retval.begin();element != retval.end(); element++) {
    s[i] = element->first;
    d[i++] = element->second;
  }
  return List::create(Named("s") = s, Named("d")= d);
  END_RCPP
}