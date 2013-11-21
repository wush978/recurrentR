#include <map>
#include <Rcpp.h>
using namespace Rcpp;

typedef std::vector< std::vector< std::vector< std::vector<double> > > > Rao;

// [[Rcpp::export]]
SEXP rao_gen(Function f, IntegerVector m, List t, NumericVector y) {
  BEGIN_RCPP
  XPtr<Rao> retval(new Rao);
  Rao& rao(*retval);
  int n = m.size();
  rao.resize(n);
  NumericVector a(1), b(1), c(1), d(1);
  double y_min;
  for(int i = 0;i < n - 1;i++) {
    NumericVector t_i(wrap(t[i]));
    rao[i].resize(n);
    for(int j = i + 1;j < n;j++) {
      NumericVector t_j(wrap(t[j]));
      rao[i][j].resize(m[i]);
      y_min = (y[i] < y[j] ? y[i] : y[j]);
      for(int k = 0;k < m[i];k++) {
        rao[i][j][k].resize(m[j]);
        for(int l = 0;l < m[j];l++) {
          a[0] = i+1;
          b[0] = j+1;
          c[0] = t_i[k];
          d[0] = t_j[l];
          if (c[0] > y_min) {
            rao[i][j][k][l] = 0;
          }
          else {
            if (d[0] > y_min) {
              rao[i][j][k][l] = 0;
            }
            else {
              rao[i][j][k][l] = as<double>(f(wrap(a), wrap(b), wrap(c), wrap(d)));
            }
          }
        }
      }
    }
  }
  return retval;
  END_RCPP
}