#include <map>
#include <Rcpp.h>
using namespace Rcpp;

typedef std::vector< std::vector< std::vector< std::vector<double> > > > Rao;

// [[Rcpp::export]]
SEXP rao_gen(Function f, IntegerVector m, List t, NumericVector y, int p) {
  BEGIN_RCPP
  List retval(p);
  std::vector<Rao*> rao;
  for(int r = 0;r < p;r++) {
    rao.push_back(new Rao);
    XPtr<Rao> temp(rao[r]);
    retval[r] = temp;
  }
  int n = m.size();
  for(int r = 0;r < p;r++) {
    rao[r]->resize(n);
  }
  NumericVector a(1), b(1), c(1), d(1);
  double y_min;
  for(int i = 0;i < n - 1;i++) {
    NumericVector t_i(wrap(t[i]));
    for(int r = 0;r < p;r++) {
      (*rao[r])[i].resize(n);  
    }
    for(int j = i + 1;j < n;j++) {
      NumericVector t_j(wrap(t[j]));
      for(int r = 0;r < p;r++) {
        (*rao[r])[i][j].resize(m[i]);
      }
      y_min = (y[i] < y[j] ? y[i] : y[j]);
      for(int k = 0;k < m[i];k++) {
        for(int r = 0;r < p;r++) {
          (*rao[r])[i][j][k].resize(m[j], 0);
        }
        for(int l = 0;l < m[j];l++) {
          c[0] = t_i[k];
          d[0] = t_j[l];
          if (c[0] > y_min) {
            continue;
          }
          else {
            if (d[0] > y_min) {
              continue;
            }
            else {
              a[0] = i+1;
              b[0] = j+1;
              NumericVector temp(f(wrap(a), wrap(b), wrap(c), wrap(d)));
              for(int r = 0;r < p;r++) {
                (*rao[r])[i][j][k][l] = temp[r];
              }
            }
          }
        }
      }
    }
  }
  return retval;
  END_RCPP
}

// [[Rcpp::export]]
SEXP rao_export(SEXP pRao) {
  BEGIN_RCPP
  Rao& rao(*XPtr<Rao>(pRao));
  return wrap(rao);
  END_RCPP
}

// [[Rcpp::export]]
SEXP S(List pRao_list, NumericVector beta) {
  BEGIN_RCPP
  std::vector<Rao*> rao;
  for(int r = 0;r < pRao_list.size();r++) {
    XPtr<Rao> pRao(wrap(pRao_list[r]));
    rao.push_back(&(*pRao));
  }
  double temp;
  NumericVector retval(rao.size());
  retval.fill(0);
  int n = rao[0]->size();
  for(int i = 0;i < n - 1;i++) {
    for(int j = i + 1;j < n;j++) {
      for(int k = 0;k < (*rao[0])[i][j].size();k++) {
        for(int l = 0;l < (*rao[0])[i][j][k].size();l++) {
          temp = 0;
          for(int r = 0;r < rao.size();r++) {
            temp += (*rao[r])[i][j][k][l] * beta[r];
          }
          temp = exp(temp);
          for(int r = 0;r < rao.size();r++) {
            retval[r] -= (*rao[r])[i][j][k][l] * temp / (1 + temp);            
          }
        }
      }
    }
  }
  return retval;
  END_RCPP
}

// [[Rcpp::export]]
SEXP dS_over_dbeta(List pRao_list, NumericVector beta) {
  BEGIN_RCPP
  std::vector<Rao*> rao;
  for(int r = 0;r < pRao_list.size();r++) {
    XPtr<Rao> pRao(wrap(pRao_list[r]));
    rao.push_back(&(*pRao));
  }
  int n = rao[0]->size();
  int p = rao.size();
  NumericMatrix retval(p, p);
  retval.fill(0);
  double temp;
  for(int i = 0;i < n - 1;i++) {
    for(int j = i + 1;j < n;j++) {
      for(int k = 0;k < (*rao[0])[i][j].size();k++) {
        for(int l = 0;l < (*rao[0])[i][j][k].size();l++) {
          temp = 0;
          for(int r = 0;r < rao.size();r++) {
            temp += (*rao[r])[i][j][k][l] * beta[r];
          }
          temp = exp(temp);
          for(int r1 = 0;r1 < rao.size();r1++) {
            for(int r2 = 0;r2 < rao.size();r2++) {
              retval(r1,r2) -= (*rao[r1])[i][j][k][l] * (*rao[r2])[i][j][k][l] * temp / pow(1 + temp, 2);
            }
          }
        }
      }
    }
  }
  return retval;
  END_RCPP
}