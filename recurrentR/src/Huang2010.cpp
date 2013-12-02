#include <map>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;

typedef std::vector< std::vector< std::vector< std::vector<double> > > > Rao;
typedef std::vector< std::vector<int> > TIndex;
inline const double X_value_get(const double* X_value, const int* dim, int i, int j, int k) {
  return X_value[i + j * dim[0] + k * dim[0] * dim[1]];
}

SEXP init_vrao(const int p, const int n, std::vector<Rao*>& rao) {
  List retval(p);
  for(int r = 0;r < p;r++) {
    rao.push_back(new Rao);
    XPtr<Rao> temp(rao[r]);
    retval[r] = temp;
  }
  return retval;
}

// [[Rcpp::export]]
SEXP rao_gen(Function f, IntegerVector m, List t, NumericVector y, int p) {
  BEGIN_RCPP
  std::vector<Rao*> rao;
  int n = m.size();
  List retval(init_vrao(p, n, rao));

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
        c[0] = t_i[k];
        if (c[0] > y_min) {
          continue;
        }
        for(int l = 0;l < m[j];l++) {
          d[0] = t_j[l];
          if (d[0] > y_min) {
            continue;
          }
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
  return retval;
  END_RCPP
}



inline const double rao_array(const double* X_value, int i, int j, int col_k, int col_l, int* dim, int r) {
  return X_value_get(X_value, dim, i, col_l, r) +
    X_value_get(X_value, dim, j, col_k, r) -
    X_value_get(X_value, dim, i, col_k, r) -
    X_value_get(X_value, dim, j, col_l, r);
}

// [[Rcpp::export]]
SEXP t_index_gen(IntegerVector m, List t, NumericVector s) {
  BEGIN_RCPP
  XPtr< TIndex > retval(new TIndex());
  TIndex& t_index(*retval);
  int n = t.size();
  { // construct t_index
    t_index.resize(n);
    for(int i = 0;i < n;i++) {
      if (m[i] == 0) continue;
      t_index[i].resize(m[i]);
      NumericVector ti(wrap(t[i]));
      for(int j = 0;j < m[i];j++) {
        NumericVector::iterator it = std::lower_bound(s.begin(), s.end(), ti[j]);
        t_index[i][j] = it - s.begin();
      }
    }
  }
  return retval;
  END_RCPP
}

//[[Rcpp::export]]
SEXP t_index_inverse_gen(SEXP Rt_index, NumericVector s) {
  BEGIN_RCPP
  const TIndex& t_index(*XPtr< TIndex >(Rt_index));
  XPtr< TIndex > pretval_i(new TIndex()), pretval_j(new TIndex());
  TIndex &retval_i(*pretval_i), &retval_j(*pretval_j);
  retval_i.resize(s.size());
  retval_j.resize(s.size());
  for(int i = 0;i < t_index.size();i++) {
    for(int j = 0;j < t_index[i].size();j++) {
      
      retval_i[t_index[i][j]].push_back(i);
      retval_j[t_index[i][j]].push_back(j);
    }
  }
  return List::create(Named("i") = pretval_i, Named("j") = pretval_j);
  END_RCPP
}

//[[Rcpp::export]]
SEXP t_index_query(SEXP Rt_index, int i, int j, bool is_R_index = true) {
  BEGIN_RCPP
  const TIndex& t_index(*XPtr< TIndex >(Rt_index));
  if (is_R_index) return wrap(t_index[i-1][j-1] + 1);
  else return wrap(t_index[i][j]);
  END_RCPP
}

// [[Rcpp::export]]
SEXP rao_gen_array(NumericVector X_value, IntegerVector m, List t, NumericVector y, SEXP Rt_index) {
  BEGIN_RCPP
  int n = y.size();
  IntegerVector dim(wrap(X_value.attr("dim")));
  int p = dim[2];
  
  const TIndex& t_index(*XPtr< TIndex >(Rt_index));

  std::vector<Rao*> rao;
  List retval(init_vrao(p, n, rao));
  for(int r = 0;r < p;r++) {
    rao[r]->resize(n);
  }
  double y_min;
  double *pX_value = REAL(wrap(X_value));
  int *pdim = INTEGER(wrap(dim));
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
        if (t_i[k] > y_min) {
          continue;
        }
        for(int l = 0;l < m[j];l++) {
          if (t_j[l] > y_min) {
            continue;
          }
          int X_value_col_k = t_index[i][k];
          int X_value_col_l = t_index[j][l];
          for(int r = 0;r < p;r++) {
            (*rao[r])[i][j][k][l] = rao_array(pX_value, i, j, X_value_col_k, X_value_col_l, pdim, r);
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

const std::vector<Rao*> get_rao_vec(List& pRao_list) {
  std::vector<Rao*> rao;
  for(int r = 0;r < pRao_list.size();r++) {
    XPtr<Rao> pRao(wrap(pRao_list[r]));
    rao.push_back(&(*pRao));
  }
  return rao;
}

// [[Rcpp::export]]
SEXP S(List pRao_list, NumericVector beta) {
  BEGIN_RCPP
  const std::vector<Rao*> rao(get_rao_vec(pRao_list));
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
  const std::vector<Rao*> rao(get_rao_vec(pRao_list));
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

void set_g(std::vector<double>& g_ij, int i, int j, const std::vector<Rao*>& rao, double* beta) {
  std::fill(g_ij.begin(), g_ij.end(), 0);
  int p = rao.size();
  double temp;
  for(int ijg1 = 0;ijg1 < (*rao[0])[i][j].size();ijg1++) {
    for(int ijg2 = 0;ijg2 < (*rao[0])[i][j][ijg1].size();ijg2++) {
      temp = 0;
      for(int r = 0;r < p;r++) {
        temp += (*rao[r])[i][j][ijg1][ijg2] * beta[r];
      }
      temp = exp(temp);
      for(int r = 0;r < p;r++) {
        g_ij[r] = - (*rao[r])[i][j][ijg1][ijg2] * temp / (1 + temp);
      }
    }
  }
}

//[[Rcpp::export]]
SEXP V1_hat(List pRao_list, NumericVector beta) {
  BEGIN_RCPP
  const std::vector<Rao*> rao(get_rao_vec(pRao_list));
  int n = rao[0]->size();
  int p = rao.size();
  double* pbeta = &beta[0];
  NumericMatrix retval(p, p);
  retval.fill(0);
  std::vector<double> g_ij, g_ik;
  g_ij.resize(p);
  g_ik.resize(p);
  for(int i = 0;i < n - 2;i++) {
    for(int j = i + 1;j < n - 1;j++) {
      for(int k = j + 1;k < n;k++) {
        set_g(g_ij, i, j, rao, pbeta);
        set_g(g_ik, i, k, rao, pbeta);
        for(int r1 = 0;r1 < p;r1++) {
          for(int r2 = 0;r2 < p;r2++) {
            retval(r1, r2) += g_ij[r1] * g_ik[r2];
          }
        }
      }
    }
  }
  return retval;
  END_RCPP
}

//[[Rcpp::export]]
SEXP d_beta(NumericVector beta, NumericVector X_value, SEXP Rt_index) {
  BEGIN_RCPP
  const TIndex& t_index(*XPtr< TIndex >(Rt_index));
  IntegerVector dim(wrap(X_value.attr("dim")));
  int *pdim = INTEGER(wrap(dim));
  double *pX_value = REAL(wrap(X_value));
  NumericVector d(dim[1]);
  d.fill(0);
  for(int i = 0;i < t_index.size();i++) {
    for(int j = 0;j < t_index[i].size();j++) {
      int index = t_index[i][j];
      double temp = 0;
      for(int r = 0;r < dim[2];r++) {
        temp += X_value_get(pX_value, pdim, i, index, r) * beta[r];
      }
      d[index] += exp(- temp);
    }
  }
  return d;
  END_RCPP
}

//[[Rcpp::export]]
SEXP R_beta(NumericVector beta, NumericVector X_value, SEXP Rt_index, IntegerVector s_upper_index) {
  BEGIN_RCPP
  const TIndex& t_index(*XPtr< TIndex >(Rt_index));
  IntegerVector dim(wrap(X_value.attr("dim")));
  int *pdim = INTEGER(wrap(dim));
  double *pX_value = REAL(wrap(X_value));
  NumericVector R(dim[1]);
  R.fill(0);
  for(int i = 0;i < t_index.size();i++) {
    for(int j = 0;j < t_index[i].size();j++) {
      int index = t_index[i][j];
      double temp = 0;
      for(int r = 0;r < dim[2];r++) {
        temp += X_value_get(pX_value, pdim, i, index, r) * beta[r];
      }
      temp = exp(- temp);
      for(int k = index;k < s_upper_index[i];k++) {
        R[k] += temp;
      }
    }
  }
  return R;
  END_RCPP
}

//[[Rcpp::export]]
SEXP V_hat_tilde_Q_gen(NumericVector beta, NumericVector X_value, List t_index_inverse) {
  BEGIN_RCPP
  IntegerVector dim(wrap(X_value.attr("dim")));
  int *pdim = INTEGER(wrap(dim));
  double *pX_value = &X_value[0];
  XPtr< TIndex > pretval_i(wrap(t_index_inverse["i"])), pretval_j(wrap(t_index_inverse["j"]));
  const TIndex &retval_i(*pretval_i), &retval_j(*pretval_j);
  // X_value_get(pX_value, pdim)
  List retval(pdim[2]);
  std::vector< NumericVector* > retval_cache;
  for(int l = 0;l < pdim[2];l++) {
    retval_cache.push_back(new NumericVector(pdim[1] + 1));
    (*retval_cache[l])[0] = 0;
  }
  for(int s = 0;s < pdim[1];s++) {
    for(int t = 0;t < retval_i[s].size();t++) {
      int i = retval_i[s][t], j = retval_j[s][t];
      double temp = 0;
      for(int r = 0;r < pdim[2];r++) {
        temp += beta[r] * X_value_get(pX_value, pdim, i, j, r);
      }
      temp = exp(-temp);
      for(int l = 0;l < pdim[2];l++) {
        (*retval_cache[l])[s + 1] = (*retval_cache[l])[s] - temp * X_value_get(pX_value, pdim, i, j, l);
      }
    }
  }
  for(int l = 0;l < pdim[2];l++) {
    retval[l] = (*retval_cache[l]);
  }
  return retval;
  END_RCPP
}