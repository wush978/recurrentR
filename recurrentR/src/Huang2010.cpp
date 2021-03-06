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

void g_ij(const std::vector<Rao*>& rao, const int n, NumericVector& beta, const int i, const int j, double* retval) {
  double temp;
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

// [[Rcpp::export]]
SEXP g_ij(List pRao_list, NumericVector beta, int Ri, int Rj) {
  BEGIN_RCPP
  if (Ri >= Rj) throw std::invalid_argument("Ri < Rj");
  const std::vector<Rao*> rao(get_rao_vec(pRao_list));
  NumericVector retval(rao.size());
  retval.fill(0);
  g_ij(rao, rao[0]->size(), beta, Ri-1, Rj-1, &retval[0]);
//  double temp;
//  int n = rao[0]->size(), i = Ri - 1, j = Rj - 1;
//  for(int k = 0;k < (*rao[0])[i][j].size();k++) {
//    for(int l = 0;l < (*rao[0])[i][j][k].size();l++) {
//      temp = 0;
//      for(int r = 0;r < rao.size();r++) {
//        temp += (*rao[r])[i][j][k][l] * beta[r];
//      }
//      temp = exp(temp);
//      for(int r = 0;r < rao.size();r++) {
//        retval[r] -= (*rao[r])[i][j][k][l] * temp / (1 + temp);            
//      }
//    }
//  }
  return retval;
  END_RCPP
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
SEXP V_hat_tilde_Q_s_gen(NumericVector beta, NumericVector X_value, SEXP Rt_index) {
  BEGIN_RCPP
  const TIndex& t_index(*XPtr< TIndex >(Rt_index));
  IntegerVector dim(wrap(X_value.attr("dim")));
  int *pdim = INTEGER(wrap(dim));
  double *pX_value = REAL(wrap(X_value));
  NumericMatrix retval(pdim[1], pdim[2]);
  retval.fill(0);
  for(int i = 0;i < t_index.size();i++) {
    for(int j = 0;j < t_index[i].size();j++) {
      int index = t_index[i][j];
      double temp = 0;
      for(int r = 0;r < dim[2];r++) {
        temp += X_value_get(pX_value, pdim, i, index, r) * beta[r];
      }
      temp = exp(- temp);
      for(int r = 0;r < dim[2];r++) {
        retval(index, r) -= X_value_get(pX_value, pdim, i, index, r) * temp;
      }
    }
  }
  for(int s = 1;s < pdim[1];s++) {
    for(int r = 0;r < pdim[2];r++) {
      retval(s, r) += retval(s-1, r);
    }
  }
  return retval;
  END_RCPP
}

//[[Rcpp::export]]
SEXP V_hat_tilde_R_s_gen(NumericVector beta, NumericVector X_value, SEXP Rt_index, IntegerVector s_upper_index) {
  BEGIN_RCPP
  const TIndex& t_index(*XPtr< TIndex >(Rt_index));
  IntegerVector dim(wrap(X_value.attr("dim")));
  int *pdim = INTEGER(wrap(dim));
  double *pX_value = REAL(wrap(X_value));
  NumericMatrix retval(pdim[1], pdim[2]);
  retval.fill(0);
  for(int i = 0;i < t_index.size();i++) {
    for(int j = 0;j < t_index[i].size();j++) {
      int index = t_index[i][j];
      double temp = 0;
      for(int r = 0;r < dim[2];r++) {
        temp += X_value_get(pX_value, pdim, i, index, r) * beta[r];
      }
      temp = exp(- temp);
      for(int k = index;k < s_upper_index[i];k++) {
        for(int r = 0;r < dim[2];r++) {
          retval(k, r) += - X_value_get(pX_value, pdim, i, index, r) * temp;
        }
      }
    }
  }
  return retval;
  END_RCPP
}

// [[Rcpp::export]]
SEXP phi_i_j_hat_s_gen(NumericVector Q_s, NumericVector R_s, NumericMatrix V_Q, NumericMatrix V_R) {
  BEGIN_RCPP
  int s_len = Q_s.size();
  int X_dim = V_Q.ncol();
  NumericMatrix retval(s_len, X_dim);
  retval.fill(0);
  std::vector<double> term_1(X_dim, 0), term_2(X_dim, 0);
  double term_2_scalar;
  for(int si = s_len - 1;si >= 0;si--) {
    std::fill(term_1.begin(), term_1.end(), 0);
    std::fill(term_2.begin(), term_2.end(), 0);
    term_2_scalar = (Q_s[si] - (si == 0 ? 0 : Q_s[si - 1])) / pow(R_s[si], 2);
    for(int r = 0;r < X_dim;r++) {
      term_1[r] = V_Q(si, r) - (si == 0 ? 0 : V_Q(si - 1, r));
      term_1[r] = term_1[r] / R_s[si];
      term_2[r] = V_R(si, r) * term_2_scalar;
      retval(si, r) = (si == s_len - 1 ? 0 : retval(si+1, r)) + term_1[r] - term_2[r];
    }
  }
  return retval;
  END_RCPP
}

// [[Rcpp::export]]
SEXP exp_X_beta_d_Lambda_0(NumericVector beta, NumericVector X_value, NumericVector s_index_upper, NumericVector dLambda_s) {
  BEGIN_RCPP
  IntegerVector dim(wrap(X_value.attr("dim")));
  int *pdim = INTEGER(wrap(dim));
  double *pX_value = REAL(wrap(X_value));
  NumericVector retval(s_index_upper.size());
  retval.fill(0);
  double X_value_beta = 0;
  for(int i = 0;i < pdim[0];i++) {
    for(int si = 0; si < s_index_upper[i];si++) {
      X_value_beta = 0;
      for(int r = 0;r < beta.size();r++) {
        X_value_beta += X_value_get(pX_value, pdim, i, si, r) * beta[r];
      }
      X_value_beta = exp(-X_value_beta);
      retval[i] += X_value_beta * dLambda_s[si];
    }
  }
  return retval;
  END_RCPP
}

void X_V_2_inv_g_ij_exp_X_beta_d_Lambda_0(
  NumericVector& beta, double* pX_value, int* pdim, NumericVector& s_index_upper, 
  NumericVector Lambda_s, NumericVector dLambda_s, NumericMatrix& V2_inv,
  const std::vector<Rao*>& rao, const std::vector< std::vector< double* > >& kappa_i_j, int ii, int j,
  double* retval) {
  double X_value_beta = 0, X_value_beta_dkappa_ij_Lambda = 0, X_value_cache, X_V2_inv_g_ij_exp_Xbeta_dLambda;
  std::vector<double> V2_inv_cache(beta.size(), 0), g_ij_cache(beta.size(), 0);
  for(int k = 0;k < pdim[0];k++) {
    for(int si = 0; si < s_index_upper[k];si++) {
      X_value_beta = 0;
      memset(&V2_inv_cache[0], 0, V2_inv_cache.size() * sizeof(double));
      for(int r = 0;r < beta.size();r++) {
        X_value_cache = X_value_get(pX_value, pdim, k, si, r);
        X_value_beta += X_value_cache * beta[r];
        if (ii != j) {
          for(int r2 = 0;r2 < beta.size();r2++) {
            V2_inv_cache[r2] += X_value_cache * V2_inv(r, r2);
          }
        }
      }
      memset(&g_ij_cache[0], 0, sizeof(double) * g_ij_cache.size());
      X_V2_inv_g_ij_exp_Xbeta_dLambda = 0;
      if (ii != j) {
        g_ij(rao, rao[0]->size(), beta, std::min(ii, j), std::max(ii, j), &g_ij_cache[0]);
        for(int r = 0;r < beta.size();r++) {
          X_V2_inv_g_ij_exp_Xbeta_dLambda += V2_inv_cache[r] * g_ij_cache[r];
        }
      }
      X_value_beta = exp(X_value_beta);
      X_value_beta_dkappa_ij_Lambda = X_value_beta * (kappa_i_j[ii][j][si] * Lambda_s[si] - (si == 0 ? 0 : kappa_i_j[ii][j][si - 1] * Lambda_s[si - 1]));
      if (ii != j) {
        X_V2_inv_g_ij_exp_Xbeta_dLambda = X_V2_inv_g_ij_exp_Xbeta_dLambda * X_value_beta * dLambda_s[si];
      }
//      if (ii == 0 & j == 3) Rprintf("k: %d si:%d %.8f\n", k, si, X_V2_inv_g_ij_exp_Xbeta_dLambda);
      retval[k] += X_V2_inv_g_ij_exp_Xbeta_dLambda + X_value_beta_dkappa_ij_Lambda;
    }
  }
}

//[[Rcpp::export]]
SEXP X_V_2_inv_g_ij_exp_X_beta_d_Lambda_0(
  NumericVector beta, NumericVector X_value, NumericVector s_index_upper, 
  NumericVector Lambda_s, NumericVector dLambda_s, NumericMatrix& V2_inv,
  List pRao_list, List Rkappa_i_j_s) {
  BEGIN_RCPP
  IntegerVector dim(wrap(X_value.attr("dim")));
  int *pdim = INTEGER(wrap(dim));
  double *pX_value = REAL(wrap(X_value));
  const std::vector<Rao*> rao(get_rao_vec(pRao_list));
  std::vector< std::vector< double* > > kappa_i_j_s;
  kappa_i_j_s.resize(pdim[0], std::vector< double* >(pdim[0], (double*) NULL));
  List retval(pdim[0]);
  for(int i = 0;i < pdim[0];i++) {
    List temp(wrap(Rkappa_i_j_s[i]));
    List retval_element(pdim[0]);
    for(int j = 0;j < pdim[0];j++) {
      kappa_i_j_s[i][j] = REAL(wrap(temp[j]));
      NumericVector retval_element_proxy(pdim[0]);
      retval_element_proxy.fill(0);
      X_V_2_inv_g_ij_exp_X_beta_d_Lambda_0(beta, pX_value, pdim, s_index_upper, 
        Lambda_s, dLambda_s, V2_inv, rao, kappa_i_j_s, i, j, &retval_element_proxy[0]
      );
      retval_element[j] = retval_element_proxy;
    }
    retval[i] = retval_element;
  }
  return retval;
  END_RCPP
}

