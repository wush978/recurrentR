pRao.gen <- function(obj) {
  key <- "pRao"
  if (!is_cache(obj, key)) {
    rao.temp <- function(i,j,u,v) {
      obj@X[[i]](v) + obj@X[[j]](u) - obj@X[[i]](u) - obj@X[[j]](v)
    }
    m <- m.gen(obj)
    obj@cache[[key]] <- rao_gen(rao.temp, m, obj@t, obj@y, obj@X_dim)
  }
  obj@cache[[key]]
}

X.value.gen <- function(obj, X.value = NULL) {
  key <- "X.value"
  if (is.null(X.value)) {
    if (!is_cache(obj, key)) {
      X.value <- array(0.0, c(obj@n, length(obj@s), obj@X_dim))
      for(i in seq_along(obj@X)) {
        for(j in seq_along(obj@s)) {
          X.value[i,j,] <- obj@X[[i]](obj@s[j])
        }
      }
      obj@cache[[key]] <- X.value
    }
  } else {
    stopifnot(dim(X.value) == c(obj@n, length(obj@s), obj@X_dim))
    obj@cache[[key]] <- X.value
  }
  return(obj@cache[[key]])
}

t_index.gen <- function(obj) {
  key <- "t_index"
  if (!is_cache(obj, key)) {
    obj@cache[[key]] <- t_index_gen(m.gen(obj), obj@t, obj@s)
  }
  obj@cache[[key]]
}

pRao.gen.array <- function(obj, X.value = NULL) {
  stopifnot(length(dim(X.value)) == 3)
  key <- "pRao"
  if (!is_cache(obj, key)) {
    obj@cache[[key]] <- rao_gen_array(X.value.gen(obj, X.value), m.gen(obj), obj@t, obj@y, t_index.gen(obj))
  }
  obj@cache[[key]]
}

S.gen <- function(obj) {
  key <- "S"
  if (!is_cache(obj, key)) {
    denom <- choose(obj@n, 2)
    obj@cache[[key]] <- function(beta) {
      S(pRao.gen(obj), beta) / denom
    }
  }
  obj@cache[[key]]
}

dS_over_dbeta.gen <- function(obj) {
  key <- "dS_over_dbeta"
  if (!is_cache(obj, key)) {
    denom <- choose(obj@n, 2)
    obj@cache[[key]] <- function(beta) {
      dS_over_dbeta(pRao.gen(obj), beta) / denom
    }
  }
  obj@cache[[key]]
}

beta.hat.gen <- function(obj) {
  key <- "beta.hat"
  if (!is_cache(obj, key)) {
    slv <- nleqslv(rep(0, obj@X_dim), S.gen(obj), dS_over_dbeta.gen(obj))
    obj@cache[[key]] <- slv$x
  }
  obj@cache[[key]]
}

V1.hat.gen <- function(obj) {
  key <- "V1.hat"
  if (!is_cache(obj, key)) {
    denom <- obj@n * choose(obj@n - 1, 2)
    obj@cache[[key]] <- 4 * V1_hat(pRao.gen(obj), beta.hat.gen(obj)) / denom
  }
  obj@cache[[key]]
}

V2.hat.gen <- function(obj) {
  key <- "V2.hat"
  if (!is_cache(obj, key)) {
    dS_over_dbeta <- dS_over_dbeta.gen(obj)
    beta.hat <- beta.hat.gen(obj)
    obj@cache[[key]] <- -dS_over_dbeta(beta.hat)
  }
  obj@cache[[key]]
}

d_beta.R <- function(obj) {
  beta.hat <- beta.hat.gen(obj)
  retval <- numeric(length(obj@s))
  for(i in seq_along(obj@t)) {
    for(j in seq_along(obj@t[[i]])) {
      index <- which(obj@t[[i]][j] == obj@s)
      retval[index] <- retval[index] + exp(- obj@X[[i]](obj@t[[i]][j]) %*% beta.hat)
    }
  }
  retval / obj@n
}

d_beta.gen <- function(obj) {
  key <- "d_beta"
  if (!is_cache(obj, key)) {
    beta.hat <- beta.hat.gen(obj)
    X.value <- X.value.gen(obj)
    t_index <- t_index.gen(obj)
    obj@cache[[key]] <- d_beta(beta.hat, X.value, t_index) / obj@n
  }
  obj@cache[[key]]
}

R_beta.R <- function(obj) {
  beta.hat <- beta.hat.gen(obj)
  retval <- numeric(length(obj@s))
  for(i in seq_along(obj@t)) {
    index_upper <- tail(which(obj@s <= obj@y[i]), 1)
    if (length(index_upper) == 0) next
    for(j in seq_along(obj@t[[i]])) {
      index_lower <- which(obj@t[[i]][j] == obj@s)
      value <- exp(- obj@X[[i]](obj@t[[i]][j]) %*% beta.hat)
      if (index_lower > index_upper) next
      for(k in index_lower:index_upper) {
        retval[k] <- retval[k] + value
      }
    }
  }
  retval / obj@n
}

s_index_upper.gen <- function(obj) {
  key <- "s_index_upper"
  if (!is_cache(obj, key)) {
    y_vs_s <- outer(obj@y, obj@s, ">=")
    obj@cache[[key]] <- apply(y_vs_s, 1, function(a) {
      retval <- which(!a)
      ifelse(length(retval) == 0, length(obj@s), min(retval) - 1)
    })
  }
  obj@cache[[key]]
}

R_beta.gen <- function(obj) {
  key <- "R_beta"
  if (!is_cache(obj, key)) {
    obj@cache[[key]] <- R_beta(beta.hat.gen(obj), X.value.gen(obj), t_index.gen(obj), s_index_upper.gen(obj)) / obj@n
  }
  obj@cache[[key]]
}

Lambda_0.hat_Huang2010.gen <- function(obj) {
  key <- "Lambda_0.hat_Huang2010"
  if (!is_cache(obj, key)) {
    d_beta <- d_beta.gen(obj)
    R_beta <- R_beta.gen(obj)
    x <- obj@s
    y <- append(rev(cumprod(rev(1 - d_beta/R_beta))), 1)
    f <- stepfun(x, y)
    obj@cache[[key]] <- f
  }
  obj@cache[[key]]
}

Lambda_0.hat_Huang2010.Y.adj.gen <- function(obj) {
  key <- "Lambda_0.hat_Huang2010.Y.adj"
  if (!is_cache(obj, key)) {
    beta.hat <- beta.hat.gen(obj)
    X.value <- X.value.gen(obj)
    Lambda_0.hat_Huang2010.s <- Lambda_0.hat_Huang2010.gen(obj)(obj@s)
    s.bar <- c(0, obj@s)
    retval <- rep(0, obj@n)
    for(i in 1:obj@n) {
      index <- which(obj@s <= obj@y[i])
      X_u_beta <- as.vector(matrix(X.value[i,index,], ncol = obj@X_dim) %*% beta.hat)
      retval[i] <- diff(c(0, Lambda_0.hat_Huang2010.s[index])) %*% exp(X_u_beta)
    }
    obj@cache[[key]] <- retval
  }
  obj@cache[[key]]
}

Lambda_0.hat_Huang2010.Y.adj.inv.gen <- function(obj) {
  key <- "Lambda_0.hat_Huang2010.Y.adj.inv"
  if (!is_cache(obj, key)) {
    Lambda_0.hat_Huang2010.Y.adj <- Lambda_0.hat_Huang2010.Y.adj.gen(obj)
    retval <- 1 / Lambda_0.hat_Huang2010.Y.adj
    retval[Lambda_0.hat_Huang2010.Y.adj == 0] <- 0
    obj@cache[[key]] <- retval
  }
  obj@cache[[key]]
}


gamma.bar.hat_Huang2010.gen <- function(obj) {
  key <- "gamma.bar.hat_Huang2010"
  if (!is_cache(obj, key)) {
    beta.hat <- beta.hat.gen(obj)
    X.value <- X.value.gen(obj)
    Lambda_0.hat_Huang2010.Y.adj.inv <- Lambda_0.hat_Huang2010.Y.adj.inv.gen(obj)
    W.bar <- cbind(1, obj@W)
    b <- m.gen(obj) * Lambda_0.hat_Huang2010.Y.adj.inv
    gamma.bar.hat <- rep(0, ncol(W.bar))
    g <- function(gamma) t(W.bar) %*% (b - exp(W.bar %*% gamma))
    g.grad <- function(gamma) - t(W.bar) %*% diag(c(exp(W.bar %*% gamma)), nrow(W.bar), nrow(W.bar)) %*% W.bar
    slv <- nleqslv(gamma.bar.hat, g, jac=g.grad)
    if (sum(abs(g(slv$x))) > obj@tol) stop("Failed to converge during solving gamma")
    obj@cache[[key]] <- slv$x
  }
  obj@cache[[key]]
}

Q.hat_Huang2010 <- function(obj) {
  beta.hat <- beta.hat.gen(obj)
  function(u) {
    retval <- 0
    for(i in seq_along(obj@t)) {
      for(j in seq_along(obj@t[[i]])) {
        if (obj@t[[i]][j] > u) next
        retval <- retval + exp(- obj@X[[i]](obj@t[[i]][j]) %*% beta.hat)
      }
    }
    as.vector(retval/ obj@n)
  }
}

Q.hat_Huang2010.c <- function(obj) {
  key <- "Q.hat_Huang2010.c"
  if (!is_cache(obj, key)) {
    N <- cumsum(d_beta.gen(obj))
    obj@cache[[key]] <- new(StepFunction, obj@s, c(0, N))
  }
  obj@cache[[key]]
}

R.hat_Huang2010.R <- function(obj) {
  beta.hat <- beta.hat.gen(obj)
  function(u) {
    retval <- 0
    for(i in seq_along(obj@t)) {
      if (u > obj@y[i]) next
      for(j in seq_along(obj@t[[i]])) {
        if (u < obj@t[[i]][j]) next
        retval <- retval + exp(- obj@X[[i]](obj@t[[i]][j]) %*% beta.hat)
      }
    }
    as.vector(retval) / obj@n
  }
}

# This implementation is incorrect. However, R.hat_Huang2010.c$sort_call(obj@s) should correct.
R.hat_Huang2010.c <- function(obj) {
  key <- "R.hat_Huang2010.c"
  if (!is_cache(obj, key)) {
    R_beta <- R_beta.gen(obj)
    obj@cache[[key]] <- new(StepFunction, obj@s, c(0, R_beta))  
  }
}

V.hat_tilde_Q.R <- function(obj) {
  X.value <- X.value.gen(obj)
  beta.hat <- beta.hat.gen(obj)
  function(u) {
    retval <- numeric(obj@X_dim)
    for(i in seq_along(obj@t)) {
      for(j in seq_along(obj@t[[i]])) {
        if (obj@t[[i]][j] > u) next
        temp <- obj@X[[i]](obj@t[[i]][j])
        retval <- retval - exp( - as.vector(temp %*% beta.hat) ) * temp
      }
    }
    retval / obj@n
  }
}

# t_index.inverse.gen <- function(obj) {
#   key <- "t_index.inverse"
#   if (!is_cache(obj, key)) {
#     obj@cache[[key]] <- t_index_inverse_gen(t_index.gen(obj), obj@s)
#   }
#   obj@cache[[key]]
# }

# V.hat_tilde_Q.gen <- function(obj) {
#   key <- "V.hat_tilde_Q"
#   if (!is_cache(obj, key)) {
#     temp <- V_hat_tilde_Q_gen(beta.hat.gen(obj), X.value.gen(obj), t_index.inverse.gen(obj))
#     obj@cache[[key]] <- lapply(temp, function(d) {
#       new(StepFunction, obj@s, d/obj@n)
#     })
#   }
#   obj@cache[[key]]
# }

V.hat_tilde_Q.s.gen <- function(obj) {
  key <- "V.hat_tilde_Q.s"
  if (!is_cache(obj, key)) {
    obj@cache[[key]] <- V_hat_tilde_Q_s_gen(beta.hat.gen(obj), X.value.gen(obj), t_index.gen(obj)) / obj@n
  }
  obj@cache[[key]]
}

V.hat_tilde_R.R <- function(obj) {
  beta.hat <- beta.hat.gen(obj)
  X.value <- X.value.gen(obj)
  function(u) {
    retval <- numeric(obj@X_dim)
    for(i in seq_along(obj@t)) {
      if (u > obj@y[i]) next
      for(j in seq_along(obj@t[[i]])) {
        if (u < obj@t[[i]][j]) next
        temp <- obj@X[[i]](obj@t[[i]][j])
        retval <- retval - as.vector(exp(- temp %*% beta.hat)) * temp
      }
    }
    retval / obj@n
  }
}

V.hat_tilde_R.s.gen <- function(obj) {
  key <- "V.hat_tilde_R"
  if (!is_cache(obj, key)) {
    obj@cache[[key]] <- V_hat_tilde_R_s_gen(beta.hat.gen(obj), X.value.gen(obj), t_index.gen(obj), s_index_upper.gen(obj)) / obj@n
  }
  obj@cache[[key]]
}

phi.i.j.hat.s.R <- function(obj) {
  key <- "phi.i.j.hat"
  if (!is_cache(obj, key)) {
    retval <- list()
    beta.hat <- beta.hat.gen(obj)
    V.hat_tilde_Q.s <- V.hat_tilde_Q.s.gen(obj)
    V.hat_tilde_R.s <- V.hat_tilde_R.s.gen(obj)
    R_beta <- R_beta.gen(obj)
    Q.hat <- Q.hat_Huang2010.c(obj)
    d_beta <- d_beta.gen(obj)
    V2.hat <- V2.hat.gen(obj)
    pRao.array <- pRao.gen.array(obj, X.value.gen(obj))
    generator <- function(i,j) {
      force(i)
      force(j)
      function(t) {
        if (i == j) return(0)
        term_1 <- rep(0, obj@X_dim)
        term_2 <- rep(0, obj@X_dim)
        for(si in seq_along(obj@s)) {
          if (obj@s[si] < t) next
          if (si == 1) {
            term_1 <- term_1 + as.vector(V.hat_tilde_Q.s[si,]) / R_beta[si]
            term_2 <- term_2 + Q.hat$sort_call(c(0, obj@s[1])) %*% V.hat_tilde_R.s[si,] / R_beta[si]^2
          } else {
            term_1 <- term_1 + as.vector(V.hat_tilde_Q.s[si,] - V.hat_tilde_Q.s[si - 1,]) / R_beta[si]
            term_2 <- term_2 + diff(Q.hat$sort_call(c(obj@s[si-1], obj@s[si]))) * V.hat_tilde_R.s[si,] / R_beta[si]^2
            
          }
#           term_1 <- term_1 + as.vector(V.hat_tilde_Q.s[si,] - ifelse(si == 1, rep(0, obj@X_dim), V.hat_tilde_Q.s[si - 1,])) / R_beta[si]
#           term_2 <- term_2 + ifelse(si == 1, , ) * V.hat_tilde_R.s[si,] / R_beta[si]^2
        }
        as.vector(solve(t(V2.hat), (term_1 - term_2)) %*% g_ij(pRao.array, beta.hat, i, j))
      }
    }
    for(i in seq_len(obj@n)) {
      retval[[i]] <- list()
      for(j in seq_len(obj@n)) {
        retval[[i]][[j]] <- generator(i,j)
      }
    }
    obj@cache[[key]] <- retval
  }
  obj@cache[[key]]
}

phi.i.j.hat.s.gen <- function(obj) {
  key <- "phi.i.j.hat.s"
  if (!is_cache(obj, key)) {
    beta.hat <- beta.hat.gen(obj)
    V.hat_tilde_Q.s <- V.hat_tilde_Q.s.gen(obj)
    V.hat_tilde_R.s <- V.hat_tilde_R.s.gen(obj)
    R_beta <- R_beta.gen(obj)
    Q.hat.s <- Q.hat_Huang2010.c(obj)$sort_call(obj@s)
    V2.hat <- V2.hat.gen(obj)
    pRao.array <- pRao.gen.array(obj, X.value.gen(obj))
    term_1.mat <- phi_i_j_hat_s_gen(Q.hat.s, R_beta, V.hat_tilde_Q.s, V.hat_tilde_R.s)
    term_1.slv <- solve(V2.hat, t(term_1.mat))
    generator <- function(i,j) {
      force(i)
      force(j)
      if (i == j) return(new(StepFunction, obj@s, rep(0, length(obj@s) + 1)))
      retval <- new(StepFunction, obj@s, c(as.vector(g_ij(pRao.array, beta.hat, min(i,j), max(i,j)) %*% term_1.slv), 0))
      retval$lower_bound <- TRUE
      retval
    }
    retval <- list()
    for(i in seq_len(obj@n)) {
      retval[[i]] <- list()
      for(j in seq_len(obj@n)) {
        retval[[i]][[j]] <- generator(i,j)
      }
    }
    obj@cache[[key]] <- retval
  }
  obj@cache[[key]]
}

psi_i.hat.R.gen <- function(obj) {
  key <- "psi_i.hat.R"
  if (!is_cache(obj, key)) {
    retval <- list()
    beta.hat <- beta.hat.gen(obj)
    R_beta <- R_beta.gen(obj)
    Q.hat <- Q.hat_Huang2010.c(obj)
    t_index <- t_index.gen(obj)
    generator <- function(i) {
      force(i)
      function(t) {
        term_1 <- 0
        term_2 <- 0
        for(j in seq_along(obj@t[[i]])) {
          si <- t_index_query(t_index, i, j)
          term_1_delta <- 1 / R_beta[si]
          while(ifelse(si <= length(obj@s), obj@s[si] <= obj@y[i], FALSE)) {
            if (si == 1) {
              term_2_delta <- exp(- obj@X[[i]](obj@t[[i]][j]) %*% beta.hat) * 
                diff(Q.hat$sort_call(c(0, obj@s[1]))) / 
                (R_beta[si]^2)
            } else {
              term_2_delta <- exp(- obj@X[[i]](obj@t[[i]][j]) %*% beta.hat) * 
                diff(Q.hat$sort_call(c(obj@s[si-1], obj@s[si]))) / 
                (R_beta[si]^2)
            }
            if (obj@s[si] > t) {
#               print(sprintf("si: %d term_2_delta: %.8f", si, term_2_delta))
              term_2 <- term_2 + term_2_delta
            }
            si <- si + 1
          }
          if (t >= obj@t[[i]][j]) next
          term_1 <- term_1 + term_1_delta
        }
#         print(sprintf("term_1: %0.8f term_2: %0.8f", term_1, term_2))
        as.vector(term_1 - term_2)
      }
    }
    for(i in seq_along(obj@t)) {
      retval[[i]] <- generator(i)
    }
    obj@cache[[key]] <- retval
  }
  obj@cache[[key]]
}

psi_i.hat.gen <- function(obj) {
  key <- "psi_i.hat"
  if (!is_cache(obj, key)) {
    retval <- list()
    beta.hat <- beta.hat.gen(obj)
    R_beta <- R_beta.gen(obj)
    Q.hat <- Q.hat_Huang2010.c(obj)
    t_index <- t_index.gen(obj)
    s_index_upper <- s_index_upper.gen(obj)
    X.value <- X.value.gen(obj)
    t_ij_vs_s <- recurrentR:::t_ij_vs_s.gen(obj)
    generator <- function(i) {
      force(i)
      if (length(obj@t[[i]]) == 0) {
        new(StepFunction, numeric(0), 0)
      } else {
        sj <- unlist(lapply(seq_along(obj@t[[i]]), function(j) t_index_query(t_index, i, j)))
        R_beta.inv <- 1/R_beta[sj] # on obj@s[si]
        retval1 <- new(StepFunction, obj@s[sj], c(rev(cumsum(rev(R_beta.inv))), 0))
        si_u <- s_index_upper[i]
        value <- diff(c(0, Q.hat$sort_call(obj@s[seq_len(si_u)]))) / 
          R_beta[seq_len(si_u)]^2
        indicator.mat <- local({
          indicator.mat <- matrix(t_ij_vs_s[[i]][seq_len(si_u),], nrow = si_u)
          indicator.mat <- t(indicator.mat)
          t(indicator.mat * as.vector(exp(- X.value[i,sj,] %*% beta.hat)))
        })
        value.expand <- indicator.mat * value
        retval2 <- new(StepFunction, obj@s[seq_len(si_u)], c(apply(apply(value.expand, 2, function(a) rev(cumsum(rev(a)))), 1, sum), 0))
#         browser()
        retval1 - retval2
      }
    }
    for(i in seq_along(obj@t)) {
      retval[[i]] <- generator(i)
    }
    obj@cache[[key]] <- retval
  }
  obj@cache[[key]]
}


kappa_i.j.gen <- function(obj) {
  key <- "kappa_i.j"
  if (!is_cache(obj, key)) {
    obj@cache[[key]] <- list()
    psi_i.hat <- psi_i.hat.gen(obj)
    phi.i.j.hat.s <- phi.i.j.hat.s.gen(obj)
    generator <- function(i,j) {
      force(i)
      force(j)
      function(t) {
        phi.i.j.hat.s[[i]][[j]]$call(t) - (psi_i.hat[[i]]$call(t) + psi_i.hat[[j]]$call(t)) / 2
      }
    }
    for(i in seq_len(obj@n)) {
      obj@cache[[key]][[i]] <- list()
      for(j in seq_len(obj@n)) {
        obj@cache[[key]][[i]][[j]] <- generator(i,j)
      }
    }
  }
  obj@cache[[key]]
}

kappa_i.j.s.gen <- function(obj) {
  key <- "kappa_i.j.s"
  if (!is_cache(obj, key)) {
    obj@cache[[key]] <- list()
    psi_i.hat <- psi_i.hat.gen(obj)
    psi_i.hat.s <- lapply(psi_i.hat, function(o) o$sort_call(obj@s))
    phi.i.j.hat.s <- phi.i.j.hat.s.gen(obj)
    generator <- function(i,j) {
      force(i)
      force(j)
      phi.i.j.hat.s[[i]][[j]]$sort_call(obj@s) - (psi_i.hat.s[[i]] + psi_i.hat.s[[j]]) / 2
    }
    for(i in seq_len(obj@n)) {
      obj@cache[[key]][[i]] <- list()
      for(j in seq_len(obj@n)) {
        obj@cache[[key]][[i]][[j]] <- generator(i,j)
      }
    }    
  }
  obj@cache[[key]]
}


Lambda_0.hat.assymptotic.var.gen <- function(obj) {
  key <- "Lambda_0.hat.assymptotic.var"  
  stopifnot(obj@n > 2)
  if (!is_cache(obj, key)) {
    Lambda_0.hat_Huang2010 <- Lambda_0.hat_Huang2010.gen(obj)
    kappa_i.j <- kappa_i.j.gen(obj)
    obj@cache[[key]] <- function(t) {
      kappa.sum <- 0
      for(i in seq_len(obj@n)) {
        for(j in 1:(obj@n - 1)) {
          if (j == i) next
          for(k in (j+1):obj@n) {
            kappa.sum <- kappa.sum + kappa_i.j[[min(i,j)]][[max(i,j)]](t) * kappa_i.j[[min(i,k)]][[max(i,k)]](t)
          }
        }
      }
      8 * (Lambda_0.hat_Huang2010(t)^2) * kappa.sum / (obj@n * (obj@n - 1) * (obj@n - 2))
    }
  }
  obj@cache[[key]]
}

xi.i.j.hat.R <- function(obj) {
  key <- "xi.i.j.hat.R"
  if (!is_cache(obj, key)) {
    W.bar <- cbind(1, obj@W)
    m <- m.gen(obj)
    X.value <- X.value.gen(obj)  
    beta.hat <- beta.hat.gen(obj)
    Lambda_0.hat.s <- Lambda_0.hat_Huang2010.gen(obj)(obj@s)
    dLambda_0.hat.s <- diff(c(0, Lambda_0.hat.s))
    V2.hat <- V2.hat.gen(obj)
    V2.hat.inv <- solve(V2.hat)
    pRao_list <- pRao.gen(obj)
    gamma.bar.hat <- gamma.bar.hat_Huang2010.gen(obj)
    kappa_i.j.s <- kappa_i.j.s.gen(obj)
    s_index_upper <- s_index_upper.gen(obj)
    int_y_Lambda_0 <- numeric(obj@n)
    for(i in seq_len(obj@n)) {
      for(si in seq_along(obj@s)) {
        s <- obj@s[si]
        if (s > obj@y[i]) next
        int_y_Lambda_0[i] <- int_y_Lambda_0[i] + exp(- X.value[i,si,] %*% beta.hat) * (dLambda_0.hat.s[si])
      }
    }
    get_term <- function(i) {
      force(i)
      as.vector((m[i] / int_y_Lambda_0[i] - exp(W.bar[i,] %*% gamma.bar.hat)) * W.bar[i,] / 2)
    }
    s.0 <- numeric(length(obj@s))
    retval <- list()
    for(i in seq_len(obj@n)) {
      retval[[i]] <- list()
      for(j in seq_len(obj@n)) {
        dkappa.Lambda_0 <- diff(c(0, kappa_i.j.s[[i]][[j]] * Lambda_0.hat.s))
        tempk <- numeric(obj@n)
        for(k in seq_len(obj@n)) {
          temp2 <- exp(X.value[k,,] %*% beta.hat) * dkappa.Lambda_0
          temp1 <- if (i == j) s.0 else (X.value[k,,] %*% V2.hat.inv %*% g_ij(pRao_list, beta.hat, min(i, j), max(i, j))) * (exp(X.value[k,,] %*% beta.hat)) * dLambda_0.hat.s
          tempk[k] <- sum(temp1[1:s_index_upper[k]] + temp2[1:s_index_upper[k]])
        }
        retval[[i]][[j]] <- apply(((- m * W.bar) / int_y_Lambda_0^2) * tempk, 2, mean) + get_term(i) + get_term(j)
      }
    }
    obj@cache[[key]] <- retval
  }
  obj@cache[[key]]
}

xi.i.j.hat.gen <- function(obj) {
  key <- "xi.i.j.hat"
  if (!is_cache(obj, key)) {
    W.bar <- cbind(1, obj@W)
    m <- m.gen(obj)
    X.value <- X.value.gen(obj)  
    beta.hat <- beta.hat.gen(obj)
    Lambda_0.hat.s <- Lambda_0.hat_Huang2010.gen(obj)(obj@s)
    dLambda_0.hat.s <- diff(c(0, Lambda_0.hat.s))
    V2.hat <- V2.hat.gen(obj)
    V2.hat.inv <- solve(V2.hat)
    pRao_list <- pRao.gen(obj)
    gamma.bar.hat <- gamma.bar.hat_Huang2010.gen(obj) 
    kappa_i.j.s <- kappa_i.j.s.gen(obj)
    s_index_upper <- s_index_upper.gen(obj)
    exp_X_beta_d_Lambda_0.n <- exp_X_beta_d_Lambda_0(beta.hat, X.value, s_index_upper, dLambda_0.hat.s)
    X_V_2_inv_g_ij_exp_X_beta_d_Lambda_0.n.n <- X_V_2_inv_g_ij_exp_X_beta_d_Lambda_0(
      beta.hat, X.value, s_index_upper, Lambda_0.hat.s, dLambda_0.hat.s, V2.hat.inv, pRao_list, kappa_i.j.s
    )
    get_term <- function(i) {
      force(i)
      as.vector((m[i] / exp_X_beta_d_Lambda_0.n[i] - exp(W.bar[i,] %*% gamma.bar.hat)) * W.bar[i,] / 2)
    }
    retval <- list()
    for(i in seq_len(obj@n)) {
      retval[[i]] <- list()
      for(j in seq_len(obj@n)) {
        retval[[i]][[j]] <- apply(((- m * W.bar) / exp_X_beta_d_Lambda_0.n^2) * X_V_2_inv_g_ij_exp_X_beta_d_Lambda_0.n.n[[i]][[j]], 2, mean) + get_term(i) + get_term(j)
      }
    }
    obj@cache[[key]] <- retval
  }
  obj@cache[[key]]
}

# dxi_i.j.hat.gen <- function(obj) {
#   key <- "dxi_i.j.hat"
#   for(!is)
# }