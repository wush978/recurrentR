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

d_beta.gen <- function(obj) {
  key <- "d_beta"
  if (!is_cache(obj, key)) {
    beta.hat <- beta.hat.gen(obj)
    X.value <- X.value.gen(obj)
    t_index <- t_index.gen(obj)
    obj@cache[[key]] <- d_beta(beta.hat, X.value, t_index)
  }
  obj@cache[[key]]
}

R_beta.gen <- function(obj) {
  key <- "R_beta"
  if (!is_cache(obj, key)) {
    beta.hat <- beta.hat.gen(obj)
    X.value <- X.value.gen(obj)
    t_index <- t_index.gen(obj)
    y_vs_s <- outer(obj@y, obj@s, ">=")
    s_index_upper <- apply(y_vs_s, 1, function(a) {
      retval <- which(!a)
      ifelse(length(retval) == 0, length(obj@s), min(retval) - 1)
    })
    obj@cache[[key]] <- R_beta(beta.hat, X.value, t_index, s_index_upper)
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
