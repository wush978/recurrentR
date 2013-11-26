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
    pRao <- pRao.gen(obj)
    X.value <- X.value.gen(obj)
  }
  obj@cache[[key]]
}