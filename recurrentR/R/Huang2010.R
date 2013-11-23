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

pRao.gen.array <- function(obj, X.value) {
  stopifnot(length(dim(X.value)) == 3)
  key <- "pRao"
  if (!is_cache(obj, key)) {
    m <- m.gen(obj)
    obj@cache[[key]] <- rao_gen_array(X.value, m, obj@t, obj@y, obj@s)
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