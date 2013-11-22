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

S.gen <- function(obj) {
  key <- "S"
  if (!is_cache(obj, key)) {
    obj@cache[[key]] <- function(beta) {
      S(pRao.gen(obj), beta)
    }
  }
  obj@cache[[key]]
}

dS_over_dbeta.gen <- function(obj) {
  key <- "dS_over_dbeta"
  if (!is_cache(obj, key)) {
    obj@cache[[key]] <- function(beta) {
      dS_over_dbeta(pRao.gen(obj), beta)
    }
  }
  obj@cache[[key]]
}