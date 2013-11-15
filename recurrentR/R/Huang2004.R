Z.hat.gen <- function(obj) {
  if (!is_cache(obj, "Z.hat")) {
    m <- m.gen(obj)
    gamma.hat <- gamma.hat.gen(obj)
    Lambda_0.hat.y.inv <- Lambda_0.hat.y.inv.gen(obj)
    obj@cache[["Z.hat"]] <- as.vector(m * Lambda_0.hat.y.inv / exp(obj@W %*% gamma.hat))
  }
  obj@cache[["Z.hat"]]
}

indicator.T.gen <- function(obj) {
  if (!is_cache(obj, "indicator.T")) {
    y_k_vs_y_i <- y_k_vs_y_i.gen(obj)
    obj@cache[["indicator.T"]] <- t(y_k_vs_y_i)
  }
  obj@cache[["indicator.T"]]
}

U.gen <- function(obj) {
  if (!is_cache(obj, "U")) {
    indicator.T <- indicator.T.gen(obj)
    term_1 <- as.vector(obj@D %*% obj@W / obj@n)
    Z.hat <- Z.hat.gen(obj)
    obj@cache[["U"]] <- function(alpha) {
      W.effect <- exp(as.vector(obj@W %*% alpha))
      Z.W.effect <- Z.hat * W.effect
      term_2_dem.mat <- t(Z.W.effect * indicator.T)
      term_2_dem.i <- as.vector(term_2_dem.mat %*% rep(1, obj@n)) #apply(term_2_dem.mat, 1, sum)
      term_2_num.i <- term_2_dem.mat %*% obj@W # row vector is num of term_2 of U.hat(beta)
      term_2 <- obj@D %*% (term_2_num.i / term_2_dem.i) / obj@n
      as.vector(term_1 - term_2)
    }
  }
  obj@cache[["U"]]
}

Gamma.gen <- function(obj) {
  if (!is_cache(obj, "Gamma")) {
    Z.hat <- Z.hat.gen(obj)
    indicator.T <- indicator.T.gen(obj)
    obj@cache[["Gamma"]] <- function(alpha) {
      W.effect <- exp(as.vector(obj@W %*% alpha))
      Z.W.effect <- Z.hat * W.effect
      term_dem.mat <- t(Z.W.effect * indicator.T)
      term_dem.i <- as.vector(term_dem.mat %*% rep(1, obj@n))#apply(term_dem.mat, 1, sum)
      term_1 <- matrix(sapply(1:obj@n, function(i) t(obj@W) %*% (term_dem.mat[i,] * obj@W)) %*% (obj@D/term_dem.i), ncol(obj@W), ncol(obj@W))
      term_2_num.i.tmp <- term_dem.mat %*% obj@W # row vector is num of term_2 of U.hat(beta)
      term_2 <- t(term_2_num.i.tmp) %*% ((obj@D/term_dem.i^2) * term_2_num.i.tmp)
      (-term_1 + term_2) / obj@n
    }    
  }
  obj@cache[["Gamma"]]
}

alpha.hat.gen <- function(obj) {
  if (!is_cache(obj, "alpha.hat")) {
    U <- U.gen(obj)
    Gamma <- Gamma.gen(obj)
    slv <- nleqslv(rep(0, ncol(obj@W)), U, Gamma)
    if (sum(abs(U(slv$x))) > obj@tol) stop("Failed to converge during solving gamma")
    obj@cache[["alpha.hat"]] <- slv$x
  }
  obj@cache[["alpha.hat"]]
}

H_0.hat.gen <- function(obj) {
  if (!is_cache(obj, "H_0.hat")) {
    Z.hat <- Z.hat.gen(obj)
    alpha.hat <- alpha.hat.gen(obj)
    indicator.T <- indicator.T.gen(obj)
    W.effect <- exp(as.vector(obj@W %*% alpha.hat))
    term_dem.mat <- t(Z.hat * W.effect * indicator.T)
    term_dem.i <- as.vector(term_dem.mat %*% rep(1, obj@n))
    obj@cache[["H_0.hat"]] <- function(t) {
      as.vector((obj@D / term_dem.i) %*% outer(seq_len(obj@n), seq_along(t), function(i, j) as.integer(obj@y[i] <= t[j])))
    }
  }
  obj@cache[["H_0.hat"]]
}

psi_3i.y.gen <- function(obj) {
  if (!is_cache(obj, "psi_3i.y")) {
    browser()
  }
  obj@cache[["psi_3i.y"]]
}
