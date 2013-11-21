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

psi_3i.y.gen <- function(obj, b) {
  if (!is_cache(obj, "psi_3i.y")) {
    m <- m.gen(obj)
    n <- obj@n
    Lambda_0.hat.y.inv <- Lambda_0.hat.y.inv.gen(obj)
    gamma.hat <- gamma.hat.gen(obj)
    bi.y <- b.hat.y.gen(obj)
    fi.hat.i <- recurrentR:::fi.hat.i.gen(obj)
    fi.seq <- fi.hat.i[-1, , drop = FALSE]
    Zi <- Z.hat.gen(obj)
    m.F.hat.y.inv <- m * Lambda_0.hat.y.inv
    indicator.T <- indicator.T.gen(obj)
    term1.2 <- obj@W %*% fi.seq + bi.y
    exp.X.b.gamma.hat <- exp(as.vector(obj@W %*% (b - gamma.hat)))
    term2 <- term1.1 <- t(m.F.hat.y.inv * exp.X.b.gamma.hat * indicator.T)
    term1 <- term1.1 %*% term1.2 / n
    term3 <- as.vector(as.vector(Zi * exp(obj@W %*% b)) %*% indicator.T) / n
    retval <- term1 + term2 - term3 
    obj@cache[["psi_3i.y"]] <- retval
  }
  obj@cache[["psi_3i.y"]]
}

psi_4i.y.gen <- function(obj, b) {
  if (!is_cache(obj, "psi_4i.y")) {
    m <- m.gen(obj)
    Lambda_0.hat.y.inv <- Lambda_0.hat.y.inv.gen(obj)
    gamma.hat <- gamma.hat.gen(obj)
    bi.y <- b.hat.y.gen(obj)
    fi.hat.i <- recurrentR:::fi.hat.i.gen(obj)
    fi.seq <- fi.hat.i[-1, , drop = FALSE]
    Zi <- Z.hat.gen(obj)
    m.F.hat.y.inv <- m * Lambda_0.hat.y.inv
    indicator.T <- indicator.T.gen(obj)
    term1.2 <- obj@W %*% fi.seq + bi.y
    exp.X.b.gamma.hat <- exp(as.vector(obj@W %*% (b - gamma.hat)))
    retval.k <- function(k) {
      term2 <- term1.1 <- t((m.F.hat.y.inv * exp.X.b.gamma.hat * obj@W[,k]) * indicator.T)
      term1 <- term1.1 %*% term1.2 / obj@n
      term3 <- as.vector(as.vector(Zi * obj@W[,k] * exp(obj@W %*% b)) %*% indicator.T) / obj@n
      retval <- term1 + term2 - matrix(term3, nrow=obj@n, ncol=obj@n)
      retval
    }
    retval <- sapply(1:ncol(obj@W), retval.k, simplify="array")
    obj@cache[["psi_4i.y"]] <- retval
  }
  obj@cache[["psi_4i.y"]]
}

psi_i.y.gen <- function(obj, b) {
  if (!is_cache(obj, "psi_i.y")) {
    m <- m.gen(obj)
    Lambda_0.hat.y.inv <- Lambda_0.hat.y.inv.gen(obj)
    gamma.hat <- gamma.hat.gen(obj)
    bi.y <- b.hat.y.gen(obj)
    fi.hat.i <- recurrentR:::fi.hat.i.gen(obj)
    fi.seq <- fi.hat.i[-1, , drop = FALSE]
    Zi <- Z.hat.gen(obj)
    m.F.hat.y.inv <- m * Lambda_0.hat.y.inv
    indicator.T <- indicator.T.gen(obj)
    psi_3i.y <- psi_3i.y.gen(obj, b)
    psi_4i.y <- psi_4i.y.gen(obj, b)
    term2 <- obj@D %*% obj@W / obj@n
    Z.exp.X.beta.I.Y.geq.s <- t(as.vector(Zi * exp(obj@W %*% b)) * indicator.T)
    term6.num <- term5.num <- term3.num.2 <- Z.exp.X.beta.I.Y.geq.s %*% obj@W
    term6.dem <- term5.dem <- term4.dem <- term3.dem <- as.vector(Z.exp.X.beta.I.Y.geq.s %*% rep(1, obj@n))
    term6 <- obj@D %*% (term6.num / term6.dem) / obj@n
    retval.i <- function(i) {
      term1 <- obj@D[i] * obj@W[i,]
      term3.num.1 <- psi_3i.y[,i]
      term3 <- obj@D %*% (term3.num.1 * term3.num.2 / term3.dem^2) #! no need to divide n because term3.dem^2 contains n^2
      term4.num <- psi_4i.y[,i,]
      term4 <- obj@D %*% (term4.num / term4.dem) #!
      term5 <- obj@D[i] * term5.num[i,] / term5.dem[i]
      as.vector(term1 - term2 + term3 - term4 - term5 + term6)
    }
    obj@cache[["psi_i.y"]] <- if (ncol(obj@W) > 1) t(sapply(1:obj@n, retval.i)) else matrix(sapply(1:obj@n, retval.i), ncol=1)
  }
  obj@cache[["psi_i.y"]]
}

Sigma.hat.gen <- function(obj, b) {
  key <- "Sigma.hat"
  if (!is_cache(obj, key)) {
    psi_i.y <- psi_i.y.gen(obj, b)
    obj@cache[[key]] <- var(psi_i.y) * (obj@n-1) / obj@n
  }
  obj@cache[[key]]
}

phi_i.y.gen <- function(obj, b) {
  key <- "phi_i.y"
  if (!is_cache(obj, key)) {
    m <- m.gen(obj)
    Lambda_0.hat.y.inv <- Lambda_0.hat.y.inv.gen(obj)
    gamma.hat <- gamma.hat.gen(obj)
    bi.y <- b.hat.y.gen(obj)
    fi.hat.i <- recurrentR:::fi.hat.i.gen(obj)
    fi.seq <- fi.hat.i[-1, , drop = FALSE]
    Zi <- Z.hat.gen(obj)
    m.F.hat.y.inv <- m * Lambda_0.hat.y.inv
    indicator.T <- indicator.T.gen(obj)
    psi_3i.y <- psi_3i.y.gen(obj, b)
    psi_4i.y <- psi_4i.y.gen(obj, b)
    Gamma.hat <- Gamma.gen(obj) 
    indicator.T <- indicator.T.gen(obj)
    psi_3i.y <- psi_3i.y.gen(obj, b)
    Z.exp.X.beta.I.Y.geq.s <- t(as.vector(Zi * exp(obj@W %*% b)) * indicator.T)
    psi_i.y <- psi_i.y.gen(obj, b)
    term4.num.1 <- Z.exp.X.beta.I.Y.geq.s %*% obj@W
    term4.dem.1 <- term3.dem <- term2.dem <- term1.dem <- as.vector(Z.exp.X.beta.I.Y.geq.s %*% rep(1, obj@n))
    D.indicator.T <- t(obj@D * t(indicator.T))
    retval.i <- function(i) {
      term1.num <- psi_3i.y[,i]
      term1 <- as.vector(obj@n * D.indicator.T %*% (term1.num / term1.dem^2))
      term2 <- obj@n * obj@D[i] * indicator.T[,i] / term2.dem[i]
      term3 <- as.vector(D.indicator.T %*% (1 / term3.dem))
      term4.1 <- D.indicator.T %*% (term4.num.1 / term4.dem.1^2)
      term4.2 <- solve(Gamma.hat(b))
      term4.3 <- t(psi_i.y) # 2 * 100
      term4 <- diag((term4.1 %*% term4.2) %*% term4.3)
      term1 + term2 + term3 + term4
    }  
    obj@cache[[key]] <- sapply(1:obj@n, retval.i)
  }
  obj@cache[[key]]
}

#'@export
Huang2004 <- function(obj) {
  y <- obj@y
  Lambda_0.hat <- Lambda_0.hat.gen(obj)
  gamma.hat <- gamma.hat.gen(obj)
  Zi <- Z.hat.gen(obj)
  Gamma.hat <- Gamma.gen(obj)
  alpha.hat <- alpha.hat.gen(obj)
  H0.hat <- H_0.hat.gen(obj)
  Sigma <- Sigma.hat.gen(obj, alpha.hat)
  Gamma <- solve(Gamma.gen(obj)(alpha.hat))
  H0.hat.y <- sapply(obj@y, H0.hat)
  phi_i.y <- phi_i.y.gen(obj, alpha.hat)
  H0.hat.y.sd <- apply(phi_i.y / obj@n, 1, sd)
  i.y <- order(y)
  return(list(alpha.hat=alpha.hat, y = y[i.y], H0.hat.y=H0.hat.y[i.y], H0.hat.y.sd=H0.hat.y.sd[i.y], alpha.hat.var=(Gamma %*% Sigma %*% t(Gamma)) / obj@n))
  
}