# Implementation of Huang2004

#'@title \eqn{\hat{Z}_i}
#'
#'@param obj
#'@param F.hat
#'@param gamma.hat numeric vector, the \eqn{\hat{\alpha}} in Huang2004.
#'@reference Huang2004
Zi.hat <- function(obj, F.hat = NULL, gamma.hat = NULL, ...) {
  if (is.null(F.hat)) F.hat <- obj$F.hat
  if (is.null(gamma.hat)) gamma.hat <- obj$U.hat()[-1]
  m <- sapply(obj@t, length)
  F.hat.y.inv <- F.hat.y.inv(obj, F.hat)
  m * F.hat.y.inv / exp(obj@X[,-1] %*% gamma.hat)
}

#'@title Formulat(3) in Huang2004
#'@param obj recurrent-data object
#'@param Zi numeric vector, \eqn{\hat{Z}_i}
#'@param ... arguments pass to \code{\link{Zi.hat}}
#'@return numeric vector with length \code{ncol(obj@X) - 1}
U.hat.gen <- function(obj, Zi = NULL, ...) {
  if (is.null(Zi)) Zi <- Zi.hat(obj, ...)
  X <- obj@X[,-1]
  y <- obj@y
  n <- length(y)
  term_1 <- obj@D %*% X / n
  indicator.T <- outer(1:n, 1:n, function(i, j) as.integer(y[i] >= y[j]))
  function(beta) {
#     term_2_dem.fun <- function(i, j) {
#       Zi[j] * exp(X[j,] %*% beta) * as.integer(y[j] >= y[i])
#     }
    term_2_dem.mat <- t(as.vector(Zi * exp(X %*% beta)) * indicator.T)
    term_2_dem.i <- as.vector(term_2_dem.mat %*% rep(1, n)) #apply(term_2_dem.mat, 1, sum)
    term_2_num.i <- term_2_dem.mat %*% X # row vector is num of term_2 of U.hat(beta)
    term_2 <- obj@D %*% (term_2_num.i / term_2_dem.i) / n
    as.vector(term_1 - term_2)
  }
}

Gamma.hat.gen <- function(obj, Zi = NULL, ...) {
  if (is.null(Zi)) Zi <- Zi.hat(obj, ...)
  X <- obj@X[,-1]
  p <- ncol(X)
  y <- obj@y
  n <- length(y)
  indicator.T <- outer(1:n, 1:n, function(i, j) as.integer(y[i] >= y[j]))
  function(beta) {
    term_dem.mat <- t(as.vector(Zi * exp(X %*% beta)) * indicator.T)
    term_dem.i <- as.vector(term_dem.mat %*% rep(1, n))#apply(term_dem.mat, 1, sum)
    term_1 <- matrix(sapply(1:n, function(i) t(X) %*% (term_dem.mat[i,] * X)) %*% (obj@D/term_dem.i), p, p)
    term_2_num.i.tmp <- term_dem.mat %*% X # row vector is num of term_2 of U.hat(beta)
    term_2 <- t(term_2_num.i.tmp) %*% ((obj@D/term_dem.i^2) * term_2_num.i.tmp)
    (-term_1 + term_2) / n
  }  
}

#'@title Borrow Strength Method
#'
#'@param obj recurrent-data object
#'@param U.hat function, \eqn{\hat{U}(\beta)} in Huang2004
#'@param Gamma.hat function, \eqn{\hat{Gamma}(\beta) := \frac{d\hat{U}}{d\beta}} in Huang2004
#'@param Zi numeric vector, \eqn{\hat{Z}_i} in Huang2004
#'@param F.hat function, \eqn{\hat{\Lambda_0}} in Huang2004
#'@param gamma.hat numeric vecotr, \eqn{\hat{alpha}} in Huang2004 and \eqn{\hat{\gamma}} in WQC2001
#'@param verbose boolean value, whether print the message of \code{nleqslv} or not. Please see \code{\link{nleqslv}}.
#'@param tol numeric value, used to determine convergence. Please see \code{\link{nleqslv}}.
#'@param ...
#'@return numeric vector, \eqn{\hat{\beta}} in Huang2004
#'@export
BorrowStrengthMethod <- function(obj, U.hat = NULL, Gamma.hat = NULL, Zi = NULL, F.hat = NULL, gamma.hat = NULL, verbose = FALSE, tol = 1e-4, ...) {
  if (is.null(U.hat)) {
    if (is.null(Zi)) {
      Zi <- Zi.hat(obj, F.hat, gamma.hat)
    }
    U.hat <- U.hat.gen(obj, Zi)    
  }
  if (is.null(Gamma.hat)) {
    if (is.null(Zi)) {
      Zi <- Zi.hat(obj, F.hat, gamma.hat)
    }
    Gamma.hat <- Gamma.hat.gen(obj, Zi)  
  }
  temp <- nleqslv(rep(0, ncol(obj@X) - 1), U.hat, jac=Gamma.hat)
  if(verbose) {
    cat(sprintf("Check if beta is solved correctly: %s \n", paste(U.hat(temp$x), collapse=",")))
    cat(sprintf("message of nleqslv: %s ", temp$message))
  }
  if (sum(abs(U.hat(temp$x))) > tol) stop("Failed to converge during solving beta")
  return(temp$x)
}

H0.hat.gen <- function(obj, beta = NULL, Zi = NULL, F.hat = NULL, gamma.hat = NULL, ...) {
  if (is.null(Zi)) Zi <- Zi.hat(obj, F.hat = F.hat, gamma.hat = gamma.hat, ...)
  if (is.null(beta)) beta <- BorrowStrengthMethod(obj, Zi=Zi)
  X <- obj@X[,-1]
  y <- obj@y
  n <- length(y)
  indicator.T <- outer(1:n, 1:n, function(i, j) as.integer(y[i] >= y[j]))
  function(t) {
    term_dem.mat <- t(as.vector(Zi * exp(X %*% beta)) * indicator.T)
    term_dem.i <- as.vector(term_dem.mat %*% rep(1, n))
    as.vector((obj@D / term_dem.i) %*% outer(seq_along(y), seq_along(t), function(i, j) as.integer(y[i] <= t[j])))
  }
}

#'@title \eqn{\psi_{3i}}
#'@return numeric matrix, \code{retval[j,i]} = \eqn{psi_{3i}(y[j], b)}
psi_3i.y.gen <- function(obj, b, F.hat.y.inv = NULL, F.hat = NULL, 
                         gamma.hat.origin = NULL, fi.seq.origin = NULL, 
                         bi = NULL, bi.y = NULL, Zi = NULL,
                         indicator.T = NULL, ...) {
#   browser()
  m <- sapply(obj@t, length)
  n <- length(m)
  X <- obj@X[,-1]
  y <- obj@y
  if (is.null(F.hat)) F.hat <- obj$F.hat
  if (is.null(F.hat.y.inv)) F.hat.y.inv <- F.hat.y.inv(obj, F.hat)
  if (is.null(gamma.hat.origin)) gamma.hat.origin <- obj$U.hat()
  gamma.hat <- gamma.hat.origin[-1]
  if (is.null(bi)) {
    b.hat <- b.hat.gen(obj)
    bi <- lapply(1:n, function(i) Vectorize(b.hat(i)))
  }
  #bi.y[i,j] == bi[[j]](obj@y[i])
  if (is.null(bi.y)) bi.y <- sapply(1:n, function(i) bi[[i]](obj@y)) 
  if (is.null(fi.seq.origin)) fi.seq.origin <- fi.hat(obj, gamma = gamma.hat.origin, ...)
  fi.seq <- fi.seq.origin[-1,]
  if (is.null(Zi)) Zi <- Zi.hat(obj, F.hat, gamma.hat = gamma.hat, ...)
  m.F.hat.y.inv <- m * F.hat.y.inv
  if (is.null(indicator.T)) indicator.T <- outer(1:n, 1:n, function(i, j) as.integer(y[i] >= y[j]))
  # term1.2.i <- function(i) as.vector(X %*% fi.seq[,i]) + bi[[i]](obj@y)
  # stopifnot(all.equal(sapply(1:n, term1.2.i), X %*% fi.seq + bi.y))
  # term1.2[,1] == term1.2.i(1)
  term1.2 <- X %*% fi.seq + bi.y
  exp.X.b.gamma.hat <- exp(as.vector(X %*% (b - gamma.hat)))
  term2 <- term1.1 <- t(m.F.hat.y.inv * exp.X.b.gamma.hat * indicator.T)
  # all.equal(as.vector(as.vector(term1.2[,1]) %*% t(term1.1) / n), term1[,1])
  term1 <- term1.1 %*% term1.2 / n
#   term2.i <- function(i) m.F.hat.y.inv[i] * exp.X.b.gamma.hat[i] * indicator.T[i,]
#   all.equal(term2.i(3), term2[3,])
  term3 <- as.vector(as.vector(Zi * exp(X %*% b)) %*% indicator.T) / n
  retval <- term1 + term2 - term3 
  return(retval)
#   browser()
#   function(i) {
#     term1.2 <- as.vector(X %*% fi.seq[,i]) + bi.y[[i]] # independent of t
#     # term1.1.t<- function(t) m.F.hat.y.inv * exp(as.vector(X %*% (b - gamma.hat))) * as.integer(y >= t)
#     exp.X.b.gamma.hat <- exp(as.vector(X %*% (b - gamma.hat)))
#     term1.1 <- m.F.hat.y.inv * exp.X.b.gamma.hat * indicator.T
#     # stopifnot(sapply(1:n, function(i) all(term1.1[,i] == term1.1.t(obj@y[i]))))  
#     term1 <- as.vector(term1.2 %*% term1.1)/n
#     term2 <- m.F.hat.y.inv[i] * exp.X.b.gamma.hat[i] * indicator.T[i,]
# #       term3.t <- Vectorize(function(t) mean(Zi * exp(as.vector(X %*% b)) * as.integer(y >= t)))
#     term3 <- as.vector(as.vector(Zi * exp(X %*% b)) %*% indicator.T) / n
# #       stopifnot(isTRUE(all.equal(term3.t(y), term3)))
#     term1 + term2 - term3
#   }
}

#! retval[j,i,k]
psi_4i.y.gen <- function(obj, b, F.hat.y.inv = NULL, F.hat = NULL, 
                         gamma.hat.origin = NULL, fi.seq.origin = NULL, 
                         bi = NULL, bi.y = NULL, Zi = NULL,
                         indicator.T = NULL, ...) {
  m <- sapply(obj@t, length)
  n <- length(m)
  X <- obj@X[,-1]
  y <- obj@y
  if (is.null(F.hat)) F.hat <- obj$F.hat
  if (is.null(F.hat.y.inv)) F.hat.y.inv <- F.hat.y.inv(obj, F.hat)
  if (is.null(gamma.hat.origin)) gamma.hat.origin <- obj$U.hat()
  gamma.hat <- gamma.hat.origin[-1]
  if (is.null(bi)) {
    b.hat <- b.hat.gen(obj)
    bi <- lapply(1:n, function(i) Vectorize(b.hat(i)))
  }
  #bi.y[i,j] == bi[[j]](obj@y[i])
  if (is.null(bi.y)) bi.y <- sapply(1:n, function(i) bi[[i]](obj@y)) 
  if (is.null(fi.seq.origin)) fi.seq.origin <- fi.hat(obj, gamma = gamma.hat.origin, ...)
  fi.seq <- fi.seq.origin[-1,]
  if (is.null(Zi)) Zi <- Zi.hat(obj, F.hat, gamma.hat = gamma.hat, ...)
  m.F.hat.y.inv <- m * F.hat.y.inv
  if (is.null(indicator.T)) indicator.T <- indicator.T <- outer(1:n, 1:n, function(i, j) as.integer(y[i] >= y[j]))
  term1.2 <- X %*% fi.seq + bi.y
  exp.X.b.gamma.hat <- exp(as.vector(X %*% (b - gamma.hat)))
  retval.k <- function(k) {
    term2 <- term1.1 <- t((m.F.hat.y.inv * exp.X.b.gamma.hat * X[,k]) * indicator.T)
    term1 <- term1.1 %*% term1.2 / n
    term3 <- as.vector(as.vector(Zi * X[,k] * exp(X %*% b)) %*% indicator.T) / n
    retval <- term1 + term2 - matrix(term3, nrow=n, ncol=n)
    retval
  }
  retval <- sapply(1:ncol(X), retval.k, simplify="array")
  retval
}

psi_i.y.gen <- function(obj, b, F.hat.y.inv = NULL, F.hat = NULL, 
                         gamma.hat.origin = NULL, fi.seq.origin = NULL, 
                         bi = NULL, bi.y = NULL, Zi = NULL,
                         indicator.T = NULL, ...) {
  m <- sapply(obj@t, length)
  n <- length(m)
  X <- obj@X[,-1]
  y <- obj@y
  if (is.null(F.hat)) F.hat <- obj$F.hat
  if (is.null(F.hat.y.inv)) F.hat.y.inv <- F.hat.y.inv(obj, F.hat)
  if (is.null(gamma.hat.origin)) gamma.hat.origin <- obj$U.hat()
  gamma.hat <- gamma.hat.origin[-1]
  if (is.null(bi)) {
    b.hat <- b.hat.gen(obj)
    bi <- lapply(1:n, function(i) Vectorize(b.hat(i)))
  }
  #bi.y[i,j] == bi[[j]](obj@y[i])
  if (is.null(bi.y)) bi.y <- sapply(1:n, function(i) bi[[i]](obj@y)) 
  if (is.null(fi.seq.origin)) fi.seq.origin <- fi.hat(obj, gamma = gamma.hat.origin, ...)
  fi.seq <- fi.seq.origin[-1,]
  if (is.null(Zi)) Zi <- Zi.hat(obj, F.hat, gamma.hat = gamma.hat, ...)
  m.F.hat.y.inv <- m * F.hat.y.inv
  if (is.null(indicator.T)) indicator.T <- indicator.T <- outer(1:n, 1:n, function(i, j) as.integer(y[i] >= y[j]))
  psi_3i.y <- psi_3i.y.gen(obj, b, F.hat.y.inv = F.hat.y.inv, F.hat = F.hat, 
                           gamma.hat.origin = gamma.hat.origin, 
                           fi.seq.origin = fi.seq.origin, 
                           bi = bi, bi.y = bi.y, Zi = Zi,
                           indicator.T = indicator.T, ...)
  psi_4i.y <- psi_4i.y.gen(obj, b, F.hat.y.inv = F.hat.y.inv, F.hat = F.hat, 
                           gamma.hat.origin = gamma.hat.origin, 
                           fi.seq.origin = fi.seq.origin, 
                           bi = bi, bi.y = bi.y, Zi = Zi,
                           indicator.T = indicator.T, ...)
  term2 <- obj@D %*% X / n
  Z.exp.X.beta.I.Y.geq.s <- t(as.vector(Zi * exp(X %*% b)) * indicator.T)
  term6.num <- term5.num <- term3.num.2 <- Z.exp.X.beta.I.Y.geq.s %*% X
  term6.dem <- term5.dem <- term4.dem <- term3.dem <- as.vector(Z.exp.X.beta.I.Y.geq.s %*% rep(1, n))
  term6 <- obj@D %*% (term6.num / term6.dem) / n
  retval.i <- function(i) {
    term1 <- obj@D[i] * X[i,]
    term3.num.1 <- psi_3i.y[,i]
    term3 <- obj@D %*% (term3.num.1 * term3.num.2 / term3.dem^2) #! no need to divide n because term3.dem^2 contains n^2
    term4.num <- psi_4i.y[,i,]
    term4 <- obj@D %*% (term4.num / term4.dem) #!
    term5 <- obj@D[i] * term5.num[i,] / term5.dem[i]
    as.vector(term1 - term2 + term3 - term4 - term5 + term6)
  }
  t(sapply(1:n, retval.i))
}

Sigma.hat <- function(obj, b, ...) {
  psi_i.y <- psi_i.y.gen(obj, b, ...)
  n <- nrow(psi_i.y)
  var(psi_i.y) * (n-1) / n
}

phi_i.y.gen <- function(obj, b, F.hat.y.inv = NULL, F.hat = NULL, 
                        gamma.hat.origin = NULL, fi.seq.origin = NULL, 
                        bi = NULL, bi.y = NULL, Zi = NULL, Gamma.hat = NULL,
                        indicator.T = NULL, ...) {
  m <- sapply(obj@t, length)
  n <- length(m)
  X <- obj@X[,-1]
  y <- obj@y
  if (is.null(F.hat)) F.hat <- obj$F.hat
  if (is.null(F.hat.y.inv)) F.hat.y.inv <- F.hat.y.inv(obj, F.hat)
  if (is.null(gamma.hat.origin)) gamma.hat.origin <- obj$U.hat()
  gamma.hat <- gamma.hat.origin[-1]
  if (is.null(bi)) {
    b.hat <- b.hat.gen(obj)
    bi <- lapply(1:n, function(i) Vectorize(b.hat(i)))
  }
  #bi.y[i,j] == bi[[j]](obj@y[i])
  if (is.null(bi.y)) bi.y <- sapply(1:n, function(i) bi[[i]](obj@y)) 
  if (is.null(fi.seq.origin)) fi.seq.origin <- fi.hat(obj, gamma = gamma.hat.origin, ...)
  fi.seq <- fi.seq.origin[-1,]
  if (is.null(Zi)) Zi <- Zi.hat(obj, F.hat, gamma.hat = gamma.hat, ...)
  if (is.null(Gamma.hat)) Gamma.hat <- Gamma.hat.gen(obj, Zi) 
  m.F.hat.y.inv <- m * F.hat.y.inv
  if (is.null(indicator.T)) indicator.T <- indicator.T <- outer(1:n, 1:n, function(i, j) as.integer(y[i] >= y[j]))
  psi_3i.y <- psi_3i.y.gen(obj, b, F.hat.y.inv = F.hat.y.inv, F.hat = F.hat, 
                           gamma.hat.origin = gamma.hat.origin, 
                           fi.seq.origin = fi.seq.origin, 
                           bi = bi, bi.y = bi.y, Zi = Zi,
                           indicator.T = indicator.T, ...)
  Z.exp.X.beta.I.Y.geq.s <- t(as.vector(Zi * exp(X %*% b)) * indicator.T)
  psi_i.y <- psi_i.y.gen(obj, b, F.hat.y.inv = F.hat.y.inv, F.hat = F.hat, 
                         gamma.hat.origin = gamma.hat.origin, 
                         fi.seq.origin = fi.seq.origin, bi = bi, bi.y = bi.y, 
                         Zi = Zi, indicator.T = indicator.T, ...)
#   term6.num <- term5.num <- term3.num.2 <- Z.exp.X.beta.I.Y.geq.s %*% X
#   term6.dem <- term5.dem <- term4.dem <- term3.dem <- as.vector(Z.exp.X.beta.I.Y.geq.s %*% rep(1, n))
#   term6 <- obj@D %*% (term6.num / term6.dem) / n
#   term5 <- obj@D[i] * term5.num[i,] / term5.dem[i]
  term4.num.1 <- Z.exp.X.beta.I.Y.geq.s %*% X
  term4.dem.1 <- term3.dem <- term2.dem <- term1.dem <- as.vector(Z.exp.X.beta.I.Y.geq.s %*% rep(1, n))
  #     obj@D * indicator.T[1,] == (obj@D * t(indicator.T))[,1]
  #     obj@D * indicator.T[2,] == (obj@D * t(indicator.T))[,2]
  #     obj@D * indicator.T[3,] == (obj@D * t(indicator.T))[,3]
  D.indicator.T <- t(obj@D * t(indicator.T))
  retval.i <- function(i) { # t = y[j] 
    term1.num <- psi_3i.y[,i]
#     y.t <- ifelse(y <= y[3], 1, 0) # indicator.T[3,] == y.t
#     obj@D * y.t == D.indicator.T[3,]
#     term1 <- n * (obj@D * y.t) %*% (term1.num / term1.dem^2) #! term3.dem^2 contains n^2
    term1 <- as.vector(n * D.indicator.T %*% (term1.num / term1.dem^2))
#     term2 == sapply(y, function(t) {
#       y.t <- ifelse(y <= t, 1, 0) # D.indicator.T[j,]
#       n * obj@D[i] * y.t[i] / term2.dem[i]
#     })
    term2 <- n * obj@D[i] * indicator.T[,i] / term2.dem[i]
#     term3[3] == (obj@D * y.t) %*% (1 / term3.dem)
    term3 <- as.vector(D.indicator.T %*% (1 / term3.dem))
    #
#     term4.1 <- (obj@D * y.t) %*% (term4.num.1 / term4.dem.1^2) # == term4.1[3,]
#     term4.2 <- solve(Gamma.hat(b))
#     term4.3 <- psi_i.y[3,]
#     term4.1 %*% term4.2 %*% term4.3
    #
    term4.1 <- D.indicator.T %*% (term4.num.1 / term4.dem.1^2)
    term4.2 <- solve(Gamma.hat(b))
    term4.3 <- t(psi_i.y)
    term4 <- diag((term4.1 %*% term4.2) %*% term4.3)
    term1 + term2 + term3 + term4
  }  
  retval <- sapply(1:n, retval.i) # retval[,i] -- retval.i(i)
  retval
}

#'@export
Huang2004 <- function(obj) {
  n <- length(obj@y)
  F.hat <- obj$F.hat
  gamma.hat.origin <- obj$U.hat()
  Zi <- Zi.hat(obj, F.hat)
  Gamma.hat <- Gamma.hat.gen(obj, Zi)
  beta.hat <- BorrowStrengthMethod(obj, Gamma.hat = Gamma.hat, Zi = Zi, 
                                   F.hat = F.hat, gamma.hat = gamma.hat.origin[-1])
  H0.hat <- H0.hat.gen(obj, beta=beta.hat, Zi=Zi)
  Sigma <- Sigma.hat(obj, beta.hat)
  Gamma <- solve(Gamma.hat(beta.hat))
  H0.hat.y <- sapply(obj@y, H0.hat)
  phi_i.y <- phi_i.y.gen(obj, beta.hat)
  H0.hat.y.sd <- apply(phi_i.y / n, 1, sd)
  i.y <- order(y)
  return(list(beta.hat=beta.hat, y = y[i.y], H0.hat.y=H0.hat.y[i.y], H0.hat.y.sd=H0.hat.y.sd[i.y], beta.hat.var=(Gamma %*% Sigma %*% t(Gamma)) / n))
}