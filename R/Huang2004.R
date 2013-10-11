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
  F.hat.y <- F.hat(obj@y)
  F.hat.y.inv <- 1/F.hat.y
  F.hat.y.inv[F.hat.y == 0] <- 0
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
    term_2_dem.i <- apply(term_2_dem.mat, 1, sum)
    term_2_num.i <- term_2_dem.mat %*% X # row vector is num of term_2 of U.hat(beta)
    term_2 <- obj@D %*% (term_2_num.i / term_2_dem.i) / n
    as.vector(term_1 - term_2)
  }
}

Gamma.hat.gen <- function(obj, Zi = NULL, ...) {
  if (is.null(Zi)) Zi <- Zi.hat(obj, ...)
  X <- obj@X[,-1]
  y <- obj@y
  n <- length(y)
  indicator.T <- outer(1:n, 1:n, function(i, j) as.integer(y[i] >= y[j]))
  function(beta) {
    term_dem.mat <- t(as.vector(Zi * exp(X %*% beta)) * indicator.T)
    term_dem.i <- apply(term_dem.mat, 1, sum)
    term_1 <- matrix(sapply(1:n, function(i)t(X) %*% (term_dem.mat[i,] * X)) %*% (obj@D/term_dem.i), 2, 2)
    term_2_num.i.tmp <- term_dem.mat %*% X # row vector is num of term_2 of U.hat(beta)
    term_2 <- t(term_2_num.i.tmp) %*% ((obj@D/term_dem.i^2) * term_2_num.i.tmp)
    (-term_1 + term_2) / n
  }  
}

BorrowStrengthMethod <- function(obj, U.hat = NULL, Gamma.hat = NULL, Zi = NULL, F.hat = NULL, gamma.hat = NULL, verbose = FALSE, tol = 1e-4, ...) {
  if (is.null(U.hat)) {
    if (is.null(Zi)) {
      if (is.null(F.hat)) F.hat <- obj$F.hat
      if (is.null(gamma.hat)) gamma.hat <- obj$U.hat()[-1]
      Zi <- Zi.hat(obj, F.hat, gamma.hat)
    }
    U.hat <- U.hat.gen(obj, Zi)    
  }
  if (is.null(Gamma.hat)) {
    if (is.null(Zi)) {
      if (is.null(F.hat)) F.hat <- obj$F.hat
      if (is.null(gamma.hat)) gamma.hat <- obj$U.hat()[-1]
      Zi <- Zi.hat(obj, F.hat, gamma.hat)
    }
    Gamma.hat <- Gamma.hat.gen(obj, Zi)  
  }
  temp <- nleqslv(rep(0, ncol(obj@X) - 1), U.hat, jac=Gamma.hat)
  if(verbose) {
    cat(sprintf("Check if gamma is solved correctly: %s \n", paste(U.hat(temp$x), collapse=",")))
    cat(sprintf("message of nleqslv: %s ", temp$message))
  }
  if (sum(abs(U.hat(temp$x))) > tol) stop("Failed to converge during solving gamma")
  return(temp$x)
}