#'@title S4 Object of Recurrent Data
#'@param y numeric vector. \code{y[i]} is the failure or censored time of instance \code{i}.
#'@param D logical vector. \code{D[i]} indicates whether the 
#'corresponding \code{y[i]} is failure time(\code{TRUE}) or not.
#'@param T_0 positive numeric value. The time period of observations is \eqn{[0, T_0]}.
#'@param t list of numeric vector. The \code{t[[i]]} is the recurrent event time, \eqn{t_{i,1}$, $t_{i,2}$, 
#'..., $t_{i,m_i}} of instance \code{i}. There are several statistical model describes the behavior of 
#'recurrent event time. Please see \link{Wang2001}, \link{Huang2004} and \link{Huang2010} for details.
#'@param W numeric matrix. \code{W[i,]} is the time independent covariates of the instance \code{i}.
#'@param X list of functions. \code{X[[i]](t)} is the time dependent covariate process of instance \code{i}. 
#'This slot is only used in \link{Huang2010}.
#'@author Wush Wu
#'@exportClass "recurrent-data"
setClass(
  "recurrent-data",
  representation(
    y = "numeric",
    D = "logical",
    t = "list",
    T_0 = "numeric",
    W = "matrix",
    X = "list",
    X_dim = "integer",
    n = "integer",
    eval = "numeric",
    s = "numeric",
    d = "integer",
    cache = "environment",
    tol = "numeric"
    ),
  validity = function(object) {
    stopifnot(length(object@T_0) == 1)
    stopifnot(length(object@y) == length(object@D))
    stopifnot(length(object@y) == nrow(object@W))
    lapply(seq_along(object@y), function(i) {
      if (length(object@t[[i]]) > 0) {
        stopifnot(all(object@t[[i]] <= object@y[i]))
        stopifnot(object@y[i] <= object@T_0)
      }
    })
    if (length(object@X) > 0) {
      stopifnot(length(object@y) == length(object@X))
      temp <- lapply(object@X, function(f) f(0))
      stopifnot(sapply(temp, length) == object@X_dim)
    } else {
      stopifnot(object@X_dim == 0)
    }
  }
)

#'@export
create_recurrent_data <- function(...) {
  UseMethod("create_recurrent_data")
}

#'@export
create_recurrent_data.numeric <- function(y, D, t, T_0, W, tol = 1e-4) {
  retval <- new("recurrent-data")
  retval@y <- y
  retval@D <- D
  retval@t <- t
  retval@T_0 <- T_0
  # correct data
  for(i in seq_along(retval@y)) {
    if (retval@y[i] > retval@T_0) {
      retval@y[i] <- retval@T_0
      retval@D[i] <- FALSE
      retval@t[[i]] <- retval@t[[i]][which(retval@t[[i]] <= retval@T_0)]
    }
  }
  retval@W <- W
  retval@n <- length(y)
  retval@X_dim <- 0L
  validObject(retval)
  retval@s <- s(retval)
  retval@d <- d(retval)
  retval@tol <- tol
  retval@cache <- new.env()
  retval
}

#'@export
create_recurrent_data.list <- function(X, y, D, t, T_0, W, tol = 1e-4, X_dim = NULL) {
  retval <- create_recurrent_data.numeric(y, D, t, T_0, W, tol)
  retval@X <- X
  if (is.null(X_dim)) X_dim <- length(X[[1]](0))
  retval@X_dim <- X_dim
  validObject(retval)
  retval
}

s <- function(obj) {
  sort(unique(unlist(sapply(obj@t, unique))))
}

d <- function(obj) {
  as.vector(table(unlist(obj@t)))
}

#'@export
cache_clean <- function(obj) {
  rm(list=ls(obj@cache), envir=obj@cache)
}