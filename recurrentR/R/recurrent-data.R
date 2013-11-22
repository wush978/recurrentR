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