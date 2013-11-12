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
    n = "integer",
    eval = "numeric",
    s = "numeric",
    d = "integer",
    cache = "environment"
    ),
  validity = function(object) {
    stopifnot(length(object@T_0) == 1)
    stopifnot(length(object@y) == length(object@D))
    stopifnot(length(object@y) == nrow(object@W))
    if (length(object@X) > 0) stopifnot(length(object@y) == length(object@X))
  }
)

#'@export
create_recurrent_data <- function(...) {
  UseMethod("create_recurrent_data")
}

#'@export
create_recurrent_data.numeric <- function(y, D, t, T_0, W) {
  retval <- new("recurrent-data")
  retval@y <- y
  retval@D <- D
  retval@t <- t
  retval@T_0 <- T_0
  retval@W <- W
  retval@n <- length(y)
  validObject(retval)
  retval@s <- s(retval)
  retval@d <- d(retval)
  retval
}

#'@export
create_recurrent_data.list <- function(X, y, D, t, T_0, W) {
  retval <- create_recurrent_data.numeric(y, D, t, T_0, W)
  retval@X <- X
  validObject(retval)
  retval
}

s <- function(obj) {
  sort(unique(unlist(sapply(obj@t, unique))))
}

d <- function(obj) {
  as.vector(table(unlist(obj@t)))
}
