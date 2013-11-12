Lambda_0.hat.gen <- function(obj) {
  if (!exists("Lambda_0.hat", envir=obj@cache, inherits=FALSE)) {
    s <- obj@s
    d <- obj@d
    N <- cumsum(d)
    y.i <- order(obj@y)
    m <- sapply(obj@t, length)
    N <- N + eval_N(s, obj@y[y.i], m[y.i])
    x <- s
    y <- append(rev(cumprod(rev(1 - d/N))), 1)
    f <- stepfun(x, y)
    obj@cache[["Lambda_0.hat"]] <- function(t, bootstrap = FALSE, B = 100, error.measurement.function = stats::sd) {
      if (!bootstrap) return(f(t))
      stop("TODO: bootstrap")
    }
  }
  obj@cache[["Lambda_0.hat"]]
}

Lambda_0.hat.y.gen <- function(obj) {
  if (!exists("Lambda_0.hat(y)", envir=obj@cache, inherits=FALSE)) {
    Lambda_0.hat <- Lambda_0.hat.gen(obj)
    obj@cache[["Lambda_0.hat(y)"]] <- Lambda_0.hat(obj@y)
  }
  obj@cache[["Lambda_0.hat(y)"]]
}

Lambda_0.hat.y.inv.gen <- function(obj) {
  if (!exists("1/Lambda_0.hat(y)", envir=obj@cache, inherits=FALSE)) {
    Lambda_0.hat.y <- Lambda_0.hat.y.gen(obj)
    retval <- 1 / Lambda_0.hat.y
    retval[Lambda_0.hat.y == 0] <- 0
    obj@cache[["1/Lambda_0.hat(y)"]] <- retval
  }
  obj@cache[["1/Lambda_0.hat(y)"]]
}

m <- function(obj) {
  if (!exists("m", envir=obj@cache, inherits=FALSE)) {
    obj@cache[["m"]] <- sapply(obj@t, length)
  }
  obj@cache[["m"]]
}

gamma.hat.gen <- function(obj) {
  if (!exists("gamma.hat", envir=obj@cache, inherits=FALSE)) {
    b <- m(obj) * Lambda_0.hat.y.inv.gen(obj)
    W.bar <- cbind(1, obj@W)
    gamma.hat <- rep(0, ncol(W.bar))
    g <- function(gamma) t(W.bar) %*% (b - exp(W.bar %*% gamma))
    g.grad <- function(gamma) - t(W.bar) %*% diag(c(exp(W.bar %*% gamma)), nrow(W.bar), nrow(W.bar)) %*% W.bar
    slv <- nleqslv(gamma.hat, g, jac=g.grad)
    if (sum(abs(g(slv$x))) > obj@tol) stop("Failed to converge during solving gamma")
    obj@cache[["gamma.hat"]] <- slv$x
  }
  obj@cache[["gamma.hat"]]
}