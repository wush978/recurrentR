Lambda_0.hat.gen <- function(obj) {
  s <- obj@s
  d <- obj@d
  N <- cumsum(d)
  y.i <- order(obj@y)
  m <- sapply(obj@t, length)
  N <- N + eval_N(s, obj@y[y.i], m[y.i])
  x <- s
  y <- append(rev(cumprod(rev(1 - d/N))), 1)
  f <- stepfun(x, y)
  retval <- function(t, bootstrap = FALSE, B = 100, error.measurement.function = stats::sd) {
    if (!bootstrap) return(f(t))
    estimate <- obj$F.hat(t)
    obj.b <- obj$bootstrap(B)
    F.hat.bootstrap.result <- sapply(obj.b, function(obj) {
      obj$F.hat(t)
    })
    error.measurement = apply(F.hat.bootstrap.result, 1, error.measurement.function)
    list(estimate=estimate, error.measurement=error.measurement)
  }
  return(retval)
}