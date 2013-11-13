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

Q.hat <- function(obj) {
  s <- obj@s
  y <- sapply(c(0, s), function(u) mean(sapply(obj@t, function(t) sum(t <= u))))
  f <- stepfun(s, y)
  return(f)
}

Q.hat.c <- function(obj) {
  s <- obj@s
  d <- obj@d
  N <- cumsum(d)
  new(StepFunction, s, c(0, N)/length(obj@y))
}

R.hat <- function(obj) {
  s <- obj@s
  y <- sapply(c(0, s), function(u) mean(sapply(1:length(obj@t), function(i) {
    sum(obj@t[[i]] <= u & u <= obj@y[i])
  })))
  f <- stepfun(s, y)
  return(f)
}

R.hat.c <- function(obj) {
  s <- obj@s
  d <- obj@d
  N <- cumsum(d)
  y.i <- order(obj@y)
  m <- sapply(obj@t, length)
  N <- N + eval_N(s, obj@y[y.i], m[y.i])
  new(StepFunction, s, c(0, N/length(obj@y)))	
}

b.hat <- function(obj, i) {
  R <- R.hat(obj)
  R.t <- R(obj@t[[i]])
  k.single <- function(u) {
    if (u > obj@y[i]) return(0)
    sum(obj@t[[i]] <= u) / R(u)^2		
  }
  k <- function(u) {
    sapply(u, k.single)
  }
  function(t) {
    step_integrate(k, Q.hat(obj), t, obj@T_0) - sum(as.numeric(t < obj@t[[i]]) / R.t)
  }
}

b.hat.gen <- function(obj) {
  R <- R.hat.c(obj)
  Q <- Q.hat.c(obj)
  return(function(i) {
    R.t <- R$sort_call(obj@t[[i]])
    x <- c(obj@t[[i]], obj@y[i])
    y <- append(0:length(obj@t[[i]]), 0)
    k.numerator <- new(StepFunction, x, y)
    k <- k.numerator / R^2
    function(t) {
      if (t == obj@T_0) return(0)
      step_integrate.StepFunction(k, Q, t, obj@T_0) - sum(as.numeric(t < obj@t[[i]]) / R.t)
    }
  })
}

b.hat.y <- function(obj) {
  R <- R.hat.c(obj)
  Q <- Q.hat.c(obj)
  retval <- matrix(0.0, obj@n, obj@n)
  for(i in seq_len(obj@n)) {
    R.t <- R$sort_call(obj@t[[i]])
    x <- c(obj@t[[i]], obj@y[i])
    y <- append(0:length(obj@t[[i]]), 0)
    k.numerator <- new(StepFunction, x, y)
    k <- k.numerator / R^2
    for(j in seq_len(obj@n)) {
      t <- obj@y[j]
      if (t == obj@T_0) next
      retval[i,j] <- step_integrate.StepFunction(k, Q, t, obj@T_0) - sum(as.numeric(t < obj@t[[i]]) / R.t)
    }
  }
  retval
}
