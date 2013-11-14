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

t_ij_vs_y_k.gen <- function(obj) {
  if (!exists("t_ij_vs_y_k", envir=obj@cache, inherits=FALSE)) {
    obj@cache[["t_ij_vs_y_k"]] <- lapply(seq_len(obj@n), function(i) {
      temp <- outer(obj@t[[i]], obj@y, ">=")
      array(as.integer(temp), dim(temp))
    })
  }
  obj@cache[["t_ij_vs_y_k"]]
}

t_ij_vs_s.gen <- function(obj) {
  if (!exists("t_ij_vs_s", envir=obj@cache, inherits=FALSE)) {
    #**slow**
    obj@cache[["t_ij_vs_s"]] <- lapply(seq_len(obj@n), function(i) {
      temp <- outer(obj@s, obj@t[[i]], ">=")
      array(as.integer(temp), dim(temp))
    })
    #**slow**
  }
  obj@cache[["t_ij_vs_s"]]
}

s_vs_y_i.gen <- function(obj) {
  if (!exists("s_vs_y_i", envir=obj@cache, inherits=FALSE)) {
    s_vs_y_i <- outer(obj@s, obj@y, "<=")
    s_vs_y_i <- array(as.integer(s_vs_y_i), dim(s_vs_y_i))
    obj@cache[["s_vs_y_i"]] <- s_vs_y_i
  }
  obj@cache[["s_vs_y_i"]]
}

y_k_vs_y_i.gen <- function(obj) {
  if (!exists("y_k_vs_y_i", envir=obj@cache, inherits=FALSE)) {
    y_k_vs_y_i <- outer(obj@y, obj@y, "<=")
    obj@cache[["y_k_vs_y_i"]] <- array(as.integer(y_k_vs_y_i), dim(y_k_vs_y_i))
  }
  obj@cache[["y_k_vs_y_i"]]
}

y_k_vs_s.gen <- function(obj) {
  if (!exists("y_k_vs_s", envir=obj@cache, inherits=FALSE)) {
    y_k_vs_s <- outer(obj@y, obj@s, "<=")
    obj@cache[["y_k_vs_s"]] <- array(as.integer(y_k_vs_s), dim(y_k_vs_s))
  }
  obj@cache[["y_k_vs_s"]]
}

b.hat.y.gen <- function(obj) {
  if (!exists("b.hat.y", envir=obj@cache, inherits=FALSE)) {
    s <- obj@s
    d <- obj@d
    t_ij_vs_y_k <- t_ij_vs_y_k.gen(obj)
    t_ij_vs_s <- t_ij_vs_s.gen(obj)
    s_vs_y_i <- s_vs_y_i.gen(obj)
    y_k_vs_y_i <- y_k_vs_y_i.gen(obj)
    y_k_vs_s <- y_k_vs_s.gen(obj)
    R <- R.hat.c(obj)
    R__t_ij <- lapply(seq_len(obj@n), function(i) R$sort_call(obj@t[[i]]))
    R__s <- R$sort_call( s )
    retval <- matrix(numeric(obj@n^2), obj@n, obj@n)
    for(i in seq_len(obj@n)) {
      if (length(obj@t[[i]]) == 0) next
      term_2.i <- - (1/R__t_ij[[i]]) %*% t_ij_vs_y_k[[i]]
      term_1.i <- apply(
  #       #*** bottleneck ***
        y_k_vs_s %*% (t_ij_vs_s[[i]] * (s_vs_y_i[,i] * d  / (obj@n * R__s^2)))
  #       #*** bottleneck ***
        , 1, sum)
      retval[i,] <- term_1.i + term_2.i
  #     ***prototype***
  #     for(k in seq_len(obj@n)) {
  #       term_2 <- -sum(t_ij_vs_y_k[[i]][,k] / R__t_ij[[i]])
  #       term_1 <- sum((s_vs_y_i[,i] * d * y_k_vs_s[k,] / R__s^2) %*% t_ij_vs_s[[i]]) / obj@n
  #       retval[i,k] <- term_1 + term_2
  #     }
    }
    obj@cache[["b.hat.y"]] <- retval
  }
  obj@cache[["b.hat.y"]]
}
