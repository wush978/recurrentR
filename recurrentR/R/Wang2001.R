R.s.gen <- function(obj) {
  key <- "R.s"
  if (!is_cache(obj, key)) {
    s <- obj@s
    d <- obj@d
    N <- cumsum(d)
    y.i <- order(obj@y)
    m <- sapply(obj@t, length)
    N <- N + eval_N(s, obj@y[y.i], m[y.i])
    obj@cache[[key]] <- N
  }
  obj@cache[[key]]
}

Lambda_0.hat.gen <- function(obj) {
  if (!exists("Lambda_0.hat", envir=obj@cache, inherits=FALSE)) {
    N <- R.s.gen(obj)
    x <- obj@s
    y <- append(rev(cumprod(rev(1 - obj@d/N))), 1)
    f <- stepfun(x, y)
    obj@cache[["Lambda_0.hat"]] <- f
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

m.gen <- function(obj) {
  if (!exists("m", envir=obj@cache, inherits=FALSE)) {
    obj@cache[["m"]] <- sapply(obj@t, length)
  }
  obj@cache[["m"]]
}

w.gen <- function(obj) {
  if (!exists("w", envir=obj@cache, inherits=FALSE)) {
    obj@cache[["w"]] <- rep(1, obj@n)
  }
  obj@cache[["w"]]
}

gamma.bar.hat.gen <- function(obj) {
  if (!exists("gamma.bar.hat", envir=obj@cache, inherits=FALSE)) {
    b <- m.gen(obj) * Lambda_0.hat.y.inv.gen(obj)
    W.bar <- cbind(1, obj@W)
    gamma.bar.hat <- rep(0, ncol(W.bar))
    g <- function(gamma) t(W.bar) %*% (b - exp(W.bar %*% gamma))
    g.grad <- function(gamma) - t(W.bar) %*% diag(c(exp(W.bar %*% gamma)), nrow(W.bar), nrow(W.bar)) %*% W.bar
    slv <- nleqslv(gamma.bar.hat, g, jac=g.grad)
    if (sum(abs(g(slv$x))) > obj@tol) stop("Failed to converge during solving gamma")
    obj@cache[["gamma.bar.hat"]] <- slv$x
    if (!is.null(colnames(obj@W))) names(obj@cache[["gamma.bar.hat"]]) <- c("(Intercept)", colnames(obj@W))
  }
  obj@cache[["gamma.bar.hat"]]
}

gamma.hat.gen <- function(obj) {
  if (!is_cache(obj, "gamma.hat")) {
    obj@cache[["gamma.hat"]] <- gamma.bar.hat.gen(obj)[-1]
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
  key <- "Q.hat.c"
  if (!is_cache(obj, key)) {
    N <- cumsum(obj@d)
    obj@cache[[key]] <- new(StepFunction, obj@s, c(0, N)/obj@n)
  }
  obj@cache[[key]]
}

R.hat <- function(obj) {
  function(u) {
    retval <- 0
    for(i in seq_along(obj@t)) {
      if (u > obj@y[i]) next
      for(j in seq_along(obj@t[[i]])) {
        if (u < obj@t[[i]][j]) next
        retval <- retval + 1
      }
    }
    retval / obj@n
  }
}

# This implemtation is incorrect. However, it gives the correct answer on `obj@s` which is the place only needed.
R.hat.c <- function(obj) {
  key <- "R.hat.c"
  if (!is_cache(obj, key)) {
    N <- R.s.gen(obj)
    obj@cache[[key]] <- new(StepFunction, obj@s, c(0, N/length(obj@y)))	
  }
  obj@cache[[key]]
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

# c.hat.gen <- function(obj) {
#   if (!exists("c.hat", envir=obj@cache, inherits=FALSE)) {
#     m <- m.gen(obj)
#     Lambda.hat.T_0 <- Lambda_0.hat.gen(obj)(obj@T_0)
#     Lambda_0.hat.y.inv <- Lambda_0.hat.y.inv.gen(obj)
#     b.hat.y <- b.hat.y.gen(obj)
#     obj@cache[["c.hat"]] <- lapply(seq_len(obj@n), function(i) {
#       force(i)
#       bi.y <- b.hat.y[i,]
#       mean(bi.y * m * Lambda_0.hat.y.inv) + m[i] * Lambda_0.hat.y.inv[i] - Lambda.hat.T_0
#     })
#   }
#   obj@cache[["c.hat"]]
# }

c.hat.i.gen <- function(obj) {
  if (!exists("c.hat.i", envir=obj@cache, inherits=FALSE)) {
    m <- m.gen(obj)
    Lambda.hat.T_0 <- Lambda_0.hat.gen(obj)(obj@T_0)
    Lambda_0.hat.y.inv <- Lambda_0.hat.y.inv.gen(obj)
    b.hat.y <- b.hat.y.gen(obj)
    retval <- numeric(obj@n)
    for(i in seq_len(obj@n)) {
      bi.y <- b.hat.y[i,]
      retval[i] <- mean(bi.y * m * Lambda_0.hat.y.inv) + m[i] * Lambda_0.hat.y.inv[i] - Lambda.hat.T_0
    }
    obj@cache[["c.hat.i"]] <- retval
  }
  obj@cache[["c.hat.i"]]
}

e.hat.i.gen <- function(obj) {
  if (!exists("e.hat.i", envir=obj@cache, inherits=FALSE)) {
    m <- m.gen(obj)
    Lambda_0.hat.y.inv <- Lambda_0.hat.y.inv.gen(obj)
    gamma.bar.hat <- gamma.bar.hat.gen(obj)
    w <- rep(1, obj@n)
    b.hat.y <- b.hat.y.gen(obj)
    W.bar <- cbind(1, obj@W)
    retval <- vector("list", obj@n)
    for(i in seq_len(obj@n)) {
      term_1 <- (w * m * b.hat.y[i,] * Lambda_0.hat.y.inv)
      term_2 <- m[i] *  Lambda_0.hat.y.inv[i]
      retval[[i]] <- - as.vector(term_1 %*% W.bar) / obj@n + (w[i] * W.bar[i,] * (term_2  - exp(W.bar[i,] %*% gamma.bar.hat)))
    }
    obj@cache[["e.hat.i"]] <- do.call(cbind, retval)
  }
  obj@cache[["e.hat.i"]]
}

Psi.bar.hat.gen <- function(obj) {
  if (!is_cache(obj, "Psi.bar.hat")) {
    w <- w.gen(obj)
    gamma.bar.hat <- gamma.bar.hat.gen(obj)
    W.bar <- cbind(1, obj@W)
    dei.dgamma <- function(i) {
      as.vector(- w[i] * exp(W.bar[i,] %*% gamma.bar.hat)) * (W.bar[i,] %*% t(W.bar[i,]))
    }
    dei.dgamma.i <- lapply(seq_len(obj@n), function(i) -dei.dgamma(i))
    obj@cache[["Psi.bar.hat"]] <- Reduce("+", dei.dgamma.i) / obj@n
  }
  obj@cache[["Psi.bar.hat"]]
}

Psi.bar.hat.inv.gen <- function(obj) {
  if (!is_cache(obj, "Psi.bar.hat.inv")) {
    Psi.bar.hat <- Psi.bar.hat.gen(obj)
    obj@cache[["Psi.bar.hat.inv"]] <- solve(Psi.bar.hat)
  }
  obj@cache[["Psi.bar.hat.inv"]]
}
    

fi.hat.i.gen <- function(obj) {
  if (!is_cache(obj, "fi.hat.i")) {
    w <- w.gen(obj)
    gamma.bar.hat <- gamma.bar.hat.gen(obj)
    Psi.bar.hat.inv <- Psi.bar.hat.inv.gen(obj)
    ei.seq <- e.hat.i.gen(obj)
#    TO CHECK
    obj@cache[["fi.hat.i"]] <- Psi.bar.hat.inv %*% ei.seq
  }
  obj@cache[["fi.hat.i"]]
}

#'@title Semi-parametric Estimator of the Cumulative Rate Function
#'@param obj A \code{recurrent-data} object.
#'@param methods One of \code{c("none", "bootstrap", "asymptotic")}. The method of evaluating standard deviation.
#'@param B An \code{integer} value. The size of bootstrap.
#'@return A \code{list} of: \enumerate{
#'\item \code{Lambda_0.hat}: Function. The estimated cumulative rate function. Note that 
#'the function is normalized such that \code{\Lambda_0.hat(T_0)} is 1.
#'\item \code{gamma.bar.hat}: Numeric Vector. The estimated regression parameter \eqn{\gamma}.
#'\item \code{Lambda_0.hat.var}: Function. The variance of \code{Lambda_0.hat}.
#'\item \code{gamma.bar.hat.var}: Numeric Matrix. The estimated variance of \code{gamma.bar.hat}.
#'}
#'@details This is an implementation of non-parametric estimator of the cumulative rate
#'function based on [Wang et. al. 2001]. The recurrent event time, \eqn{t_{i,1}}
#', \eqn{t_{i,2}}, ..., \eqn{t_{i,m_i}}, are the realization of a poisson process 
#'\eqn{N_i(.)} whose itensity is moded as \deqn{\lambda_i(t) = \lambda_0(t) z_i exp(W_i \gamma)},
#'where:
#'\enumerate{
#'\item \eqn{z_i} is a nonnegative-valued latent variable such that \eqn{E(z_i | W_i) = E(z_i)}.
#'\item The baseline intensity function \eqn{\lambda_0(t)} is a probability function: \enumerate{
#'  \item \eqn{\lambda_0(t) \neq 0}
#'  \item \eqn{\Lambda_0(T_0) = \int_0^{T_0} \lambda_0(u) du = 1}
#'  \item \eqn{\gamma} is a \eqn{R^{1 \times q}} vector.
#'  \item Condition on \eqn{W_i, z_i}, \eqn{N_i(.)} and \eqn{y_i} are independent.
#'  }
#'}
#'@author Wush Wu
#'@references M-C., Wang, J. Qin, and C.-T. Chiang. 2001. “Analyzing 
#'Recurrent Event Data With Informative Censoring.” 
#'Journal of the American Statistical Association 96: 1057–1065. 
#'http://EconPapers.repec.org/RePEc:bes:jnlasa:v:96:y:2001:m:september:p:1057-1065.
#'@examples 
#'\dontrun{
#'library(survrec)
#'data(MMC)
#'obj <- create_recurrent_data.data.frame(
#'  MMC, id = "id", time = "time", time_type = "relatively",
#'  indicator = "event", indicator_value = list("recurrent" = 1, "censor" = 0),
#'  covariate = "group"
#')
#'wang_2001 <- Wang2001(obj)
#'# Plot estimated Lambda_0.hat
#'Lambda_0.hat <- wang_2001$Lambda_0.hat
#'curve(Lambda_0.hat, 0, obj@@T_0)
#'}
#'@export
Wang2001 <- function(obj, methods = c("none", "bootstrap", "asymptotic"), B = 100) {
  Lambda_0.hat <- Lambda_0.hat.gen(obj)
  gamma.bar.hat <- gamma.bar.hat.gen(obj)
  if (methods[1] == "none") {
    return(list(
      Lambda_0.hat = Vectorize(Lambda_0.hat),
      gamma.bar.hat = gamma.bar.hat
      ))
  }
  if (methods[1] == "bootstrap") {
    gamma.bar.hat.Bootstrap <- Lambda_0.hat.Bootstrap <- vector("list", B)
    for(i in seq_len(B)) {
      index.resampled <- sample(seq_len(obj@n), obj@n, TRUE)
      obj.resampled <- create_recurrent_data.numeric(
        y=obj@y[index.resampled],
        D=obj@D[index.resampled],
        t=obj@t[index.resampled],
        T_0=obj@T_0,
        W=obj@W[index.resampled, , drop = FALSE],
        tol=obj@tol
        )
      temp <- Wang2001(obj.resampled, "none")
      Lambda_0.hat.Bootstrap[[i]] <- temp$Lambda_0.hat
      gamma.bar.hat.Bootstrap[[i]] <- temp$gamma.bar.hat
    }
    return(list(
      Lambda_0.hat = Vectorize(Lambda_0.hat),
      gamma.bar.hat = gamma.bar.hat,
      Lambda_0.hat.var = Vectorize(function(t) {
        var(sapply(Lambda_0.hat.Bootstrap, function(f) f(t)))
      }),
      gamma.bar.hat.var = var(do.call(rbind, gamma.bar.hat.Bootstrap))
      ))
  }
  if (methods[1] == "asymptotic") {
    b.hat <- b.hat.gen(obj)
    b <- lapply(seq_len(obj@n), b.hat)
    psi.inv <- Psi.bar.hat.inv.gen(obj)
    ei.seq <- e.hat.i.gen(obj)
#     fi.seq <- fi.hat.i.gen(obj)
    return(list(
      Lambda_0.hat = Vectorize(Lambda_0.hat),
      gamma.bar.hat = gamma.bar.hat,
      Lambda_0.hat.var = Vectorize(function(t) {
        Lambda_0.hat(t)^2 * mean(sapply(b, function(f) f(t)^2)) / obj@n
      }),
      gamma.bar.hat.var = psi.inv %*% var(t(ei.seq)) %*% psi.inv / obj@n
    ))
  }
}
