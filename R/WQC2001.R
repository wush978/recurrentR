s <- function(obj) {
	sort(unique(unlist(sapply(obj@t, unique))))
}

# d.single <- function(obj, s) {
# 	sum(sapply(obj@t, function(t) {
# 		sum(t == s)
# 	}))
# }

d <- function(obj) {
# 	sapply(s, function(s) d.single(obj, s))
	as.vector(table(unlist(obj@t)))
}

# N.single <- function(obj, s) {
# 	sum((obj@y >= s) * sapply(obj@t, function(t) {
# 		sum(t <= s)
# 	}))
# }

# N <- function(obj) {
# 	cumsum(d(obj))
# }

F.hat <- function(obj) {
	s <- s(obj)
	d <- d(obj)
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

Lambda.hat <- function(obj, bootstrap = FALSE) {
	F.hat <- obj$F.hat
	m <- sapply(obj@t, length)
	y.inv <- 1 / F.hat(obj@y)
	if (any(is.nan(y.inv))) stop("TODO: ill condition")
	if (ncol(obj@X) == 0) {
		Lambda.hat.T_0 <- mean(m * y.inv)
	} else {
		gamma <- obj$U.hat()
		Lambda.hat.T_0 <- exp(gamma[1])
	}
	return(function(t, bootstrap = FALSE, B = 100, error.measurement.function = stats::sd) {
		if (!bootstrap) return(Lambda.hat.T_0 * F.hat(t))
		estimate <- obj$Lambda.hat(t)
		obj.b <- obj$bootstrap(B)
		Lambda.hat.bootstrap.result <- sapply(obj.b, function(obj) {
			obj$Lambda.hat(t)
		})
		error.measurement = apply(Lambda.hat.bootstrap.result, 1, error.measurement.function)
		list(estimate=estimate, error.measurement=error.measurement)
	})
}

#'@title Solve 
U.hat.solve <- function(X, b, tol, verbose = FALSE, gamma = NULL) {
# 	gamma <- solve(t(X) %*% X, t(X) %*% (b-1))
	if (is.null(gamma)) gamma <- rep(0, ncol(X))
	g <- function(gamma) {
		t(X) %*% (b - exp(X %*% gamma))
	}
	g.gradient <- function(gamma) {
		- t(X) %*% diag(c(exp(X %*% gamma)), nrow(X), nrow(X)) %*% X
	}
	temp <- nleqslv(gamma, g, jac=g.gradient)
	if(verbose) {
		cat(sprintf("Check if gamma is solved correctly: %s \n", paste(g(temp$x), collapse=",")))
		cat(sprintf("message of nleqslv: %s ", temp$message))
	}
	if (sum(abs(g(temp$x))) > tol) stop("Failed to converge during solving gamma")
	return(temp$x)
}

sapply_pb <- function(X, FUN, ...)
{
	env <- environment()
	pb_Total <- length(X)
	counter <- 0
	pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
	
	wrapper <- function(...){
		curVal <- get("counter", envir = env)
		assign("counter", curVal +1 ,envir=env)
		setTxtProgressBar(get("pb", envir=env), curVal +1)
		FUN(...)
	}
	res <- sapply(X, wrapper, ...)
	close(pb)
	res
}


U.hat <- function(obj) {
	function(bootstrap = FALSE, B = 100, error.measurement.function = stats::sd, tol=1e-4) {
		m <- sapply(obj@t, length)
		y <- obj@y
		F.hat <- obj$F.hat
		F.hat.y <- F.hat(y)
    if (sum(F.hat.y == 0) > 0) {
      stopifnot(m[F.hat.y == 0] == 0)
    }
		b <- m / F.hat.y
    b[is.nan(b)] <- 0
		estimate <- U.hat.solve(obj@X, b, tol)
		if (!bootstrap) {
			return(estimate)
		}
		obj.b <- obj$bootstrap(B)
		U.hat.bootstrap.result <- sapply_pb(obj.b, function(obj) {
			obj$U.hat(B = B)
		})
		error.measurement = apply(U.hat.bootstrap.result, 1, error.measurement.function)
		list(estimate=estimate, error.measurement=error.measurement)
	}
}

Q.hat <- function(obj) {
	s <- s(obj)
	y <- sapply(c(0, s), function(u) mean(sapply(obj@t, function(t) sum(t <= u))))
	f <- stepfun(s, y)
	return(f)
}

Q.hat.c <- function(obj) {
	s <- s(obj)
	d <- d(obj)
	N <- cumsum(d)
	new(StepFunction, s, c(0, N)/length(obj@y))
}

R.hat <- function(obj) {
	s <- s(obj)
	y <- sapply(c(0, s), function(u) mean(sapply(1:length(obj@t), function(i) {
		sum(obj@t[[i]] <= u & u <= obj@y[i])
	})))
	f <- stepfun(s, y)
	return(f)
}

R.hat.c <- function(obj) {
	s <- s(obj)
	d <- d(obj)
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


c.hat.gen <- function(obj, F.hat = NULL, bi.gen = NULL) {
	m <- sapply(obj@t, length)
	if (is.null(bi.gen)) bi.gen <- b.hat.gen(obj)
	if (is.null(F.hat)) F.hat <- obj$F.hat
	Lambda.hat.T_0 <- obj$Lambda.hat(obj@T_0)
	F.hat.y <- F.hat(obj@y)
	return(function(i) {
		bi <- bi.gen(i)
		bi.y <- sapply(obj@y, function(y) bi(y))
		mean(bi.y * m / F.hat.y) + m[i] / F.hat.y[i] - Lambda.hat.T_0
	})	
}

d.hat.gen <- function(obj, F.hat = NULL, Lambda.hat = NULL, bi.gen = NULL, ci.gen = NULL) {
	if (is.null(F.hat)) F.hat <- obj$F.hat
	if (is.null(Lambda.hat)) Lambda.hat <- obj$Lambda.hat
	if (is.null(bi.gen)) bi.gen <- b.hat.gen(obj)
	if (is.null(ci.gen)) ci.gen <- c.hat.gen(obj, F.hat = F.hat, bi.gen = bi.gen)
	return(function(i) {
		bi <- bi.gen(i)
		ci <- ci.gen(i)
		return(function(t) {
			F.hat(t) * (ci + Lambda.hat(obj@T_0) * bi(t))
		})
	})
}

e.hat.gen <- function(obj, F.hat = NULL, bi.gen = NULL, w = NULL, gamma = NULL) {
	if (is.null(F.hat)) F.hat <- obj$F.hat
	if (is.null(bi.gen)) bi.gen <- b.hat.gen(obj)
	if (is.null(w)) w <- rep(1, length(obj@y)) else stopifnot(length(w) == length(obj@y))
	if (is.null(gamma)) gamma <- obj$U.hat() else stopifnot(length(gamma) == ncol(obj@X))
	m <- sapply(obj@t, length)
	F.hat.y <- F.hat(obj@y)
	return(function(i) {
		bi <- bi.gen(i)
		term_1 <- (w * m * sapply(obj@y, bi) / F.hat.y)
    term_1[F.hat.y == 0] <- 0
		term_2 <- ifelse(F.hat.y[i] != 0, m[i] /  F.hat.y[i], 0)
		as.vector(term_1 %*% obj@X) / length(obj@y) + (w[i] * obj@X[i,] * (term_2  - exp(obj@X[i,] %*% gamma)))
	})	
}

fi.hat <- function(obj, w = NULL, gamma = NULL, psi.inv = NULL, ei.seq = NULL) {
  if(is.null(w)) w <- rep(1, length(obj@y)) else stopifnot(length(w) == length(obj@y))
  if (is.null(gamma)) gamma <- obj$U.hat() else stopifnot(length(gamma) == ncol(obj@X))
  if (is.null(psi.inv)) {
    dei.dgamma <- dei.dgamma.gen(obj, w, gamma)
    dei.dgamma.i <- lapply(1:length(obj@y), function(i) -dei.dgamma(i))
    psi <- Reduce("+", dei.dgamma.i) / length(obj@y)
    psi.inv <- solve(psi)
  }
  if (is.null(ei.seq)) {
    ei <- e.hat.gen(obj, gamma = gamma)
    ei.seq <- sapply(seq_along(obj@y), ei)
  }
  (psi.inv %*% ei.seq)
}
  
dei.dgamma.gen <- function(obj, w = NULL, gamma = NULL) {
	if(is.null(w)) w <- rep(1, length(obj@y)) else stopifnot(length(w) == length(obj@y))
	if (is.null(gamma)) gamma <- obj$U.hat() else stopifnot(length(gamma) == ncol(obj@X))
	function(i) {
		as.vector(- w[i] * exp(obj@X[i,] %*% gamma)) * (obj@X[i,] %*% t(obj@X[i,]))
	}
}

#'@title Asymptotic Variance Estimator
#'
#'@return list. 
#'
#'@export
asymptotic.var <- function(obj, w = NULL, gamma = NULL) {
	n <- length(obj@y)
	if (is.null(w)) w <- rep(1, n) else stopifnot(length(w) == n)
	if (is.null(gamma)) gamma <- obj$U.hat() else stopifnot(length(gamma) == ncol(obj@X))
	F.hat <- obj$F.hat
	b.hat.gen <- recurrentR:::b.hat.gen(obj)
	c.hat.gen <- recurrentR:::c.hat.gen(obj, F.hat=F.hat, b.hat.gen)
	d.hat.gen <- recurrentR:::d.hat.gen(obj, F.hat=F.hat, bi.gen=b.hat.gen, ci.gen=c.hat.gen)
	d <- list()
	for(i in seq_along(obj@y)) {
		d[[i]] <- d.hat.gen(i)
	}
	if (ncol(obj@X) == 0 | nrow(obj@X) == 0) return(list(Lambda.hat.var=function(t) {
		mean(sapply(1:n, function(i) d[[i]](t))^2 / n)
	}))
  dei.dgamma <- dei.dgamma.gen(obj, w, gamma)
  dei.dgamma.i <- lapply(1:length(obj@y), function(i) -dei.dgamma(i))
  psi <- Reduce("+", dei.dgamma.i) / n
  psi.inv <- solve(psi)
  ei <- e.hat.gen(obj)
  ei.seq <- sapply(seq_along(obj@y), ei)
	gamma.var.hat <- psi.inv %*% var(t(ei.seq)) %*% psi.inv / n
	fi.seq <- fi.hat(obj, w, gamma, psi.inv, ei.seq)[1,]
	b <- list()
	for(i in seq_along(obj@y)) {
		b[[i]] <- b.hat.gen(i)
	}
	colnames(gamma.var.hat) <- rownames(gamma.var.hat)
	return(list(Lambda.hat.var = function(t) {
		(F.hat(t)^2 * exp(2 * gamma[1]) * mean((sapply(b, function(b) b(t)) + fi.seq)^2)) / n
	}, gamma.var.hat = gamma.var.hat))
}

