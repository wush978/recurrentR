s <- function(obj) {
	sort(unique(unlist(sapply(obj@t, unique))))
}

d.single <- function(obj, s) {
	sum(sapply(obj@t, function(t) {
		sum(t == s)
	}))
}

d <- function(obj, s) {
	sapply(s, function(s) d.single(obj, s))
}

N.single <- function(obj, s) {
	sum((obj@y >= s) * sapply(obj@t, function(t) {
		sum(t <= s)
	}))
}

N <- function(obj, s) {
	sapply(s, function(s) N.single(obj, s))
}

F.hat <- function(obj) {
	s <- s(obj)
	d <- d(obj, s)
	N <- N(obj, s)
	retval.single <- function(t) {
		prod(ifelse(s > t, 1 - d/N, 1))
	}
	retval <- function(t, bootstrap = FALSE, B = 100, error.measurement.function = base::sd) {
		if (!bootstrap) return(sapply(t, function(t) retval.single(t)))
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
	Lambda.hat.T_0 <- mean(m / F.hat(obj@y))
	return(function(t, bootstrap = FALSE, B = 100, error.measurement.function = base::sd) {
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
U.hat.solve <- function(X, b, tol) {
	gamma <- solve(t(X) %*% X, t(X) %*% (b-1))
	f1 <- function(gamma) {
		t(X) %*% (b - exp(X %*% gamma))
	}
	f2 <- function(gamma) {
		t(X) %*% diag(c(exp(X %*% gamma)), nrow(X), nrow(X)) %*% X
	}
# 	f2.gamma <- f2(gamma)
# 	gamma.new <- solve(f2.gamma, f2.gamma %*% gamma + f1(gamma))
# 	while (sum(abs(f1(gamma.new))) > tol) {
# 		gamma <- gamma.new
# 		f2.gamma <- f2(gamma)
# 		gamma.new <- solve(f2.gamma, f2.gamma %*% gamma + f1(gamma), tol=tol)
# 	}
	temp <- nleqslv(gamma, f1, jac=f2)
	if (length(temp$msg) > 0 && nchar(temp$msg) > 0) warning(temp$msg)
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
	function(bootstrap = FALSE, B = 100, error.measurement.function = sd, tol=1e-4) {
		m <- sapply(obj@t, length)
		y <- obj@y
		F.hat <- obj$F.hat
		b <- m / obj$F.hat(y)
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