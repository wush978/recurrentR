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
	retval <- function(t, bootstrap.se = FALSE, B = 100, error.measurement.function = sd) {
		if (!bootstrap.se) return(sapply(t, function(t) retval.single(t)))
		estimate <- sapply(t, function(t) retval.single(t))
		obj.b <- obj$bootstrap(B)
		F.hat.bootstrap.result <- sapply(obj.b, function(obj) {
			obj$F.hat(t)
		})
		error.measurement = apply(F.hat.bootstrap.result, 1, error.measurement.function)
		list(estimate=estimate, error.measurement=error.measurement)
	}
	return(retval)
}

Lambda.hat <- function(obj, bootstrap.se = FALSE) {
	F.hat <- obj$F.hat
	m <- sapply(obj@t, length)
	Lambda.hat.T_0 <- mean(m / F.hat(obj@y))
	return(function(t) {
		Lambda.hat.T_0 * F.hat(t)
	})
}