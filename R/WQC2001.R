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
	retval <- function(t) {
		sapply(t, function(t) retval.single(t))
	}
	return(retval)
}

Lambda.hat <- function(obj) {
	F.hat <- obj$F.hat
	m <- sapply(obj@t, length)
	Lambda.hat.T_0 <- mean(m / F.hat(obj@y))
	return(function(t) {
		Lambda.hat.T_0 * F.hat(t)
	})
}