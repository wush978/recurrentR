Delta_i <- function(obj) {
	as.integer(obj@y == obj@T_0)
}

Z <- function(obj, F.hat = NULL, gamma = NULL) {
	if (is.null(F.hat)) F.hat <- obj$F.hat
	if (is.null(gamma)) gamma <- obj$U.hat()
	m <- sapply(obj@y, length)
	y.inv <- F.hat()
}


