Delta_i <- function(obj) {
	as.integer(obj@y == obj@T_0)
}

Z_i.hat <- function(obj, F.hat = NULL, gamma = NULL) {
	if (is.null(F.hat)) F.hat <- obj$F.hat
	if (is.null(gamma)) gamma <- obj$U.hat()
	m <- sapply(obj@y, length)
	y.inv <- 1 / F.hat(obj@y)
	return(m / y.inv * exp(obj@X %*% gamma))
}

Delta_i <- function(obj) {
	stopifnot(length(obj@D) == length(obj@y))
	as.integer(obj@D < obj@y)
}

#'@title Borrow-Strength Method
BSM <- function(obj) {
	Zi <- Z_i.hat(obj)
	D <- Delta_i(obj)
	a_i = D * as.integer(obj@y <= obj@T_0)
	n <- length(obj@y)
	y.index <- order(obj@y, decreasing=TRUE)
	X <- obj@X[y.index,]
	U <- function(Beta) {
		Si <- (Zi * exp(obj@X %*% Beta))[y.index]
		SSi <- cumsum(Si) # denominator of right part of eq. (2) in Huang 2004
		as.vector((D - sapply(1:n, function(i) {
			sum((Si[i] / SSi)[i:n])
		})) %*% X)
	}
	dU <- function(Beta) {
		Si <- (Zi * exp(obj@X %*% Beta))[y.index]
		SSi <- cumsum(Si) # denominator of right part of eq. (2) in Huang 2004
		DSSi <- .Call("DSSj", Si, X, PACKAGE="recurrentR")
	}
}