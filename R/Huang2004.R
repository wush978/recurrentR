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
BSM <- function(obj, tol, verbose = TRUE) {
	Zi <- Z_i.hat(obj)
	D <- Delta_i(obj)
	a_i = D * as.integer(obj@y <= obj@T_0)
	n <- length(obj@y)
	y.index <- order(obj@y, decreasing=TRUE)
	X <- obj@X[y.index, -1]
	Si <- function(Beta) as.vector(Zi * exp(X %*% Beta))
	DSi <- function(Beta) t(X) %*% diag(Si(Beta))
	SSi <- function(Beta) cumsum(Si(Beta))
	DSSi <- function(Beta) t(apply(DSi(Beta), 1, cumsum))
	U <- function(Beta) {
		S <- Si(Beta)
		SS <- SSi(Beta)
		retval <- Beta
		retval[] <- 0
		for(i in 1:n) {
			retval <- retval + (D[i] - S[i] * sum((D / SS)[i:n])) * X[i,]
		}
		retval
	}
	
	DU <- function(Beta) {
		S <- Si(Beta)
		DS <- DSi(Beta)
		SS <- SSi(Beta)
		DSS <- DSSi(Beta)
		retval <- matrix(0.0, ncol=length(Beta), nrow=length(Beta))
		for(i in 1:n) {
			temp <- Beta
			temp[] <- 0
			for(j in i:n) {
				temp <- temp - D[j] * (DS[,i] * SS[j] - S[i] * DSS[,j]) / SS[j]^2
			}
			retval <- retval + matrix(temp, ncol = 1) %*% X[i,]
		}
		retval
	}
	
	temp <- nleqslv(x, fn=U, jac=DU)
	if(verbose) {
		cat(sprintf("Check if gamma is solved correctly: %s \n", paste(g(temp$x), collapse=",")))
		cat(sprintf("message of nleqslv: %s ", temp$message))
	}
	if (sum(abs(g(temp$x))) > tol) stop("Failed to converge during solving gamma")
	return(temp$x)
}