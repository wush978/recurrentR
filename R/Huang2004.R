Z_i.hat <- function(obj, F.hat = NULL, gamma = NULL) {
	if (is.null(F.hat)) F.hat <- obj$F.hat
	if (is.null(gamma)) gamma <- obj$U.hat()
	m <- sapply(obj@y, length)
	y.inv <- 1 / F.hat(obj@y)
	return(m / y.inv * exp(obj@X %*% gamma))
}

#'@title Indicator of censoring type.
#'@description
#'\eqn{D_i = 1} if the failure time is observed.
Delta_i <- function(obj) {
	stopifnot(length(obj@D) == length(obj@y))
	as.integer(obj@D < obj@y)
}

#'@title Moment Estimators associated with subject-specific latent variable
#'@description
#'\eqn{S_i = Z_i e^{X_i \beta}}
#'The \code{Zi} and \code{X} are subject-specific latent varialbe and time-independent covariate
#'associated with ordered \code{y}
#'@return numeric vector
Si <- function(Beta, Zi, X) as.vector(Zi * exp(X %*% Beta))

#'@title Moment Estimators associated with subject-specific latent variable
#'@description
#'\eqn{\frac{d}{d\beta} S_i = S_i X_i^T}.
#'The \code{Zi} and \code{X} are subject-specific latent varialbe and time-independent covariate
#'associated with ordered \code{y}
#'@seealso \code{\link{Si}}
#'@return numeric vector
DSi <- function(Beta, Zi, X) t(X) %*% diag(Si(Beta, Zi, X))

#'@title Moment Estimators associated with subject-specific latent variable
#'@description
#'\eqn{SS_i = \sum_{j=1}^{i} S_j}.
#'The \code{Zi} and \code{X} are subject-specific latent varialbe and time-independent covariate
#'associated with ordered \code{y}
#'@seealso \code{\link{Si}}
#'@return numeric vector
SSi <- function(Beta, Zi, X) cumsum(Si(Beta, Zi, X))

#'@title Moment Estimators associated with subject-specific latent variable
#'@description
#'\eqn{\frac{d}{d\beta}SS_i = \sum_{j=1}^{i} S_j X_j^T}
#'The \code{Zi} and \code{X} are subject-specific latent varialbe and time-independent covariate
#'associated with ordered \code{y}
#'@seealso \code{\link{SSi}}
#'@return numeric vector
DSSi <- function(Beta, Zi, X) t(apply(DSi(Beta, Zi, X), 1, cumsum))


#'@title Borrow-Strength Method
#'
#'@description
#'Estimate regression coefficient of the hazard function of the failure time.
BSM <- function(obj, tol = 1e-7, verbose = FALSE) {
	y.index <- order(obj@y, decreasing=TRUE)
	y <- obj@y[y.index]
	Zi <- Z_i.hat(obj)[y.index]
	D <- Delta_i(obj)[y.index]
	n <- length(obj@y)
	stopifnot(ncol(obj@X) > 1) # TODO
	X <- obj@X[y.index, -1]
	y.rank <- rank(y, ties.method="min")
	browser()
	U <- function(Beta) {
		S <- Si(Beta)
		SS <- SSi(Beta)
		retval <- Beta
		retval[] <- 0
		for(i in 1:n) {
			temp <- D[i] * X[i,]
			for(j in 1:n) {
				if (y[j] < y[i]) next
				temp <- temp - D[i] * S[j]*X[j,] / SS[i]
			}
			retval <- retval + temp
		}
		retval
	}
	
	DU <- function(Beta) {
		S <- Si(Beta, Zi, X)
		DS <- DSi(Beta, Zi, X)
		SS <- SSi(Beta, Zi, X)
		DSS <- DSSi(Beta, Zi, X)
		retval <- matrix(0.0, ncol=length(Beta), nrow=length(Beta))
		for(i in 1:n) {
			temp <- retval
			temp[] <- 0
			for(j in 1:n) {
				if (y[j] < y[i]) next
				temp <- temp - D[i] * matrix((DS[,j] * SS[i] - S[j] * DSS[,i]) / SS[i]^2, ncol=1) %*% X[j,]
			}
			retval <- retval + temp
		}
		retval
	}
	
	temp <- nleqslv(rep(0, ncol(obj@X) - 1), fn=U, jac=DU)
	if(verbose) {
		cat(sprintf("Check if gamma is solved correctly: %s \n", paste(U(temp$x), collapse=",")))
		cat(sprintf("message of nleqslv: %s ", temp$message))
	}
	if (sum(abs(U(temp$x))) > tol) stop("Failed to converge during solving gamma")
	return(temp$x)
}

inverse <- function(x) x[length(x):1]

H0.hat <- function(obj, ...) {
	Beta <- BSM(obj, ...)
	y.index <- order(obj@y, decreasing=TRUE)
	y <- obj@y[y.index]
	Zi <- Z_i.hat(obj)[y.index]
	D <- Delta_i(obj)[y.index]
	n <- length(obj@y)
	stopifnot(ncol(obj@X) > 1) # TODO
	X <- obj@X[y.index, -1]
	S <- Si(Beta, Zi, X)
	SS <- SSi(Beta, Zi, X)
	browser()	
}
