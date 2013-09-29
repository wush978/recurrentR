#'@title Estimator of random effect \eqn{Z_i}
#'
#'@description \deqn{\hat{Z}_i = \frac{m_i}{\hat{\Lambda}_0(Y_i) e^{X_i \hat{\alpha}}}}
#'@param obj recurrent-data object
#'@param F.hat cache of \code{obj$F.hat}
#'@param F.hat cache of \code{obj$U.hat()}.
Z_i.hat <- function(obj, F.hat = NULL, gamma = NULL) {
	if (is.null(F.hat)) F.hat <- obj$F.hat
	if (is.null(gamma)) gamma <- obj$U.hat()
	m <- sapply(obj@t, length)
	F.hat.y <- F.hat(obj@y)
	if (sum(F.hat.y == 0) > 0) stop("TODO: ill condition")
	return(m / (F.hat.y * exp(obj@X %*% gamma)))
}

#'@title Indicator of censoring type.
#'@description
#'\eqn{D_i = 1} if the failure time is observed.
Delta_i <- function(obj) {
	as.integer(obj@D)
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

U.gen <- function(Zi, y, X, D, n) {
  function(Beta) {
    S <- Si(Beta, Zi, X)
    SS <- SSi(Beta, Zi, X)
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
    retval / n
  }
}

Gamma.hat.gen <- function(Zi, y, X, D, n) {
  function(Beta) {
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
    retval / n
  }  
}


#'@title Borrow-Strength Method
#'
#'@description
#'Estimate regression coefficient of the hazard function of the failure time.
#'
#'@param obj recurrent-data object
BSM <- function(obj, tol = 1e-7, verbose = FALSE) {
	y.index <- order(obj@y, decreasing=TRUE)
	y <- obj@y[y.index]
	Zi <- Z_i.hat(obj)[y.index]
	D <- Delta_i(obj)[y.index]
	n <- length(obj@y)
	stopifnot(ncol(obj@X) > 1) # TODO
	X <- obj@X[y.index, -1]
	y.rank <- rank(y, ties.method="min")
	U <- U.gen(Zi, y, X, D, n)
	DU <- Gamma.hat.gen(Zi, y, X, D, n)
	temp <- nleqslv(rep(0, ncol(obj@X) - 1), fn=U, jac=DU)
	if(verbose) {
		cat(sprintf("Check if gamma is solved correctly: %s \n", paste(U(temp$x), collapse=",")))
		cat(sprintf("message of nleqslv: %s ", temp$message))
	}
	if (sum(abs(U(temp$x))) > tol) stop("Failed to converge during solving gamma")
	return(temp$x)
}

inverse <- function(x) x[length(x):1]

#'@title Brewslow estimation
#'
#'@description An implementation of brewslow estimation, the formula (4) in 
#'Huang 2004.
#'
#'@param obj recurrent-data object
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
  g <- function(s) {
    index <- which(y >= s)
    1/sum(S[index])
  }
  f <- function(s) {
    index <- which(y <= s)
    sum(D[index] * sapply(y[index], g))
  }
  x <- inverse(unique(y))
  stepfun(x, append(0, sapply(x, f)))
}

phi_3.gen <- function(obj) {
  y <- obj@y
  m <- sapply(obj@t, length)
  X <- obj@X[,-1]
  F.hat.y <- obj$F.hat(y)
  term_1 <-  m / F.hat.y
  term_1[F.hat.y == 0] <- 0
  alpha <- BSM(obj)
  b.i <- lapply(1:length(y), b.hat.gen(obj))
  fi.seq <- fi.hat(obj)[-1,]
  Z_i <- Z_i.hat(obj)
  retval <- function(i, t, b) {
    term_2 <- term_1 * exp(X %*% (b - alpha)) * ifelse(obj@y > t, 1, 0)
    retval <- mean(term_2 * (X %*% fi.seq[,i] + Vectorize(b.i[[i]])(y)))
    retval + term_2[i] - mean(Z_i * exp(X %*% b) * ifelse(obj@y > t, 1, 0))
  }
}

phi_4.gen <- function(obj) {
  y <- obj@y
  m <- sapply(obj@t, length)
  X <- obj@X[,-1]
  F.hat.y <- obj$F.hat(y)
  term_1 <-  m / F.hat.y
  term_1[F.hat.y == 0] <- 0
  alpha <- BSM(obj)
  b.i <- lapply(1:length(y), b.hat.gen(obj))
  fi.seq <- fi.hat(obj)[-1,]
  Z_i <- Z_i.hat(obj)
  function(i, t, b) {
    term_2 <- term_1 * exp(X %*% (b - alpha)) * ifelse(obj@y > t, 1, 0)
    coef <- as.vector(term_2 * (X %*% fi.seq[,i] + Vectorize(b.i[[i]])(y)))
    retval <- (coef %*% X)/nrow(X)
    retval + term_2[i] * X[i,] - as.vector(Z_i * exp(X %*% b) * ifelse(obj@y > t, 1, 0)) %*% X / nrow(X)
  }
}

phi.gen <- function(obj) {
  y <- obj@y
  m <- sapply(obj@t, length)
  X <- obj@X[,-1]
  F.hat.y <- obj$F.hat(y)
  term_1 <-  m / F.hat.y
  alpha <- BSM(obj)
  b.i <- lapply(1:length(y), b.hat.gen(obj))
  fi.seq <- fi.hat(obj)[-1,]
  Z_i <- Z_i.hat(obj)
  phi_3 <- phi_3.gen(obj)
  phi_4 <- phi_4.gen(obj)
  D <- Delta_i(obj)
  function(i, beta) {
    term1 <- sapply(y, function(s) {
      term.0 <- as.vector(Z_i * exp(X %*% beta) * ifelse(y >= s, 1, 0))
      term.den <- mean(term.0)^2
      term.num <- term.0 %*% X / nrow(X)
      retval <- as.vector(phi_3(i, s, beta) * term.num / term.den)
    })
    term2 <- sapply(y, function(s) {
      term.0 <- as.vector(Z_i * exp(X %*% beta) * ifelse(y >= s, 1, 0))
      term.den <- mean(term.0)
      retval <- as.vector(phi_4(i, s, beta) / term.den)
    })
    term3 <- sapply(y, function(s) {
      term.0 <- as.vector(Z_i * exp(X %*% beta) * ifelse(y >= s, 1, 0))
      term.den <- mean(term.0)
      term.num <- term.0 %*% X / nrow(X)
      retval <- as.vector(term.num / term.den)
    })
    retval <- as.vector(D[i] * X[i,]) - as.vector(D %*% X / nrow(X))
    retval <- retval + as.vector(term1 %*% D) / nrow(X)
    retval <- retval - as.vector(term2 %*% D) / nrow(X)
    retval <- retval + as.vector(term3 %*% D) / nrow(X)
    a <- rep(0, nrow(X)) ; a[i] <- D[i]
    retval <- retval - as.vector(term3 %*% a)
    retval
  }  
}


Sigma.hat.gen <- function(obj) {
  n <- length(obj@y)
  phi <- phi.gen(obj)
  function(beta) {
    phi_i <- sapply(1:n, function(i) phi(i, beta))
    phi_star <- apply(phi_i, 1, mean)
    term <- phi_i - phi_star
    term %*% t(term) / n
  }
}

beta.hat <- function(obj) {
  
}

phi_i.gen <- function(obj) {
  
}

H0.hat.var <- function(obj) {
  
}