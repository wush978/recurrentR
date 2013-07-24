#'@title Generate Homogeneous Poisson Process
#'
#'@param lambda, numeric value.
#'@param T_0, positive numeric value.
#'
#'@description Generate a homogeneous poisson process with constant intensity \code{lambda} 
#'in time interval \code{[0,T_0]}.
#'
#'@export
gen_homo_poisson <- function(lambda, T_0) {
	result <- c()
	t <- 0
	while(t <= T_0) {
		u <- runif(1)
		t <- t - log(u) / lambda
		if (t <= T_0) result <- append(result, t)
	}
	result
}


#'@title Generate Inhomogeneous Poisson Process
#'
#'@param lambda, a positive function.
#'@param T_0, positive numeric value.
#'@param lambda_u, numeric value. The upper bound of \code{lambda(t)} in \code{[0, T_0]}.
#'If \code{lambda_u} is \code{NULL}, then \code{\link{optimize}} will be used to find the upper bound.
#'
#'@description Generate an inhomogeneous poisson process with intensity function \code{lambda(t)} 
#'in time interval \code{[0,T_0]}.
#'
#'@export
gen_inhomo_poisson <- function(lambda, T_0, lambda_u = NULL) {
	result <- c()
	lambda_u <- optimize(lambda, c(0, T_0), maximum = TRUE)$objective
	t <- 0
	while(t <= T_0) {
		u <- runif(1)
		t <- t - log(u) / lambda_u
		if (t > T_0) break
		u2 <- runif(1)
		if (u2 <= lambda(t)/lambda_u) result <- append(result, t)
	}
	result
}
