

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
	if (is.null(lambda_u)) {
		lambda_u <- optimize(lambda, c(0, T_0), maximum = TRUE)$objective
	}
	result <- gen_homo_poisson(lambda_u, T_0)
	result[runif(length(result)) < sapply(result, function(r) lambda(r) / lambda_u)]
}
