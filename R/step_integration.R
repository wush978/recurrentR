#'@title Integration on step function
#'
#'@param g Integrated function of the integral
#'@param f Step function
#'@param a Numeric value
#'@param b Numeric value
#'
#'@description 
#'Evaluate \deqn{\int_a^b {g(x) df(x)}} where \eqn{f} is a step function
#'
#'@return Numeric value
#'
#'@export
step_integrate <- function(g, f, a, b) {
	stopifnot("stepfun" %in% class(f))
	e <- environment(f)
	x <- e[["x"]]
	y <- c(e[["f"]], e[["y"]])
	p <- diff(y)
	index <- which(x >= a & x <= b)
	sum(g(x[index]) * p[index])
}