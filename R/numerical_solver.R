#'@title Solving Estimation Equation
#'
#'@description Let g(gamma) = X^t (y - exp(X gamma)) where X is n x p, y is n x 1, gamma is p x 1.
#'This function applies multivariate newton's method to solve gamma given X and y.
#'
#'@param X, numeric matrix. A n x p matrix.
#'@param y, numeric vector. A n x 1 vector.
#'@export
eq_solver <- function(X, y, x0 = NULL, steptol = 1e-06, iterlim = 100) {
	g <- function(gamma) t(X) %*% (y - exp(X %*% gamma))
	Jg <- function(gamma) -t(X) %*% diag(exp(c(X %*% gamma))) %*% X
	if (is.null(x0)) {
		x0 <- rep(0, ncol(X))
	}
	x <- solve(Jg(x0), Jg(x0) %*% x0 - g(x0))
	iter <- 1
	while (sum(abs(x - x0)) > steptol && iter < iterlim) {
		x0 <- x
		x <- solve(Jg(x0), Jg(x0) %*% x0 - g(x0))
		iter <- iter + 1
	}
	list(g = g(x), x = x, iter=iter, tol=sum(abs(x - x0)))
}
