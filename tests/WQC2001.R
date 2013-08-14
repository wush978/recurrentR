library(recurrentR)
set.seed(1)
# lambda <- function(x) exp(sin(x))
lambda <- function(x) exp(-x/10)
T_0 <- rpois(1, 40)

curve(lambda, from=0, to = T_0)

y <- rpois(nrow(iris), T_0)
y <- as.numeric(ifelse(y < T_0, y, T_0))
t <- sapply(y, function(y) {
	# 	browser()
	lambda_i <- function(t) exp(rnorm(1)) * lambda(t)
	retval <- gen_inhomo_poisson(lambda_i, y - 1, lambda_i(0))
	if (is.null(retval)) return(vector("numeric", 0))
	return(retval)
	# 	return(as.numeric(ceiling(retval)))
})
obj <- new("recurrent-data", model.matrix(~.,iris), y, t, data.frame(), T_0)
x.eval <- seq(0, obj@T_0, length=100)
y.eval <- obj$Lambda.hat(x.eval)
system.time({
	v <- a.var.hat(obj)
	v.eval <- sapply(x.eval, function(x) v(x))
})
u.eval <- y.eval + 2*sqrt(v.eval)
l.eval <- y.eval - 2*sqrt(v.eval)
plot(x.eval, y.eval, type="l", xlab="x", ylab="Lambda.hat", ylim=c(min(l.eval), max(u.eval)))
lines(x.eval, l.eval, lty=2)
lines(x.eval, u.eval, lty=2)
Lambda.hat <- obj$Lambda.hat
system.time({
	b.eval <- Lambda.hat(x.eval, bootstrap=TRUE)
})
u.eval <- y.eval + 2 * b.eval$error.measurement
l.eval <- y.eval - 2 * b.eval$error.measurement
lines(x.eval, u.eval, col=2, lty=2)
lines(x.eval, l.eval, col=2, lty=2)
legend("bottomright", c("Lambda.hat", "asymptotic pointwise C.I.", "bootstrap pointwise C.I."), lty=c(1, 2, 2), col=c(1, 1, 2))

