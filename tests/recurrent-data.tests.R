library(recurrentR)

set.seed(1)

# Data Generation

lambda <- function(x) exp(sin(x))
T_0 <- rpois(1, 40)

curve(lambda, from=0, to = T_0)

y <- rpois(nrow(iris), T_0)
y <- as.numeric(ifelse(y < T_0, y, T_0))
t <- sapply(y, function(y) {
	as.numeric(ceiling(gen_inhomo_poisson(lambda, y - 1, 4)))
})
obj <- new("recurrent-data", iris, y, t, data.frame(), T_0)
rm(list=c("y", "t", "T_0"), envir=globalenv())

# Estimate WQC2001 Model 1(No covariates)

F.hat <- obj$F.hat
Lambda.hat <- obj$Lambda.hat
curve(F.hat, 0, obj@T_0)
curve(Lambda.hat, 0, obj@T_0)

curve(lambda, 0, obj@T_0)
Lambda.single <- function(t) integrate(lambda, 0, t)$value
Lambda <- function(T_0) {
	return(function(t) sapply(t, Lambda.single))
}
answer <- Lambda(obj@T_0)
curve(answer, 0, obj@T_0, col=2)
curve(Lambda.hat, 0, obj@T_0, add=TRUE, lty=1)

# Bootstrap Estimate

x.eval <- seq(from=0, to=obj@T_0, length=100)
obj$F.hat(x.eval, TRUE, error.measurement.function=sd)
obj$F.hat(x.eval, TRUE, error.measurement.function=function(a) c(quantile(a, 0.975), quantile(a, 0.025)))

result <- obj$Lambda.hat(x.eval, TRUE, error.measurement.function=sd)
curve(answer, 0,obj@T_0, col=2)
lines(x.eval, result$estimate, col=1)
lines(x.eval, result$estimate + 2 * result$error.measurement, lty = 2, col = 1)
lines(x.eval, result$estimate - 2 * result$error.measurement, lty = 2, col = 1)

result <- obj$Lambda.hat(x.eval, TRUE, error.measurement.function=function(a) c(quantile(a, 0.975), quantile(a, 0.025)))
curve(answer, 0,obj@T_0, col=2)
lines(x.eval, result$estimate, col=1)
lines(x.eval, result$error.measurement[1,], lty = 2, col = 1)
lines(x.eval, result$error.measurement[2,], lty = 2, col = 1)
