library(recurrentR)
library(RPPGen)

lambda <- function(x) exp(-x/10)
par(mfrow=c(1, 2))
curve(lambda, 0, 40)
Lambda <- function(x) 10 * (1 - exp(-x/10))
curve(Lambda, 0, 40)
T_0 <- rpois(1, 40)
h <- function(t) rep(1/T_0, length(t))
gen_z <- function() runif(1, 0.5, 1.5)
n <- 150
beta <- c(1, -1)
gamma <- c(2, -2)
X <- cbind(sin(1:n), sample(c(0, 1), n, TRUE))

y <- rpois(n, T_0)
y <- as.numeric(ifelse(y < T_0, y, T_0))
t <- vector("list", n)
D <- numeric(n)
for(i in seq_along(y)) {
	z <- gen_z() 
	lambda_i <- function(t) z * lambda(t) * exp(as.vector(X[i,] %*% beta))
# 	h_i <- function(t) z * h(t) * exp(as.vector(X[i,] %*% gamma))
	D[i] <- rexp(1, z * exp(as.vector(X[i,]%*%gamma)) / T_0)
	t[[i]] <- gen_inhomo_poisson(lambda_i, y[i], lambda_i(0))
}
obj <- new("recurrent-data", X, y, t, data.frame(), T_0, D)

recurrentR:::H0.hat(obj, tol=1e-07, TRUE)