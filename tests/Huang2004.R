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
gamma <- c(2, -3)
X <- cbind(sin(1:n), sample(c(0, 1), n, TRUE))

y <- rpois(n, T_0)
y <- as.numeric(ifelse(y < T_0, y, T_0))
t <- vector("list", n)
D <- numeric(n)
Z.true <- sapply(1:n, function(i) gen_z())
for(i in seq_along(y)) {
	z <- Z.true[i]
	lambda_i <- function(t) z * lambda(t) * exp(as.vector(X[i,] %*% beta))
# 	h_i <- function(t) z * h(t) * exp(as.vector(X[i,] %*% gamma))
  h_i <- 1/T_0 * z * exp(as.vector(X[i,]%*%gamma))
	D[i] <- rexp(1,  h_i)
	t[[i]] <- gen_inhomo_poisson(lambda_i, min(D[i], y[i]), lambda_i(0))
}
stopifnot(all(sapply(t, function(v) ifelse(length(v) > 0, max(v), 0)) < D))
stopifnot(all(sapply(t, function(v) ifelse(length(v) > 0, max(v), 0)) < y))
D.index <- D < y
y[D < y] <- D[D < y]
obj <- new("recurrent-data", X, y, t, data.frame(), T_0, D.index)

# validate in WQC2001

function() {
  
  c.hat.gen <- recurrentR:::c.hat.gen(obj)
  ci <- sapply(1:length(t), c.hat.gen)
  # 	obj$Lambda.hat(T_0)
  
  x.eval <- seq(0, obj@T_0, length=100)
  y.eval <- obj$Lambda.hat(x.eval)
  gamma <- obj$U.hat()
  system.time({
  	v <- asymptotic.var(obj, gamma = gamma)
  	v.eval <- sapply(x.eval, function(x) v$Lambda.hat.var(x))
  })
  u.eval <- y.eval + 2*sqrt(v.eval)
  l.eval <- y.eval - 2*sqrt(v.eval)
  plot(x.eval, Lambda(x.eval), type="l", xlab="x", ylab="", 
  		 col = 3, lwd = 2, ylim=c(min(l.eval), max(u.eval)))
  lines(x.eval, y.eval, lty=1, col = 1)
  lines(x.eval, l.eval, lty=2)
  lines(x.eval, u.eval, lty=2)
}

# F.hat(obj@y)
# 
# lm(y ~ x -1, )
data.frame(y = Z.true, x = recurrentR:::Z_i.hat(obj))
recurrentR:::BSM(obj)
H0 <- recurrentR:::H0.hat(obj, tol=1e-07, TRUE)
curve(H0, 0, obj@T_0)
H <- function(t) t / T_0
curve(H, 0, obj@T_0, add=TRUE, col=3)