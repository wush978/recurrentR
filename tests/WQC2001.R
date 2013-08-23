library(recurrentR)
# set.seed(1)
# lambda <- function(x) exp(sin(x))

do_exp <- function(lambda, Lambda, T_0, gen_z) {
	curve(lambda, from=0, to = T_0)
	
	y <- rpois(nrow(iris), T_0)
	y <- as.numeric(ifelse(y < T_0, y, T_0))
	t <- sapply(y, function(y) {
		# 	browser()
		z <- gen_z()
		lambda_i <- function(t) z * lambda(t)
		retval <- gen_inhomo_poisson(lambda_i, y - 1, lambda_i(0))
		if (is.null(retval)) return(vector("numeric", 0))
		return(retval)
		# 	return(as.numeric(ceiling(retval)))
	})
	obj <- new("recurrent-data", model.matrix(~.,iris), y, t, data.frame(), T_0)
	
	c.hat.gen <- recurrentR:::c.hat.gen(obj)
	ci <- sapply(1:length(t), c.hat.gen)
	mean(ci)
	mean(ci^2)
	Lambda(T_0)
	obj$Lambda.hat(T_0)
	
	x.eval <- seq(0, obj@T_0, length=100)
	y.eval <- obj$Lambda.hat(x.eval)
	system.time({
		v <- a.var.hat(obj)
		v.eval <- sapply(x.eval, function(x) v(x))
	})
	u.eval <- y.eval + 2*sqrt(v.eval)
	l.eval <- y.eval - 2*sqrt(v.eval)
	plot(x.eval, Lambda(x.eval), type="l", xlab="x", ylab="", 
		col = 3, lwd = 2, ylim=c(min(l.eval), max(u.eval)))
	lines(x.eval, y.eval, lty=1, col = 1)
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
	
	m <- 100
	real.sim <- matrix(0, length(x.eval), m)
	pb <- txtProgressBar(style=3, max=m)
	for(i in 1:m) {
		t.sim <- sapply(y, function(y) {
			# 	browser()
			z <- gen_z()
			lambda_i <- function(t) z * lambda(t)
			retval <- gen_inhomo_poisson(lambda_i, y - 1, lambda_i(0))
			if (is.null(retval)) return(vector("numeric", 0))
			return(retval)
			# 	return(as.numeric(ceiling(retval)))
		})
		obj.sim <- new("recurrent-data", model.matrix(~.,iris), y, t.sim, data.frame(), T_0)
		real.sim[,i] <- obj.sim$Lambda.hat(x.eval)
		setTxtProgressBar(pb, i)
	}
	close(pb)
	
	lines(x.eval, apply(real.sim, 1, function(a) quantile(a, 0.025)), lty=2, col=3, lwd=2)
	lines(x.eval, apply(real.sim, 1, function(a) quantile(a, 0.975)), lty=2, col=3, lwd=2)
	
	legend("bottomright", 
		c("Lambda.hat", "Lambda", "asymptotic pointwise C.I.", "bootstrap pointwise C.I.", "parametric bootstrap pointwise C.I."), 
		lty=c(1, 1, 2, 2, 2), col=c(1, 3, 1, 2, 3), lwd=c(1, 2, 1, 1, 2))
}

par(mfrow=c(2, 2))
lambda <- function(x) exp(-x/10)
Lambda <- function(x) 10 * (1 - exp(-x/10))
T_0 <- rpois(1, 40)
gen_z <- function() runif(1, 0.5, 1.5)
do_exp(lambda, Lambda, T_0, gen_z)

lambda <- function(x) rep(1, length(x))
Lambda <- function(x) x
T_0 <- rpois(1, 40)
gen_z <- function() rexp(1)
do_exp(lambda, Lambda, T_0, gen_z)
