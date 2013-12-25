Simulation Study
===


```r
library(recurrentR)
```

```
## Loading required package: Rcpp Loading required package: nleqslv
```

```r
library(RPPGen)
library(plotrix)
```


# [@Wang2001]

## Parameters and Generators


```r
a <- 0.5
b <- 0.3
p <- 0.5
frailty.dist = c("Pois", "Gamma")
T_0 <- 10
n <- 400
lambda_0 <- local({
    lambda_0 <- function(u) rep(1/10, length(u))
    Total <- T_0/10
    function(u) {
        lambda_0(u)/Total
    }
})
Lambda_0 <- function(u) u/10
lambda_y <- sqrt(200)
k_y <- 2
h_0 <- function(u) k_y * u^(k_y - 1)/lambda_y^k_y
H_0 <- function(t) (t/lambda_y)^k_y

x <- rbinom(n, 1, 0.5)
Ez <- switch(frailty.dist[1], Pois = 10, Gamma = 2/5)
alpha.true <- b
gamma.true <- c(log(Ez), a)
gen_huang_2004 <- function() {
    x <- rbinom(n, 1, p)
    switch(frailty.dist[1], Pois = {
        z <- rpois(n, 10)
    }, Gamma = {
        z <- rgamma(n, 2, scale = 5)
    })
    y <- numeric(n)
    t <- vector("list", n)
    for (i in seq_len(n)) {
        y[i] <- rweibull(1, shape = k_y, scale = lambda_y * (z[i] * exp(b * 
            x[i]))^(-1/k_y))
        y[i] <- min(T_0, y[i])
        lambda_i <- function(t) z[i] * exp(a * x[i]) * lambda_0(t)
        t[[i]] <- gen_inhomo_poisson(lambda_i, y[i], z[i] * exp(a * x[i]))
    }
    create_recurrent_data.numeric(y, y < T_0, t, T_0, matrix(x, nrow = n))
}
```


## Asymptotic Variance


```r
obj <- gen_huang_2004()
wang_2001 <- Wang2001(obj)
huang_2004 <- Huang2004(obj = obj)
curve(Lambda_0, 0, 10, col = 2, ylim = c(0, 1.2))
Lambda_0.hat <- function(x) wang_2001$Lambda_0.hat(x)
curve(Lambda_0.hat, add = TRUE, 0, 10, col = 1)
Lambda_0.hat.upper <- Vectorize(function(x) wang_2001$Lambda_0.hat(x) + 2 * 
    sqrt(wang_2001$Lambda_0.hat.var(x)))
curve(Lambda_0.hat.upper, add = TRUE, col = 1, 0, 10, lty = 2)
Lambda_0.hat.lower <- Vectorize(function(x) wang_2001$Lambda_0.hat(x) - 2 * 
    sqrt(wang_2001$Lambda_0.hat.var(x)))
curve(Lambda_0.hat.lower, add = TRUE, col = 1, 0, 10, lty = 2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r

curve(H_0, 0, T_0, col = 2, ylim = c(0, H_0(T_0) * 1.1))
lines(huang_2004$y, huang_2004$H0.hat.y)
lines(huang_2004$y, huang_2004$H0.hat.y + qnorm(0.975) * huang_2004$H0.hat.y.sd, 
    lty = 2)
lines(huang_2004$y, huang_2004$H0.hat.y - qnorm(0.975) * huang_2004$H0.hat.y.sd, 
    lty = 2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 

```r

plotCI(1:3, c(wang_2001$gamma.bar.hat, huang_2004$alpha.hat), uiw = qnorm(0.975) * 
    c(sqrt(diag(wang_2001$gamma.bar.hat.var)), sqrt(diag(huang_2004$alpha.hat.var))))
points(1:3, c(gamma.true, alpha.true), col = 2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-33.png) 


## Parametric Bootstrap


```r
B <- 20
gamma.hat <- matrix(NA, B, 2)
x.eval <- seq(from = 0, to = T_0, length.out = 100)
Lambda_0.hat <- matrix(NA, B, length(x.eval))
if (interactive()) pb <- txtProgressBar(max = B, style = 3)
for (i in seq_len(B)) {
    temp <- Wang2001(obj = gen_huang_2004())
    gamma.hat[i, ] <- temp$gamma.bar.hat
    Lambda_0.hat[i, ] <- temp$Lambda_0.hat(x.eval)
    gc()
    if (interactive()) 
        setTxtProgressBar(pb, i)
}
if (interactive()) close(pb)
curve(Lambda_0, 0, 10, col = 2, ylim = c(0, 1.2))
y.eval <- apply(Lambda_0.hat, 2, mean)
y.eval.sd <- apply(Lambda_0.hat, 2, sd)
lines(x.eval, y.eval)
lines(x.eval, y.eval + qnorm(0.975) * y.eval.sd, lty = 2)
lines(x.eval, y.eval - qnorm(0.975) * y.eval.sd, lty = 2)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-41.png) 

```r
plotCI(1:2, apply(gamma.hat, 2, mean), uiw = apply(gamma.hat, 2, sd), main = "Parametric Bootstrap", 
    xlab = "", ylab = "", xaxt = "n")
points(1:2, gamma.true, col = 2)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-42.png) 

