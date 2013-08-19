p <- rexp(100)
p <- p/sum(p)

n <- 10^5
g <- function(x) sin(x)
f <- stepfun(1:100, c(0, cumsum(p)))

stopifnot(all.equal(sum(sin(1:100) * p), step_integrate(g, f, 0, 100)))
plot(density(sapply(1:1000, function(i) mean(g(sample(1:100, n, TRUE, prob=p))))))
abline(v=step_integrate(g, f, 0, 100))