library(recurrentR)

x <- 1:10
y <- rnorm(11)

f1 <- new(StepFunction, x, y)
f2 <- stepfun(x, y)

for(i in 0:11) {
	stopifnot(f1$call(i) == f2(i))
	stopifnot(f1$call(i + 0.5) == f2(i + 0.5))
}
