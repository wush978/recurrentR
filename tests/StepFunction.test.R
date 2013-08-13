library(recurrentR)
set.seed(1)
	
x <- 1:10
y <- rnorm(11)

f1 <- new(recurrentR:::StepFunction, x, y)
f2 <- stepfun(x, y)

for(i in 0:11) {
	stopifnot(f1$call(i) == f2(i))
	stopifnot(f1$call(i + 0.5) == f2(i + 0.5))
}


x <- 0:11
stopifnot(f1$sort_call( x ) == f2(x))
f1$sort_call(x)^2
f1$sort_call(x)/2
(f1^2)$sort_call(x)
stopifnot((f1^2)$sort_call( x ) == f1$sort_call( x )^2)