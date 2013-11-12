library(recurrentR)
library(microbenchmark)

obj <- local({
  n <- 100
  x <- sample(0:1, n, TRUE)
  T_0 <- 40
  y <- rexp(n, 1/T_0)
  D <- y < T_0
  y <- ifelse(D, y, T_0)
  t <- lapply(y, function(y) {
    n <- rpois(1, y * 0.4)
    if (n == 0) return(numeric(0))
    sort(runif(n, 0, y))
  })
  create_recurrent_data.numeric(y, D, t, T_0, matrix(x, ncol = 1))
})

Lambda_0.hat <- recurrentR:::Lambda_0.hat.gen(obj)
Lambda_0.hat.y <- Lambda_0.hat(obj@y)
gamma.hat <- recurrentR:::gamma.hat.gen(obj)
