library(recurrentR)
data(huang_2010)
obj.raw <- attr(huang_2010, "restore")(huang_2010[[2]]) 
obj <- create_recurrent_data(
  matrix(obj.raw$W, ncol=1), y=obj.raw$Y, t=obj.raw$t, 
  W=lapply(obj.raw$X, function(x) function(t) c(x * log(t), x * sin(t))), 10, obj.raw$Y < 10
)

rho.gen <- function(obj, i, k) {
  function(t, u) {
    obj@W[[i]](u) + obj@W[[k]](t) - obj@W[[i]](t) - obj@W[[k]](u)
  }
}

h <- function(obj, i, j, beta) {
  browser()
  upper_bound <- min(obj@y[i], obj@y[j])
  index_i <- which(obj@t[[i]] < upper_bound)
  index_j <- which(obj@t[[j]] < upper_bound)
  rho <- Vectorize(rho.gen(obj, i, j))
  retval <- 1
  for(i1 in index_i) {
    for(i2 in index_j) {
      retval <- retval * 1/(1 + exp(as.vector(rho(obj@t[[i]][i1], obj@t[[j]][i2])) %*% beta))
    }
  }
  retval
}

h(obj, 1, 2, 0.5)