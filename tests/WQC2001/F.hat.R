library(recurrentR)

F.hat.prototype <- function(obj) {
  s <- recurrentR:::s(obj)
  d <- recurrentR:::d(obj)
  N <- sapply(s, function(s) {
    retval <- 0
    for(i in seq_along(obj@y)) {
      if(s > obj@y[i]) next
      retval <- retval + sum(obj@t[[i]] <= s)
    }
    retval
  })
  function(t) {
    index <- which(s > t)
    if (length(index) == 1) return(1)
    prod(1 - d[index] / N[index])
  } 
}

data(obj.list)
for(i in seq_along(obj.list)) {
  obj <- obj.list[[i]]
  F.hat0 <- F.hat.prototype(obj)
  F.hat1 <- obj$F.hat
  F.hat0.y <- sapply(obj@y, function(y) F.hat0(y))
  F.hat1.y <- F.hat1(obj@y)
  stopifnot(F.hat0.y == F.hat1.y)
  which(F.hat0.y != F.hat1.y)
  F.hat0.y[81]
  F.hat1.y[81]
}
