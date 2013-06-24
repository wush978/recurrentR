library(recurrentR)
y <- list(c(1,2,3), c(2,5,7))
obj <- new("recurrent-data", y, data.frame())
obj$d(1.0)

F.hat(obj, 7)