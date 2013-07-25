gen_bootstrap <- function(obj) {
	n <- length(obj@y)
	index_set <- 1:n
	return(function(B) {
		retval <- vector("list", B)
		for(i in seq_along(retval)) {
			index <- sample(index_set, n, TRUE)
			X <- ifelse(nrow(obj@X) > 0, obj@X[index,], obj@X)
			y <- obj@y[index]
			t <- obj@t[index]
			W <- ifelse(nrow(obj@W) > 0, obj@W[index,], obj@W)
			T_0 <- obj@T_0
			retval[[i]] <- new("recurrent-data", X, y, t, W, T_0)
		}
		return(retval)
	})
}