
#'@title Recurrent Data
#'
#'@description The S4 class for recurrent data.
#'
#'@param X, data.frame. The time independent covariate of subjects.
#'@param y, numeric vector. The censor time.
#'@param t, list. The time of recurrent events.
#'@param W, data.frame. The time dependent covariates.
#'@param T_0, numeric value. The time of termination of experients.
#'
#'@exportClass recurrent-data
setClass(
	"recurrent-data",
	representation(
		X = "matrix",
		y = "numeric",
		t = "list",
		W = "data.frame",
		T_0 = "numeric"
	)
)

setMethod(
	"initialize", 
	"recurrent-data",
	function(.Object, X, y, t, W, T_0) {
		.Object@X <- X
		.Object@y <- y
		.Object@t <- t
		.Object@W <- W
		.Object@T_0 <- T_0
		.Object
	}
	)


setMethod(
	"$",
	signature(x = "recurrent-data"),
	function (x, name) {
		switch(name,
			"F.hat" = F.hat(x),
			"Lambda.hat" = Lambda.hat(x),
			"bootstrap" = gen_bootstrap(x)
		)
	}
)
