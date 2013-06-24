setClass(
  "recurrent-data", 
  representation=
    representation(
      y = "list",
      X = "data.frame",
      attr = "list"
      )
  )

setMethod("initialize",
          signature(.Object = "recurrent-data"),
          function (.Object, y, X) 
          {
            .Object@y <- y
            .Object@X <- X
            attr <- list()
            attr$s <- sort(unlist(y))
            attr$d <- table(attr$s)
            attr$s <- unique(attr$s)
            attr$N <- cumsum(attr$d[length(attr$d):1])[length(attr$d):1]
            attr$N <- attr$N 
            .Object@attr <- attr
            return(.Object)
          }
)
