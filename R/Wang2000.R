#'@export
F.hat <- function(obj, t) {
  index <- obj@attr$s > t
  d <- obj@attr$d[index]
  N <- obj@attr$N[index]
  prod(1 - d/N)
}