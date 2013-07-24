\name{gen_inhomo_poisson}
\alias{gen_inhomo_poisson}
\title{Generate Inhomogeneous Poisson Process}
\usage{
  gen_inhomo_poisson(lambda, T_0, lambda_u = NULL)
}
\arguments{
  \item{lambda,}{a positive function.}

  \item{T_0,}{positive numeric value.}

  \item{lambda_u,}{numeric value. The upper bound of
  \code{lambda(t)} in \code{[0, T_0]}. If \code{lambda_u}
  is \code{NULL}, then \code{\link{optimize}} will be used
  to find the upper bound.}
}
\description{
  Generate an inhomogeneous poisson process with intensity
  function \code{lambda(t)} in time interval
  \code{[0,T_0]}.
}

