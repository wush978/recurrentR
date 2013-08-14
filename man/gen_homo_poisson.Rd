\name{gen_homo_poisson}
\alias{gen_homo_poisson}
\title{Generate Homogeneous Poisson Process}
\usage{
  gen_homo_poisson(lambda, T_0)
}
\arguments{
  \item{lambda}{Numeric value, the intensity.}

  \item{T_0}{Positive numeric value, the upper bound.}
}
\description{
  Generate a homogeneous poisson process with constant
  intensity \code{lambda} in time interval \code{[0,T_0]}.
}

