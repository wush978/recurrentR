\name{BorrowStrengthMethod}
\alias{BorrowStrengthMethod}
\title{Borrow Strength Method}
\usage{
  BorrowStrengthMethod(obj, U.hat = NULL, Gamma.hat = NULL,
    Zi = NULL, F.hat = NULL, gamma.hat = NULL,
    verbose = FALSE, tol = 1e-04, ...)
}
\arguments{
  \item{obj}{recurrent-data object}

  \item{U.hat}{function, \eqn{\hat{U}(\beta)} in Huang2004}

  \item{Gamma.hat}{function, \eqn{\hat{Gamma}(\beta) :=
  \frac{d\hat{U}}{d\beta}} in Huang2004}

  \item{Zi}{numeric vector, \eqn{\hat{Z}_i} in Huang2004}

  \item{F.hat}{function, \eqn{\hat{\Lambda_0}} in
  Huang2004}

  \item{gamma.hat}{numeric vecotr, \eqn{\hat{alpha}} in
  Huang2004 and \eqn{\hat{\gamma}} in WQC2001}

  \item{verbose}{boolean value, whether print the message
  of \code{nleqslv} or not. Please see
  \code{\link{nleqslv}}.}

  \item{tol}{numeric value, used to determine convergence.
  Please see \code{\link{nleqslv}}.}

  \item{...}{}
}
\description{
  Borrow Strength Method
}

