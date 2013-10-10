\name{DSSi}
\alias{DSSi}
\title{Moment Estimators associated with subject-specific latent variable}
\usage{
  DSSi(Beta, Zi, X)
}
\value{
  numeric vector
}
\description{
  \eqn{\frac{d}{d\beta}SS_i = \sum_{j=1}^{i} S_j X_j^T} The
  \code{Zi} and \code{X} are subject-specific latent
  varialbe and time-independent covariate associated with
  ordered \code{y}
}
\seealso{
  \code{\link{SSi}}
}

