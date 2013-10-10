\name{Z_i.hat}
\alias{Z_i.hat}
\title{Estimator of random effect \eqn{Z_i}}
\usage{
  Z_i.hat(obj, F.hat = NULL, gamma = NULL)
}
\arguments{
  \item{obj}{recurrent-data object}

  \item{F.hat}{cache of \code{obj$F.hat}}

  \item{F.hat}{cache of \code{obj$U.hat()}.}
}
\description{
  \deqn{\hat{Z}_i = \frac{m_i}{\hat{\Lambda}_0(Y_i) e^{X_i
  \hat{\alpha}}}}
}

