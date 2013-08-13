\name{eq_solver}
\alias{eq_solver}
\title{Solving Estimation Equation}
\usage{
  eq_solver(X, y, x0 = NULL, steptol = 1e-06,
    iterlim = 100)
}
\arguments{
  \item{X,}{numeric matrix. A n x p matrix.}

  \item{y,}{numeric vector. A n x 1 vector.}
}
\description{
  Let g(gamma) = X^t (y - exp(X gamma)) where X is n x p, y
  is n x 1, gamma is p x 1. This function applies
  multivariate newton's method to solve gamma given X and
  y.
}

