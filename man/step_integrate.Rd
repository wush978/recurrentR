\name{step_integrate}
\alias{step_integrate}
\title{Integration on step function}
\usage{
  step_integrate(g, f, a, b)
}
\arguments{
  \item{g}{Integrated function of the integral}

  \item{f}{Step function}

  \item{a}{Numeric value}

  \item{b}{Numeric value}
}
\value{
  Numeric value
}
\description{
  Evaluate \deqn{\int_a^b {g(x) df(x)}} where \eqn{f} is a
  step function
}

