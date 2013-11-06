\name{create_recurrent_data}
\alias{create_recurrent_data}
\title{Creator of recurrent-data}
\usage{
  create_recurrent_data(X, y, t, W, T_0, D)
}
\arguments{
  \item{X}{data.frame. The time independent covariate of
  subjects.}

  \item{y}{numeric vector. The censor time.}

  \item{t}{list. The time of recurrent events.}

  \item{W}{data.frame. The time dependent covariates.}

  \item{T_0}{numeric value. The time of termination of
  experients.}

  \item{D}{numeric vector. The failure time.}
}
\description{
  Generate the recurrent-data object.
}

