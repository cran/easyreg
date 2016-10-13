\name{data4}
\alias{data4}
\docType{data}
\title{
data4: fictional example
}
\description{
Quantitative treatments in several latin squares design.
}
\usage{data(data4)}
\format{
  A data frame with 50 observations on the following 5 variables.
  \describe{
    \item{\code{treatment}}{a numeric vector}
    \item{\code{square}}{a numeric vector}
    \item{\code{animal}}{a factor with levels \code{a1} \code{a2} \code{a3} \code{a4} \code{a5}}
    \item{\code{period}}{a factor with levels \code{p1} \code{p2} \code{p3} \code{p4} \code{p5}}
    \item{\code{milk_fat}}{a numeric vector}
  }
}

\examples{
data(data4)
summary(data4)
}

