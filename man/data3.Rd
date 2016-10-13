\name{data3}
\alias{data3}
\docType{data}
\title{
data3: fictional example
}
\description{
Quantitative treatments in latin square design.
}
\usage{data(data3)}
\format{
  A data frame with 25 observations on the following 4 variables.
  \describe{
    \item{\code{treatment}}{a numeric vector}
    \item{\code{animal}}{a factor with levels \code{a1} \code{a2} \code{a3} \code{a4} \code{a5}}
    \item{\code{period}}{a factor with levels \code{p1} \code{p2} \code{p3} \code{p4} \code{p5}}
    \item{\code{milk_fat}}{a numeric vector}
  }
}

\examples{
data(data3)
summary(data3)
}

