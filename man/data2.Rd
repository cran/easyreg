\name{data2}
\alias{data2}
\docType{data}
\title{
data2: Kaps and Lamberson (2009): page 434
}
\description{
Quantitative treatments in randomizad block design.
}
\usage{data(data2)}
\format{
  A data frame with 25 observations on the following 3 variables.
  \describe{
    \item{\code{protein_level}}{a numeric vector}
    \item{\code{litter}}{a factor with levels \code{l1} \code{l2} \code{l3} \code{l4} \code{l5}}
    \item{\code{feed_conversion}}{a numeric vector}
  }
}
\references{
KAPS, M. and LAMBERSON, W. R. Biostatistics for Animal Science: an introductory text. 2nd Edition. CABI Publishing, Wallingford, Oxfordshire, UK, 2009. 504p.
}
\examples{
data(data2)
summary(data2)
}

